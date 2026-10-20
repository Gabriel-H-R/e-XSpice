/**
 * @author Gabriel Henrique Silva
 * @file SmallSignalLinearization.cpp
 * @brief Operating point + (G, C) snapshot + ordering, shared by .AC and .NOISE.
 */

#include "SmallSignalLinearization.h"

#include "../core/Device.h"
#include "../core/Node.h"
#include "../parallel/MPIManager.h"
#include "DaeSystem.h"

#include <algorithm>
#include <iostream>

namespace e_XSpice {

bool SmallSignalLinearization::prepare(Circuit& circuit, const Options& opt) {
    ok_ = false;
    circuit_ = &circuit;

    system_.reset(new DaeSystem(circuit));
    DaeSystem& system = *system_;

    // Frequency distribution replicates the same (G, C) on every rank, so the
    // inner device/row distribution is disabled in that mode (the caller tells
    // us which one it is doing).
    if (opt.use_distributed_layout) {
        MPIManager& mpi = MPIManager::getInstance();
        system.setMPIManager(&mpi);
    }
    system.setup();
    system.setMaxIterations(opt.max_iterations);
    system.setRelTol(opt.rel_tol);
    system.setAbsTol(opt.abs_tol);
    system.setGminSchedule(1e-3, 1e-12);
    system.setSolverType(SolverType::BICGSTAB);
    system.setVerbose(opt.verbose);
    system.beginSolve(/*initial_solve=*/true);
    system.resetSolveStats();

    const NLSolveResult nl = system.solveNonlinear(/*time=*/0.0,
                                                   /*q_coeff=*/0.0,
                                                   /*use_history=*/false);
    if (!nl.converged) {
        message_ = "operating point did not converge: " + nl.message;
        return false;
    }
    system.extractSolution();
    newton_iterations_ = nl.total_iterations;

    // One refresh so (G, C) belong to the CONVERGED point (the Newton loop
    // returns with the Jacobians of the previous iterate), then the pair is
    // snapshotted: dF/dx and dQ/dx share the CSR pattern, so the frequency
    // sweep just rescales C by w.
    system.refreshLinearization(/*include_gmin=*/false);
    const DaeMatrix& G = system.jacobianF();
    const DaeMatrix& C = system.jacobianQ();

    n_full_ = system.size();
    num_node_rows_ = std::max(0, system.numNodeRows());

    // ---- ground elimination ------------------------------------------------
    gnd_row_ = -1;
    {
        const auto& nodes = circuit.getNodes();
        const auto& node_map = system.nodeMap();
        for (const auto& n : nodes) {
            if (!n->isGround()) continue;
            const int id = n->getId();
            if (id >= 0 && id < static_cast<int>(node_map.size())) {
                gnd_row_ = node_map[static_cast<size_t>(id)];
            }
            break;
        }
    }

    ac_of_row_.assign(static_cast<size_t>(n_full_), -1);
    row_of_ac_.clear();
    row_of_ac_.reserve(static_cast<size_t>(n_full_));
    for (int r = 0; r < n_full_; ++r) {
        if (r == gnd_row_) continue;
        ac_of_row_[static_cast<size_t>(r)] = static_cast<int>(row_of_ac_.size());
        row_of_ac_.push_back(r);
    }
    nn_ = static_cast<int>(row_of_ac_.size());
    if (nn_ <= 0) {
        message_ = "nothing to solve: the circuit has no unknowns";
        return false;
    }

    // Filtered CSR (same pattern, ground row/column dropped) plus the matching
    // value vectors: one O(nnz) pass, done once for the whole sweep.
    row_ptr_.assign(static_cast<size_t>(nn_) + 1, 0);
    col_idx_.clear();
    g_vals_.clear();
    c_vals_.clear();
    {
        const std::vector<int>& rp = G.rowPointers();
        const std::vector<int>& ci = G.colIndices();
        const std::vector<double>& gv = G.values();
        const std::vector<double>& cv = C.values();
        col_idx_.reserve(gv.size());
        g_vals_.reserve(gv.size());
        c_vals_.reserve(gv.size());
        for (int i = 0; i < nn_; ++i) {
            const int r = row_of_ac_[static_cast<size_t>(i)];
            row_ptr_[static_cast<size_t>(i)] = static_cast<int>(col_idx_.size());
            for (int p = rp[r]; p < rp[r + 1]; ++p) {
                const int c = ci[static_cast<size_t>(p)];
                if (c == gnd_row_) continue;      // -g·x(gnd) term disappears
                col_idx_.push_back(c);
                g_vals_.push_back(gv[static_cast<size_t>(p)]);
                c_vals_.push_back(cv[static_cast<size_t>(p)]);
            }
        }
        row_ptr_[static_cast<size_t>(nn_)] = static_cast<int>(col_idx_.size());
    }

    // ---- pivot-safe ordering ----------------------------------------------
    // The frequency system must have a usable diagonal in every column before a
    // no-pivot factorization can touch it. The ordering is decided once, from
    // the matrix at the FIRST frequency (the DAE pattern declares a diagonal
    // slot for every row, so the values are what tell a real self-term from an
    // empty one) and then frozen for the whole sweep. Branch rows are preferred
    // partners: their couplings are the ±1 of a branch equation, so they make
    // a well-scaled pivot.
    {
        std::vector<int> priority(static_cast<size_t>(nn_), 0);
        for (int k = 0; k < (n_full_ - num_node_rows_); ++k) {
            const int a = ac_of_row_[static_cast<size_t>(num_node_rows_ + k)];
            if (a >= 0) priority[static_cast<size_t>(a)] = 1;
        }
        // Any positive reference frequency works; the ordering decision only
        // looks at WHICH entries exist and how big they are relative to each
        // other, which does not change with w in a circuit matrix. A DC
        // analysis passes 0, because its matrix is G alone.
        const double kOmegaRef = opt.reference_omega;
        std::vector<double> ref_re(col_idx_.size());
        std::vector<double> ref_im(col_idx_.size());
        for (size_t p = 0; p < col_idx_.size(); ++p) {
            ref_re[p] = g_vals_[p];
            ref_im[p] = kOmegaRef * c_vals_[p];
        }
        perm_ = StaticRowPermutation::build(row_ptr_, col_idx_, ref_re, ref_im,
                                            priority);
        if (!perm_.ok) {
            message_ = "frequency system cannot be ordered for factorization: " +
                       perm_.reason;
            return false;
        }
        if (perm_.weak_diagonals > 0 && opt.verbose) {
            std::cout << "[SmallSignal] pivot-safe ordering: "
                      << perm_.weak_diagonals
                      << " column(s) without a self-term reordered ("
                      << perm_.swaps << " row swap(s))" << std::endl;
        }
    }

    // Apply the row permutation to the pattern AND the value vectors.
    {
        std::vector<int> new_rp(static_cast<size_t>(nn_) + 1, 0);
        std::vector<int> new_ci;
        std::vector<double> g2, c2;
        new_ci.reserve(col_idx_.size());
        g2.reserve(g_vals_.size());
        c2.reserve(c_vals_.size());
        for (int k = 0; k < nn_; ++k) {
            new_rp[static_cast<size_t>(k)] = static_cast<int>(g2.size());
            const int src = perm_.row_of_pos[static_cast<size_t>(k)];
            for (int p = row_ptr_[static_cast<size_t>(src)];
                 p < row_ptr_[static_cast<size_t>(src) + 1]; ++p) {
                new_ci.push_back(col_idx_[static_cast<size_t>(p)]);
                g2.push_back(g_vals_[static_cast<size_t>(p)]);
                c2.push_back(c_vals_[static_cast<size_t>(p)]);
            }
        }
        new_rp[static_cast<size_t>(nn_)] = static_cast<int>(g2.size());
        row_ptr_.swap(new_rp);
        col_idx_.swap(new_ci);
        g_vals_.swap(g2);
        c_vals_.swap(c2);
    }

    ok_ = true;
    return true;
}

std::vector<double> SmallSignalLinearization::permuteRhs(
    const std::vector<double>& full) const {
    std::vector<double> out(static_cast<size_t>(nn_), 0.0);
    for (int k = 0; k < nn_; ++k) {
        // position k now holds the equation of full row
        // row_of_ac[perm.row_of_pos[k]]: the RHS follows the rows.
        const int r = row_of_ac_[static_cast<size_t>(
            perm_.row_of_pos[static_cast<size_t>(k)])];
        out[static_cast<size_t>(k)] =
            (r >= 0 && r < static_cast<int>(full.size()))
                ? full[static_cast<size_t>(r)]
                : 0.0;
    }
    return out;
}

std::vector<double> SmallSignalLinearization::dropGroundRow(
    const std::vector<double>& full) const {
    // Keeps the ORIGINAL row order (only the ground row is removed): this
    // is the indexing a TRANSPOSED solve expects for its right-hand side.
    std::vector<double> out(static_cast<size_t>(nn_), 0.0);
    int k = 0;
    for (int r = 0; r < n_full_; ++r) {
        if (r == gnd_row_) continue;
        out[static_cast<size_t>(k++)] =
            (r < static_cast<int>(full.size())) ? full[static_cast<size_t>(r)] : 0.0;
    }
    return out;
}

std::vector<double> SmallSignalLinearization::unpermuteSolution(
    const std::vector<double>& perm) const {
    std::vector<double> out(static_cast<size_t>(n_full_), 0.0);
    unpermuteSolution(perm, out);
    return out;
}

void SmallSignalLinearization::unpermuteSolution(
    const std::vector<double>& perm, std::vector<double>& out) const {
    out.assign(static_cast<size_t>(n_full_), 0.0);
    for (int k = 0; k < nn_; ++k) {
        const int r = row_of_ac_[static_cast<size_t>(k)];
        if (r < 0) continue;
        out[static_cast<size_t>(r)] =
            (k < static_cast<int>(perm.size())) ? perm[static_cast<size_t>(k)] : 0.0;
    }
}

} // namespace e_XSpice
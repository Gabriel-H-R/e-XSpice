/**
 * @author Gabriel Henrique Silva
 * @file DaeSystem.cpp
 * @brief Shared DAE/NL solver core for all analysis engines.
 */

#include "DaeSystem.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <cstdio>
#include <iostream>

namespace e_XSpice {

DaeSystem::DaeSystem(Circuit& circuit)
    : circuit_(circuit) {}

// ===============================================================
// SETUP
// ===============================================================

void DaeSystem::setup() {
    // ---- 1) node map: ground is never materialized ----
    const auto& nodes = circuit_.getNodes();
    num_nodes_ = static_cast<int>(nodes.size());
    node_map_.assign(circuit_.getNumNodes(), -1);
    for (int i = 0; i < num_nodes_; ++i) {
        node_map_[nodes[i]->getId()] = i;
    }

    // ---- 2) extra unknowns (branch currents) ----
    devices_.assign(circuit_.getDevices().begin(), circuit_.getDevices().end());
    num_extra_ = 0;
    for (Device* d : devices_) {
        num_extra_ += std::max(0, d->numExtraUnknowns());
    }
    n_ = num_nodes_ + num_extra_;

    // ---- 3) topology to devices ----
    int extra = num_nodes_;
    for (Device* d : devices_) {
        if (d->numExtraUnknowns() > 0) {
            d->setBranchIndex(extra);
            extra += d->numExtraUnknowns();
        } else {
            d->setBranchIndex(-1);
        }
        d->setTopology(node_map_);
    }

    // ---- 4) per-device state history ----
    for (Device* d : devices_) {
        d->registerStates(store_);
    }

    // ---- 5) fixed sparsity pattern (union of the dFdx and dQdx stamps) ----
    // Both matrices share the SAME pattern, hence the SAME nnz ordering, so
    // resolving slots once yields positions valid in both matrices.
    std::vector<std::vector<int>> rows(devices_.size());
    std::vector<std::vector<int>> cols(devices_.size());

    jac_f_ = DaeMatrix(n_);
    jac_q_ = DaeMatrix(n_);

    // Every row gets a diagonal position (node-gmin augmentation + LU
    // robustness require it, even when no device writes it).
    for (int i = 0; i < n_; ++i) {
        jac_f_.declareEntry(i, i);
        jac_q_.declareEntry(i, i);
    }

    for (size_t i = 0; i < devices_.size(); ++i) {
        devices_[i]->registerJacobianPattern(rows[i], cols[i]);
        for (size_t k = 0; k < rows[i].size(); ++k) {
            if (rows[i][k] < 0 || cols[i][k] < 0) continue;  // ground/absent
            jac_f_.declareEntry(rows[i][k], cols[i][k]);
            jac_q_.declareEntry(rows[i][k], cols[i][k]);
        }
    }
    jac_f_.finalizePattern();
    jac_q_.finalizePattern();

    std::vector<JacobianEntry> valid_entries;
    std::vector<int> valid_positions;
    for (size_t i = 0; i < devices_.size(); ++i) {
        // Entries referencing a ground row/column (declared as -1) do not
        // correspond to a matrix position: they get slot -1 and are skipped
        // by emit(). The returned slots keep the ALIGNMENT of the device's
        // registration order.
        valid_entries.clear();
        valid_positions.clear();
        valid_entries.reserve(rows[i].size());
        for (size_t k = 0; k < rows[i].size(); ++k) {
            if (rows[i][k] >= 0 && cols[i][k] >= 0) {
                valid_entries.emplace_back(rows[i][k], cols[i][k]);
                valid_positions.push_back(static_cast<int>(k));
            }
        }
        std::vector<int> slots(rows[i].size(), -1);
        if (!valid_entries.empty()) {
            const std::vector<int> resolved = jac_f_.resolveSlots(valid_entries);
            for (size_t k = 0; k < resolved.size(); ++k) {
                slots[valid_positions[k]] = resolved[k];
            }
        }
        devices_[i]->setJacobianSlots(slots);
    }

    resolveDiagonalSlots();

    // ---- 6) work buffers ----
    F_.assign(n_, 0.0);
    Q_.assign(n_, 0.0);
    B_.assign(n_, 0.0);
    rhs_.assign(n_, 0.0);
    x_.assign(n_, 0.0);
    x_prev_.assign(n_, 0.0);
}

void DaeSystem::resolveDiagonalSlots() {
    std::vector<JacobianEntry> diag(n_);
    for (int i = 0; i < n_; ++i) diag[i] = JacobianEntry(i, i);
    diagonal_slots_ = jac_f_.resolveSlots(diag);
}

// ===============================================================
// ASSEMBLY
// ===============================================================

void DaeSystem::assembleVectors(double time, double q_coeff) {
    DaeSolution sol(x_, circuit_, node_map_);
    const size_t nd = devices_.size();

    std::fill(F_.begin(), F_.end(), 0.0);
    std::fill(Q_.begin(), Q_.end(), 0.0);
    std::fill(B_.begin(), B_.end(), 0.0);

    // Small circuits: the serial path avoids OpenMP overhead.
    if (nd <= 32) {
        for (Device* d : devices_) {
            d->loadF(sol, time, F_);
            if (q_coeff != 0.0) d->loadQ(sol, time, Q_);
            d->loadB(time, B_);
        }
        return;
    }

    // Large circuits: one scratch set per thread reduced once at the end.
    #pragma omp parallel
    {
        std::vector<double> f(n_, 0.0), q(n_, 0.0), b(n_, 0.0);
        #pragma omp for schedule(dynamic, 16) nowait
        for (int i = 0; i < static_cast<int>(nd); ++i) {
            devices_[i]->loadF(sol, time, f);
            if (q_coeff != 0.0) devices_[i]->loadQ(sol, time, q);
            devices_[i]->loadB(time, b);
        }
        #pragma omp critical
        {
            for (int i = 0; i < n_; ++i) {
                F_[i] += f[i];
                Q_[i] += q[i];
                B_[i] += b[i];
            }
        }
    }
}

void DaeSystem::assembleJacobians(double time, double q_coeff) {
    DaeSolution sol(x_, circuit_, node_map_);
    const size_t nd = devices_.size();

    jac_f_.zeroValues();
    jac_q_.zeroValues();

    if (nd <= 32) {
        for (Device* d : devices_) {
            d->loadJacobianF(jac_f_, sol, time);
            d->loadJacobianQ(jac_q_, sol, time);
        }
    } else {
        // Fully parallel: each device writes only its own NNZ slots.
        #pragma omp parallel for schedule(dynamic, 16)
        for (int i = 0; i < static_cast<int>(nd); ++i) {
            devices_[i]->loadJacobianF(jac_f_, sol, time);
            devices_[i]->loadJacobianQ(jac_q_, sol, time);
        }
    }

    // Effective Jacobian  J = dFdx + s·dQdx. Both matrices share the same
    // pattern (hence identical nnz ordering), so the fold is a single linear
    // pass over the value arrays — no maps, no reordering.
    if (q_coeff != 0.0) {
        std::vector<double>& jf = jac_f_.values();
        const std::vector<double>& jq = jac_q_.values();
        const int nnz = jac_f_.nonZeros();
        #pragma omp parallel for simd schedule(static)
        for (int k = 0; k < nnz; ++k) {
            jf[k] += q_coeff * jq[k];
        }
    }
}

void DaeSystem::applyGminResidual() {
    if (node_gmin_ == 0.0) return;
    // Homotopy: a gmin resistor from every node to ground adds gmin·v to the
    // KCL residual (must be applied BEFORE the RHS is assembled).
    for (int i = 0; i < num_nodes_; ++i) {
        F_[i] += node_gmin_ * x_[i];
    }
}

void DaeSystem::applyGminDiagonal() {
    if (node_gmin_ == 0.0) return;
    for (int i = 0; i < num_nodes_; ++i) {
        jac_f_.accumulateSlot(diagonal_slots_[i], node_gmin_);
    }
}

void DaeSystem::buildRHS(double q_coeff, bool use_history) {
    // RHS = −F − s·Q(x) + B, plus the charge history terms in transient.
    for (int i = 0; i < n_; ++i) {
        rhs_[i] = B_[i] - F_[i] - q_coeff * Q_[i];
    }
    if (use_history) {
        const double qdot_coeff = trapezoidal_ ? 1.0 : 0.0;
        for (const Device* d : devices_) {
            d->loadChargeHistory(store_, q_coeff, qdot_coeff, rhs_);
        }
    }
}

bool DaeSystem::checkConvergence(const std::vector<double>& x_new) const {
    // SPICE-style: |Δx| ≤ abstol + reltol·|x_new| for every unknown, and no
    // nonlinear device may reject the iteration (limiting active).
    for (const Device* d : devices_) {
        if (!d->acceptSolution()) return false;
    }
    for (int i = 0; i < n_; ++i) {
        const double delta = std::abs(x_new[i] - x_prev_[i]);
        if (delta > abstol_ + reltol_ * std::abs(x_new[i])) {
            return false;
        }
    }
    return true;
}

bool DaeSystem::solveLinearSystem(std::vector<double>& dx) {
    LinearSolver solver(solver_type_);
    solver.setMaxIterations(2000);
    solver.setTolerance(1e-12);
    const SolverResult r = solver.solve(jac_f_, rhs_, dx);
    if (r.converged) return true;

    // Robust fallback: the Krylov solver can break down on strongly
    // non-symmetric MNA Jacobians; a direct solve always succeeds/fails
    // cleanly. Used for small/medium systems (and to unstick the analysis).
    if (solver_type_ != SolverType::DIRECT_LU && n_ <= 20000) {
        LinearSolver direct(SolverType::DIRECT_LU);
        const SolverResult rd = direct.solve(jac_f_, rhs_, dx);
        return rd.converged;
    }
    return false;
}

// ===============================================================
// NEWTON LOOP
// ===============================================================

NLSolveResult DaeSystem::solveNonlinear(double time,
                                        double q_coeff,
                                        bool use_history,
                                        const std::vector<double>& x0) {
    NLSolveResult result;
    const auto t0 = std::chrono::high_resolution_clock::now();

    if (!x0.empty()) x_ = x0;
    std::vector<double> dx(n_, 0.0);

    // One full Newton attempt with a fixed node-gmin (homotopy level).
    auto attempt = [&](double gmin) -> bool {
        node_gmin_ = gmin;
        for (int iter = 1; iter <= max_iterations_; ++iter) {
            assembleVectors(time, q_coeff);
            applyGminResidual();          // F += gmin·x  (before the RHS)
            assembleJacobians(time, q_coeff);
            applyGminDiagonal();          // J += gmin     (after the stamp)
            buildRHS(q_coeff, use_history);

            if (!solveLinearSystem(dx)) {
                result.message = "linear solver failed";
                if (verbose_) std::cerr << "[DaeSystem] linear solve failed (gmin=" << gmin << ")\n";
                return false;
            }

            x_prev_ = x_;
            for (int i = 0; i < n_; ++i) x_[i] += dx[i];
            ++result.total_iterations;
            if (verbose_ && very_verbose_) {
                std::cerr << "[DaeSystem] gmin=" << gmin << " iter=" << iter
                          << " dx=[" << dx[0];
                for (int i = 1; i < n_ && i < 4; ++i) std::cerr << "," << dx[i];
                std::cerr << "]\n";
            }

            if (checkConvergence(x_)) {
                result.iterations = iter;
                return true;
            }
        }
        return false;
    };

    // 1) Direct Newton with the base node-gmin (Xyce runs its DCOP with a
    //    small gmin on every node by default).
    bool ok = attempt(gmin_final_);
    if (!ok) {
        // 2) GMIN homotopy: relax from a heavily conductive circuit down to
        //    the base level, warm-starting each solve from the previous one.
        std::vector<double> levels;
        for (double g = gmin_start_; g > gmin_final_; g *= 0.1) {
            levels.push_back(g);
        }
        levels.push_back(gmin_final_);

        bool final_ok = false;
        for (double g : levels) {
            if (!attempt(g)) break;   // stop at the first failing level
            ++result.gmin_steps;
            final_ok = true;
        }
        ok = final_ok;
        if (!ok) result.message = "Newton did not converge";
    }

    // Charge/state of the ACCEPTED point (post-update x).
    if (ok) {
        commitStates();
        for (Device* d : devices_) d->confirmAccept();
    }

    node_gmin_ = 0.0;
    result.converged = ok;    result.elapsed_seconds =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
    if (ok) result.message = "converged";
    return result;
}

void DaeSystem::commitStates() {
    DaeSolution sol(x_, circuit_, node_map_);
    for (Device* d : devices_) {
        d->storeState(sol, store_);
    }
}

void DaeSystem::acceptTimeStep(double q_coeff, bool trapezoidal) {
    store_.acceptStep(q_coeff, trapezoidal);
    for (Device* d : devices_) {
        d->commitStates(store_);
    }
}

void DaeSystem::extractSolution() {
    const auto& nodes = circuit_.getNodes();
    for (int i = 0; i < num_nodes_; ++i) {
        nodes[i]->setVoltage(x_[i]);
    }
}

void DaeSystem::saveAcceptedState() {
    x_saved_ = x_;
    store_.saveAccepted();
}

void DaeSystem::restoreAcceptedState() {
    x_ = x_saved_;
    store_.restoreAccepted();
    extractSolution();
}

void DaeSystem::beginSolve(bool initial_solve) {
    for (Device* d : devices_) {
        d->beginSolve(initial_solve);
    }
}

void DaeSystem::seedInitialGuess(std::vector<double>& x) const {
    for (const Device* d : devices_) {
        d->seedInitialGuess(x);
    }
}

void DaeSystem::setICOperatingPoint(bool on) {
    for (Device* d : devices_) {
        d->setICOperatingPoint(on);
    }
}

void DaeSystem::requestBreakpoints(std::vector<double>& times, double tstop) const {
    for (const Device* d : devices_) {
        d->requestBreakpoints(times, tstop);
    }
}

double DaeSystem::maxTimeStepHint(double t) const {
    double h = 0.0;   // 0 = unconstrained
    for (const Device* d : devices_) {
        const double hint = d->maxTimeStepHint(t);
        if (hint > 0.0 && (h == 0.0 || hint < h)) h = hint;
    }
    return h;
}

} // namespace e_XSpice

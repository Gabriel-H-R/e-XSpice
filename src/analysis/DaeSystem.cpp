/**
 * @author Gabriel Henrique Silva
 * @file DaeSystem.cpp
 * @brief Shared DAE/NL solver core for all analysis engines.
 */

#include "DaeSystem.h"
#include "../parallel/MPIManager.h"
#include "../parallel/DistributedKrylov.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <omp.h>
#include <iostream>
#include <numeric>

namespace e_XSpice {

DaeSystem::DaeSystem(Circuit& circuit)
    : circuit_(circuit) {}

// ===============================================================
// SETUP
// ===============================================================

void DaeSystem::setup() {
    // ---- 1) node map: ground is never materialized ----
    // The map is indexed by NODE ID (which SPICE netlists may space
    // sparsely: 1 2 3 11 12 ...), so it must be sized by the LARGEST id,
    // not by the number of nodes. Sizing it by count is a silent OOB write
    // on the first sparse id (heap corruption) — fixed here.
    const auto& nodes = circuit_.getNodes();
    num_nodes_ = static_cast<int>(nodes.size());
    int max_id = 0;
    for (const auto& n : nodes) {
        max_id = std::max(max_id, n->getId());
    }
    node_map_.assign(static_cast<size_t>(max_id) + 1, -1);
    for (int i = 0; i < num_nodes_; ++i) {
        node_map_[nodes[i]->getId()] = i;
    }

    // ---- 2) extra unknowns (branch currents) ----
    // Observer views (raw pointers): hot paths never copy an RCP.
    devices_.clear();
    devices_.reserve(circuit_.getDevices().size());
    for (const auto& d : circuit_.getDevices()) {
        devices_.push_back(d.get());
    }
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

    // ---- 7) distribution plan (identical on every rank, no comm) ----
    int rank = 0, size = 1;
    if (mpi_ && mpi_->isInitialized()) {
        rank = mpi_->getRank();
        size = mpi_->getSize();
    }
    plan_ = buildDistributionPlan(devices_, n_, rank, size,
                                  distMinDevices_, distMinRows_);
}

int DaeSystem::numLocalDevices() const {
    if (!plan_.assemblyDistributed) return static_cast<int>(devices_.size());
    int count = 0;
    for (size_t i = 0; i < devices_.size(); ++i) {
        if (ownsDevice(i)) ++count;
    }
    return count;
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

    // Evaluates the owned devices (all of them when serial). Each rank
    // accumulates only its own subset; syncAssemblyBuffers() fuses them.
    auto evalDevice = [&](size_t i) {
        Device* d = devices_[i];
        d->loadF(sol, time, F_);
        if (q_coeff != 0.0) d->loadQ(sol, time, Q_);
        d->loadB(time, B_);
    };

    // Small circuits: the serial path avoids OpenMP overhead.
    if (nd <= 32) {
        for (size_t i = 0; i < nd; ++i) {
            if (ownsDevice(i)) evalDevice(i);
        }
    } else {
        // Large circuits: one scratch set per thread reduced once at the end.
        #pragma omp parallel
        {
            std::vector<double> f(n_, 0.0), q(n_, 0.0), b(n_, 0.0);
            #pragma omp for schedule(dynamic, 16) nowait
            for (int i = 0; i < static_cast<int>(nd); ++i) {
                if (!ownsDevice(static_cast<size_t>(i))) continue;
                Device* d = devices_[i];
                d->loadF(sol, time, f);
                if (q_coeff != 0.0) d->loadQ(sol, time, q);
                d->loadB(time, b);
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

    syncVectorBuffers();
}

void DaeSystem::syncVectorBuffers() {
    // Fuse the per-rank partial F/Q/B into the identical global vectors.
    // Called EXACTLY ONCE per assembly (in assembleVectors).
    //
    // The three vectors travel in a SINGLE Allreduce (packed scratch buffer)
    // instead of three: collectives are latency-bound, so one round trip
    // instead of three is ~3x less synchronization overhead per Newton
    // iteration. The pack/unpack memcpys are bandwidth-cheap by comparison.
    if (!plan_.assemblyDistributed || !mpi_) return;
    if (packed_.size() != F_.size() + Q_.size() + B_.size()) {
        packed_.resize(F_.size() + Q_.size() + B_.size(), 0.0);
    }
    auto* dst = packed_.data();
    std::copy(F_.begin(), F_.end(), dst); dst += F_.size();
    std::copy(Q_.begin(), Q_.end(), dst); dst += Q_.size();
    std::copy(B_.begin(), B_.end(), dst);
    mpi_->allReduceSumInPlace(packed_);
    const auto* src = packed_.data();
    std::copy(src, src + F_.size(), F_.begin()); src += F_.size();
    std::copy(src, src + Q_.size(), Q_.begin()); src += Q_.size();
    std::copy(src, src + B_.size(), B_.begin());
}

void DaeSystem::syncJacobianBuffer() {
    // Fuse the per-rank partial Jacobian (already folded J = dFdx + s·dQdx).
    // Called EXACTLY ONCE per assembly (in assembleJacobians).
    if (!plan_.assemblyDistributed || !mpi_) return;
    mpi_->allReduceSumInPlace(jac_f_.values());
}

void DaeSystem::assembleJacobians(double time, double q_coeff) {
    DaeSolution sol(x_, circuit_, node_map_);
    const size_t nd = devices_.size();

    jac_f_.zeroValues();
    jac_q_.zeroValues();

    if (nd <= 32) {
        for (size_t i = 0; i < nd; ++i) {
            if (!ownsDevice(i)) continue;
            devices_[i]->loadJacobianF(jac_f_, sol, time);
            devices_[i]->loadJacobianQ(jac_q_, sol, time);
        }
    } else {
        // Fully parallel: each device writes only its own NNZ slots.
        #pragma omp parallel for schedule(dynamic, 16)
        for (int i = 0; i < static_cast<int>(nd); ++i) {
            if (!ownsDevice(static_cast<size_t>(i))) continue;
            devices_[i]->loadJacobianF(jac_f_, sol, time);
            devices_[i]->loadJacobianQ(jac_q_, sol, time);
        }
    }

    // Effective Jacobian  J = dFdx + s·dQdx. Both matrices share the same
    // pattern (hence identical nnz ordering), so the fold is a single linear
    // pass over the value arrays — no maps, no reordering. The fold happens
    // BEFORE the Allreduce so only one Jacobian crosses the network.
    if (q_coeff != 0.0) {
        std::vector<double>& jf = jac_f_.values();
        const std::vector<double>& jq = jac_q_.values();
        const int nnz = jac_f_.nonZeros();
        #pragma omp parallel for simd schedule(static)
        for (int k = 0; k < nnz; ++k) {
            jf[k] += q_coeff * jq[k];
        }
    }

    syncJacobianBuffer();
}

void DaeSystem::refreshLinearization(bool include_gmin) {
    assembleJacobians(currentTime_, 0.0);
    // The node gmin is the DC CONVERGENCE aid of the homotopy: it is what
    // makes a floating node solvable at the operating point. It is NOT part of
    // the device physics, so a consumer that wants the linearization itself
    // (the AC analysis) leaves it out — a gmin diagonal would make every node
    // diagonal "present" at 1e-12, which is a numerical trap for any
    // no-pivot ordering and changes the small-signal answer. A consumer that
    // wants to reproduce the Newton matrix can ask for it.
    if (include_gmin) applyGminDiagonal();
    // assembleJacobians only reduces dF/dx (the transient path consumes the
    // folded J = dF/dx + s·dQ/dx). The AC analysis consumes the pair
    // SEPARATELY as (G, C), so dQ/dx has to be reduced here too.
    if (plan_.assemblyDistributed && mpi_) {
        mpi_->allReduceSumInPlace(jac_q_.values());
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

void DaeSystem::applySolutionLimits() {
    // Only ever called on the replicated (serial) path — see solveNonlinear.
    for (Device* d : devices_) {
        d->applySolutionLimit(x_);
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

bool DaeSystem::checkDevicesAccept() const {
    // Every OWNED device accepts the current evaluation (no limiter fired
    // while assembling it). The single AND fuses the votes so every rank
    // decides identically.
    bool local = true;
    for (size_t i = 0; i < devices_.size(); ++i) {
        if (ownsDevice(i) && !devices_[i]->acceptSolution()) {
            local = false;
            break;
        }
    }
    if (plan_.assemblyDistributed && mpi_) {
        return mpi_->allReduceAnd(local);
    }
    return local;
}

bool DaeSystem::checkUpdateConverged(const std::vector<double>& x_new) const {
    // SPICE-style update test: |Δx| ≤ abstol + reltol·|x_new| for every
    // unknown, and no nonlinear device may reject the iteration (limiting
    // active). x itself is identical on every rank, so the delta test needs
    // no communication but is folded into the same vote.
    if (!checkDevicesAccept()) return false;
    bool local = true;
    for (int i = 0; i < n_; ++i) {
        const double delta = std::abs(x_new[i] - x_prev_[i]);
        if (delta > abstol_ + reltol_ * std::abs(x_new[i])) {
            local = false;
            break;
        }
    }
    if (plan_.assemblyDistributed && mpi_) {
        return mpi_->allReduceAnd(local);
    }
    return local;
}

bool DaeSystem::checkResidualConverged(double q_coeff) const {
    // SPICE-style residual test on the FRESHLY assembled rhs_ (B − F −
    // s·Q + history at the candidate point): max|rhs| must be small
    // relative to the equation scale max(|B| + |F| + s·|Q|). Catches the
    // false convergences a pure update test admits on (near-)singular
    // Jacobians — tiny dx coexisting with a big residual.
    double res = 0.0, scale = 0.0;
    for (int i = 0; i < n_; ++i) {
        res = std::max(res, std::abs(rhs_[i]));
        scale = std::max(scale,
                         std::abs(B_[i]) + std::abs(F_[i]) +
                         std::abs(q_coeff) * std::abs(Q_[i]));
    }
    bool local = (res <= abstol_ + reltol_ * scale);
    if (plan_.assemblyDistributed && mpi_) {
        return mpi_->allReduceAnd(local);
    }
    return local;
}

void DaeSystem::dampNewtonStep(std::vector<double>& dx, double trust) const {
    // Adaptive trust-region safeguard, node voltages only (branch currents
    // are Lagrange multipliers enforcing stiff constraints — capping them
    // breaks enforcement; the distributed path skips damping entirely
    // since dx is row-partitioned there).
    if (plan_.solveDistributed) return;
    double peak = 0.0;
    for (int i = 0; i < num_nodes_; ++i) {
        peak = std::max(peak, std::abs(dx[i]));
    }
    if (peak > trust) {
        const double s = trust / peak;
        for (int i = 0; i < num_nodes_; ++i) dx[i] *= s;
    }
}

bool DaeSystem::solveLinearSystem(std::vector<double>& dx) {
    ++solveStats_.linearSolves;
    if (plan_.solveDistributed && mpi_) {
        ++solveStats_.distributedSolves;
        SolverResult r = distributedBiCGSTAB(jac_f_, rhs_, dx,
                                             plan_.rowBegin, plan_.rowEnd,
                                             *mpi_, 2000, 1e-12);
        solveStats_.krylovIterations += r.iterations;
        if (r.converged) return true;
        // Fall through to the replicated direct fallback (the fused Jacobian
        // is identical on every rank, so no communication is needed).
    } else if (n_ <= DirectSolveMaxRows_ && solver_type_ != SolverType::DIRECT_LU) {
        // Small/medium systems: go straight to the direct solve. A Krylov
        // solver on a nearly singular Newton Jacobian (floating node held
        // only by gmin, zero-capacitance companion, ...) can "converge" to
        // a point with arbitrary nullspace components — e.g. a 3.9 V first
        // step on a row whose RHS is ~0 — while LU returns the exact
        // minimum-norm-free deterministic answer (dx ≈ 0 there). Direct is
        // also cheaper here than up-to-2000 Krylov iterations per Newton
        // step. Krylov stays the workhorse above this size (and when the
        // user explicitly selects DIRECT_LU, that is honored as-is).
        LinearSolver direct(SolverType::DIRECT_LU);
        const SolverResult rd = direct.solve(jac_f_, rhs_, dx);
        if (rd.converged) ++solveStats_.directFallbacks;
        return rd.converged;
    } else {
        LinearSolver solver(solver_type_);
        solver.setMaxIterations(2000);
        solver.setTolerance(1e-12);
        const SolverResult r = solver.solve(jac_f_, rhs_, dx);
        solveStats_.krylovIterations += r.iterations;
        if (r.converged) return true;
    }

    // Robust fallback: the Krylov solver can break down on strongly
    // non-symmetric MNA Jacobians; a direct solve always succeeds/fails
    // cleanly. Used for small/medium systems (and to unstick the analysis).
    if (solver_type_ != SolverType::DIRECT_LU && n_ <= 20000) {
        LinearSolver direct(SolverType::DIRECT_LU);
        const SolverResult rd = direct.solve(jac_f_, rhs_, dx);
        if (rd.converged) ++solveStats_.directFallbacks;
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
    currentTime_ = time;

    if (!x0.empty()) x_ = x0;
    std::vector<double> dx(n_, 0.0);
    // Sane entry point (user guess, previous accepted step, or zeros):
    // the homotopy below restarts from here if direct Newton diverges.
    const std::vector<double> x_sane_entry_ = x_;

    // One full Newton attempt with a fixed node-gmin (homotopy level).
    auto attempt = [&](double gmin) -> bool {
        node_gmin_ = gmin;
        // Adaptive trust region (per attempt): start at the ceiling; grow
        // ×2 on residual progress, shrink ×0.5 when the residual grows.
        // Lets high-voltage circuits take big early steps while herding
        // singular-direction excursions (floating nodes, exponential
        // overshoots) down to steps the linearization can honor.
        double trust = MaxNewtonStep_;
        double prevRes = std::numeric_limits<double>::infinity();
        for (int iter = 1; iter <= max_iterations_; ++iter) {
            assembleVectors(time, q_coeff);
            applyGminResidual();          // F += gmin·x  (before the RHS)
            assembleJacobians(time, q_coeff);
            applyGminDiagonal();          // J += gmin     (after the stamp)
            buildRHS(q_coeff, use_history);

            double resNow = 0.0;
            for (int i = 0; i < n_; ++i) {
                resNow = std::max(resNow, std::abs(rhs_[i]));
            }
            if (resNow < prevRes) {
                trust = std::min(MaxNewtonStep_, trust * 2.0);
            } else {
                trust = std::max(MinNewtonStep_, trust * 0.5);
            }
            prevRes = resNow;

            if (!solveLinearSystem(dx)) {
                result.message = "linear solver failed";
                if (verbose_) std::cerr << "[DaeSystem] linear solve failed (gmin=" << gmin << ")\n";
                return false;
            }

            x_prev_ = x_;
            dampNewtonStep(dx, trust);
            for (int i = 0; i < n_; ++i) x_[i] += dx[i];
            ++result.total_iterations;
            if (std::getenv("EXSPICE_TRACE_NEWTON")) {
                std::fprintf(stderr, "[NWT] iter=%d x=[", iter);
                for (int i = 0; i < n_; ++i) {
                    std::fprintf(stderr, "%s%.9g", i ? "," : "", x_[i]);
                }
                std::fprintf(stderr, "] dx=[");
                for (int i = 0; i < n_; ++i) {
                    std::fprintf(stderr, "%s%.9g", i ? "," : "", dx[i]);
                }
                double fn = 0.0;
                for (int i = 0; i < n_; ++i) fn += rhs_[i] * rhs_[i];
                std::fprintf(stderr, "] |rhs|=%.9g\n", std::sqrt(fn));
                if (iter <= 2) {
                    std::fprintf(stderr, "[NWT] rhs=[");
                    for (int i = 0; i < n_; ++i) {
                        std::fprintf(stderr, "%s%.9g", i ? "," : "", rhs_[i]);
                    }
                    std::fprintf(stderr, "]\n[NWT] F=[");
                    for (int i = 0; i < n_; ++i) {
                        std::fprintf(stderr, "%s%.9g", i ? "," : "", F_[i]);
                    }
                    std::fprintf(stderr, "]\n");
                }
            }
            if (verbose_ && very_verbose_) {
                std::cerr << "[DaeSystem] gmin=" << gmin << " iter=" << iter
                          << " dx=[" << dx[0];
                for (int i = 1; i < n_ && i < 4; ++i) std::cerr << "," << dx[i];
                std::cerr << "]\n";
            }

            if (checkUpdateConverged(x_)) {
                // Update test passed — verify it is a TRUE solution: with a
                // (near-)singular Jacobian a tiny dx can coexist with a big
                // residual (observed: "converged" to a KCL-violating point
                // on a floating inverter output). Reassemble the residual at
                // the new point and require it small relative to the
                // equation scale, SPICE-style. Costs one extra assembly per
                // solve (only on the accepting iteration, plus rare false
                // alarms) — negligible next to the Newton loop.
                assembleVectors(time, q_coeff);
                applyGminResidual();
                buildRHS(q_coeff, use_history);
                if (checkResidualConverged(q_coeff) && checkDevicesAccept()) {
                    result.iterations = iter;
                    return true;
                }
                // Not a true solution yet: Jacobians are stale now (they
                // belong to the pre-step point), so reassemble everything
                // next pass through the loop head. The vectors/RHS above
                // are already fresh — skip straight to the Jacobian stage.
                assembleJacobians(time, q_coeff);
                applyGminDiagonal();
                buildRHS(q_coeff, use_history);
                if (!solveLinearSystem(dx)) {
                    result.message = "linear solver failed";
                    return false;
                }
                x_prev_ = x_;
                dampNewtonStep(dx, trust);
                for (int i = 0; i < n_; ++i) x_[i] += dx[i];
                ++result.total_iterations;
                if (!plan_.assemblyDistributed) applySolutionLimits();
                continue;
            }

            // SPICE/Xyce limiting feedback: devices that rejected the
            // iteration write their limited unknowns back into x, so the
            // next iteration restarts AT the limited point instead of
            // walking the exponential up by vt per iteration. Serial path
            // only: the distributed path must keep x identical on every
            // rank (per-device limiting flags are rank-local by design).
            if (!plan_.assemblyDistributed) applySolutionLimits();
        }
        return false;
    };

    // 1) Direct Newton with the base node-gmin (Xyce runs its DCOP with a
    //    small gmin on every node by default).
    bool ok = attempt(gmin_final_);
    if (!ok) {
        // 2) GMIN homotopy: relax from a heavily conductive circuit down to
        //    the base level, warm-starting each solve from the previous one.
        //    CRITICAL: restart from the PRE-ATTEMPT point, not from the
        //    diverged garbage the failed direct Newton left behind — warming
        //    homotopy from garbage (e.g. an output node at 104 V) tracks the
        //    wrong solution branch, while the sane entry point tracks the
        //    physical one. Limiter histories are re-anchored to match.
        x_ = x_sane_entry_;
        {
            const DaeSolution sol(x_, circuit_, node_map_);
            for (Device* d : devices_) {
                d->resetIterationState(sol);
            }
        }
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
        for (size_t i = 0; i < devices_.size(); ++i) {
            if (ownsDevice(i)) devices_[i]->confirmAccept();
        }
    }

    node_gmin_ = 0.0;
    result.converged = ok;    result.elapsed_seconds =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
    if (ok) result.message = "converged";
    return result;
}

void DaeSystem::commitStates() {
    DaeSolution sol(x_, circuit_, node_map_);
    // CRITICAL: clear the whole candidate history first. Each rank writes
    // only its OWNED slots below; without this, unowned slots keep the
    // previous step's fused values and the SUM fusion below would add stale
    // history on top of fresh history (double counting that grows every
    // step and destabilizes the trapezoidal qdot feedback).
    std::fill(store_.nextStateData().begin(), store_.nextStateData().end(),
              0.0);
    for (size_t i = 0; i < devices_.size(); ++i) {
        if (ownsDevice(i)) devices_[i]->storeState(sol, store_);
    }
    syncStateHistory();
}

void DaeSystem::syncStateHistory() {
    // Each rank wrote Q(x) only for its owned devices (other slots are zero
    // locally); one SUM fuses the identical candidate history everywhere.
    if (!plan_.assemblyDistributed || !mpi_) return;
    mpi_->allReduceSumInPlace(store_.nextStateData());
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
    // Refresh every device's cached outputs from the identical global
    // solution (unowned devices would otherwise report stale values).
    DaeSolution sol(x_, circuit_, node_map_);
    for (Device* d : devices_) {
        d->refreshOutputs(sol, currentTime_);
    }
}

void DaeSystem::saveAcceptedState() {
    x_saved_ = x_;
    store_.saveAccepted();
}

void DaeSystem::restoreAcceptedState() {
    x_ = x_saved_;
    store_.restoreAccepted();
    const DaeSolution sol(x_, circuit_, node_map_);
    for (Device* d : devices_) {
        d->resetIterationState(sol);
    }
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

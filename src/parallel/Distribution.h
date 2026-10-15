/**
 * @author Gabriel Henrique Silva
 * @file Distribution.h
 * @brief Work distribution plan for SPMD parallel simulation (MPI).
 *
 * Design (modular, no MPI dependency in this module):
 * ---------------------------------------------------
 * The simulator runs in SPMD mode: every rank builds the SAME circuit,
 * the SAME fixed sparsity pattern and holds row-REPLICATED assembly
 * buffers (F, Q, B, J values, state history). Parallelism comes from
 * partitioning the WORK, not the data:
 *
 *  1. DEVICE-PARTITIONED ASSEMBLY: each rank evaluates a disjoint subset
 *     of the devices (loadF/loadQ/loadB/loadJacobians/storeState) into its
 *     local buffers, then a single SUM Allreduce per buffer fuses the
 *     identical global system on every rank. Because the CSR pattern is
 *     deterministic (same declaration order everywhere), the value arrays
 *     are element-wise combinable with zero reordering.
 *
 *  2. ROW-DISTRIBUTED KRYLOV SOLVE: each rank owns a contiguous block of
 *     matrix rows. The solution vector stays replicated (Newton updates are
 *     identical everywhere), so a matvec needs only an Allgatherv of the
 *     locally-computed block plus Allreduced dots — no ghost exchange,
 *     no distributed matrix storage.
 *
 *  3. IDENTICAL CONTROL FLOW: Newton/LTE/dt decisions derive only from
 *     globally-synchronized quantities, so every rank takes the same
 *     branches with no extra coordination (limiting votes use one AND).
 *
 * Load model (per device evaluation cost, normalized to a resistor):
 *   R/I/V-source = 1, C/L = 2 (state + history), diode = 12 (exp + pnjlim).
 * Assignment is greedy LPT (longest-processing-time first): devices sorted
 * by cost, each placed on the currently least-loaded rank. Deterministic,
 * O(d log d + d·p), runs identically on every rank.
 *
 * Thresholds: distribution only pays off above a minimum problem size
 * (collectives have latency). Below the thresholds every rank simply runs
 * the serial path redundantly — always correct, zero communication.
 */

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>

namespace e_XSpice {

class Device;

/**
 * @brief Distribution plan: identical on every rank, no communication
 *        needed to build it.
 */
struct DistributionPlan {
    bool assemblyDistributed = false;  // partition devices + Allreduce buffers
    bool solveDistributed = false;     // row-block distributed Krylov solve

    std::vector<int> deviceOwner;      // per device index -> owning rank
    int rowBegin = 0;                  // this rank's Krylov row block [begin,end)
    int rowEnd = 0;

    int rank = 0;
    int size = 1;

    bool ownsDevice(int deviceIndex) const {
        return !assemblyDistributed ||
               (deviceIndex >= 0 &&
                deviceIndex < static_cast<int>(deviceOwner.size()) &&
                deviceOwner[deviceIndex] == rank);
    }
};

/**
 * @brief Builds the distribution plan (pure function — no MPI calls).
 *
 * @param devices     all devices in deterministic order (same on every rank)
 * @param numRows     system dimension n
 * @param rank        this rank (0 when serial)
 * @param size        number of ranks (1 when serial)
 * @param minDevices  minimum device count to enable distributed assembly
 * @param minRows     minimum system dimension to enable distributed solve
 */
DistributionPlan buildDistributionPlan(const std::vector<Device*>& devices,
                                       int numRows,
                                       int rank, int size,
                                       int minDevices, int minRows);

/**
 * @brief Reads the EXSPICE_DISTRIBUTE environment override:
 *        "1" forces distribution on, "0" forces serial, unset/else = auto.
 * @return 1 (force on), 0 (force off), -1 (auto)
 */
int distributionEnvOverride();

} // namespace e_XSpice

#endif // DISTRIBUTION_H

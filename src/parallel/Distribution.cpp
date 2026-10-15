/**
 * @author Gabriel Henrique Silva
 * @file Distribution.cpp
 * @brief Work distribution plan for SPMD parallel simulation (MPI).
 */

#include "Distribution.h"
#include "../core/Device.h"

#include <algorithm>
#include <cstdlib>
#include <numeric>

namespace e_XSpice {

namespace {

// Estimated per-device evaluation cost (resistor = 1 unit).
int deviceCost(const Device* d) {
    switch (d->getType()) {
        case DeviceType::DIODE:         return 12;  // exp + pnjlim + charge
        case DeviceType::CAPACITOR:     return 2;
        case DeviceType::INDUCTOR:       return 2;
        case DeviceType::MOSFET:         return 20;  // future-proof
        case DeviceType::BJT:            return 16;  // future-proof
        case DeviceType::RESISTOR:
        case DeviceType::VOLTAGE_SOURCE:
        case DeviceType::CURRENT_SOURCE:
        default:                         return 1;
    }
}

// Contiguous block partition (same layout as MPIManager::getLocalRange,
// duplicated here to keep this module MPI-free).
void blockRange(int total, int size, int rank, int& start, int& end) {
    if (total <= 0 || size <= 0) { start = end = 0; return; }
    const int base = total / size;
    const int rem  = total % size;
    if (rank < rem) {
        start = rank * (base + 1);
        end   = start + base + 1;
    } else {
        start = rank * base + rem;
        end   = start + base;
    }
}

} // namespace

int distributionEnvOverride() {
    const char* env = std::getenv("EXSPICE_DISTRIBUTE");
    if (!env) return -1;
    if (env[0] == '1') return 1;
    if (env[0] == '0') return 0;
    return -1;
}

DistributionPlan buildDistributionPlan(const std::vector<Device*>& devices,
                                       int numRows,
                                       int rank, int size,
                                       int minDevices, int minRows) {
    DistributionPlan plan;
    plan.rank = rank;
    plan.size = size;

    const int nd = static_cast<int>(devices.size());
    plan.deviceOwner.assign(nd, 0);

    if (size <= 1 || nd == 0) {
        plan.rowBegin = 0;
        plan.rowEnd = numRows;
        return plan;   // serial: ownsDevice() is true for everything
    }

    const int env = distributionEnvOverride();
    const bool forceOn  = (env == 1);
    const bool forceOff = (env == 0);

    plan.assemblyDistributed =
        forceOn || (!forceOff && nd >= minDevices && size > 1);
    plan.solveDistributed =
        forceOn || (!forceOff && numRows >= minRows && size > 1);

    if (!plan.assemblyDistributed && !plan.solveDistributed) {
        plan.rowBegin = 0;
        plan.rowEnd = numRows;
        return plan;
    }

    // ---- greedy LPT device assignment (deterministic on every rank) ----
    if (plan.assemblyDistributed) {
        std::vector<int> order(nd);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            const int ca = deviceCost(devices[a]);
            const int cb = deviceCost(devices[b]);
            return (ca != cb) ? (ca > cb) : (a < b);
        });

        std::vector<long long> load(size, 0);
        for (int idx : order) {
            int best = 0;
            for (int r = 1; r < size; ++r) {
                if (load[r] < load[best]) best = r;
            }
            plan.deviceOwner[idx] = best;
            load[best] += deviceCost(devices[idx]);
        }
    }
    // else: deviceOwner stays 0 (unused — ownsDevice() short-circuits).

    // ---- contiguous Krylov row block for this rank ----
    blockRange(numRows, size, rank, plan.rowBegin, plan.rowEnd);

    return plan;
}

} // namespace e_XSpice

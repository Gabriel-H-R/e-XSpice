/**
 * @author Gabriel Henrique Silva
 * @file Inductor.cpp
 * @brief Inductor model — flux-based DAE formulation with branch current.
 */

#include "Inductor.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"
#include <stdexcept>

namespace e_XSpice {

Inductor::Inductor(const std::string& name,
                   Node* node_plus,
                   Node* node_minus,
                   double inductance,
                   double initial_current,
                   bool ic_given)
    : Device(name, DeviceType::INDUCTOR)
    , inductance_(inductance)
    , initial_current_(initial_current)
    , ic_given_(ic_given)
{
    if (inductance <= 0.0) {
        throw std::invalid_argument("Inductor inductance must be positive");
    }

    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    if (node_plus)  node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

// ---------------------------------------------------------------
// Setup
// ---------------------------------------------------------------

int Inductor::registerStates(StateStore& store) {
    setStateIndex(store.registerStates(1));
    return 1;
}

void Inductor::registerJacobianPattern(std::vector<int>& rows,
                                       std::vector<int>& cols) const {
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const int ib = branch_index_;

    // Fixed 3-entity stamp (Xyce jacStamp_BASE):
    //   row pos: {bra}          row neg: {bra}          row bra: {pos, neg, bra}
    rows.push_back(ip); cols.push_back(ib);  // 0: [p][bra]
    rows.push_back(in); cols.push_back(ib);  // 1: [n][bra]
    rows.push_back(ib); cols.push_back(ip);  // 2: [bra][p]
    rows.push_back(ib); cols.push_back(in);  // 3: [bra][n]
    rows.push_back(ib); cols.push_back(ib);  // 4: [bra][bra]
}

void Inductor::seedInitialGuess(std::vector<double>& x) const {
    if (branch_index_ >= 0 && branch_index_ < static_cast<int>(x.size())) {
        // IC-OP: pin the inductor current at IC (Xyce seeds the same way).
        x[branch_index_] = isICOperatingPoint() ? initial_current_ : 0.0;
    }
}

// ---------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------

void Inductor::loadF(const DaeSolution& x, double time,
                     std::vector<double>& F) const {
    (void)time;
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const int ib = branch_index_;

    if (isICOperatingPoint() && ic_given_) {
        // Current source of value IC; branch row pinned to itself.
        if (ip >= 0) F[ip] += initial_current_;
        if (in >= 0) F[in] -= initial_current_;
        F[ib] += 0.0;
        last_current_ = initial_current_;
        return;
    }

    // KCL rows carry the branch current; KVL row carries the voltage drop.
    const double iL = x.data()[ib];
    last_current_ = iL;

    vcap_cache_ = (ip >= 0 ? x.data()[ip] : 0.0)
                - (in >= 0 ? x.data()[in] : 0.0);

    if (ip >= 0) F[ip] += iL;
    if (in >= 0) F[in] -= iL;
    F[ib] += -vcap_cache_;
}

void Inductor::loadQ(const DaeSolution& x, double time,
                     std::vector<double>& Q) const {
    (void)time;
    if (isICOperatingPoint()) return;

    const double iL = x.data()[branch_index_];
    last_current_ = iL;
    flux_cache_ = inductance_ * iL;
    Q[branch_index_] += flux_cache_;
}

void Inductor::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                             double time) const {
    (void)x; (void)time;
    if (isICOperatingPoint() && ic_given_) {
        // Pin i_L = IC at the operating point.
        emit(J, 4, 1.0);
        return;
    }

    emit(J, 0,  1.0);   // [p][bra]
    emit(J, 1, -1.0);   // [n][bra]
    emit(J, 2, -1.0);   // [bra][p]
    emit(J, 3,  1.0);   // [bra][n]
}

void Inductor::loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                             double time) const {
    (void)x; (void)time;
    if (isICOperatingPoint()) return;

    emit(J, 4, inductance_);   // [bra][bra] += L
}

void Inductor::storeState(const DaeSolution& x, StateStore& store) const {
    const double iL  = x.data()[branch_index_];
    const double f0  = inductance_ * iL;
    flux_cache_ = f0;
    store.setDeviceQ(state_index_, &f0, 1);
}

void Inductor::loadChargeHistory(const StateStore& store, double q_coeff,
                                 double qdot_coeff,
                                 std::vector<double>& rhs) const {
    if (q_coeff == 0.0 || branch_index_ < 0) return;

    const double f0n = store.deviceQ(state_index_, 0);        // L·i(t_n)
    const double idn = qdot_coeff != 0.0
                       ? store.deviceQDot(state_index_, 0)   // dΦ/dt(t_n)
                       : 0.0;
    rhs[branch_index_] += q_coeff * f0n + idn;
}

void Inductor::commitStates(const StateStore& store) {
    // The inductor current is the branch unknown itself; the state slot
    // (flux L·i) is only used by the integrator. Nothing to cache here.
    (void)store;
}

double Inductor::getPower() const {
    return vcap_cache_ * last_current_;
}

} // namespace e_XSpice

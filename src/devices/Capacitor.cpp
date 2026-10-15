/**
 * @author Gabriel Henrique Silva
 * @file Capacitor.cpp
 * @brief Capacitor model — charge-based DAE formulation.
 */

#include "Capacitor.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"
#include <stdexcept>

namespace e_XSpice {

Capacitor::Capacitor(const std::string& name,
                     Node* node_plus,
                     Node* node_minus,
                     double capacitance,
                     double initial_voltage,
                     bool ic_given)
    : Device(name, DeviceType::CAPACITOR)
    , capacitance_(capacitance)
    , initial_voltage_(initial_voltage)
    , ic_given_(ic_given)
{
    if (capacitance <= 0.0) {
        throw std::invalid_argument("Capacitor capacitance must be positive");
    }

    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    if (node_plus)  node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

// ---------------------------------------------------------------
// Setup
// ---------------------------------------------------------------

int Capacitor::registerStates(StateStore& store) {
    setStateIndex(store.registerStates(1));
    return 1;
}

void Capacitor::registerJacobianPattern(std::vector<int>& rows,
                                        std::vector<int>& cols) const {
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const int ib = branch_index_;

    // Fixed 3x3 stamp (same shape as Xyce jacStamp_IC): the union of the
    // dQdx 2x2 block and the IC-OP branch equation.
    //   0: [p][p]  1: [p][n]  2: [p][bra]
    //   3: [n][p]  4: [n][n]  5: [n][bra]
    //   6: [bra][p] 7: [bra][n] 8: [bra][bra]
    rows.push_back(ip); cols.push_back(ip);  // 0
    rows.push_back(ip); cols.push_back(in);  // 1
    rows.push_back(ip); cols.push_back(ib);  // 2
    rows.push_back(in); cols.push_back(ip);  // 3
    rows.push_back(in); cols.push_back(in);  // 4
    rows.push_back(in); cols.push_back(ib);  // 5
    rows.push_back(ib); cols.push_back(ip);  // 6
    rows.push_back(ib); cols.push_back(in);  // 7
    rows.push_back(ib); cols.push_back(ib);  // 8
}

// ---------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------

void Capacitor::refreshVcap(const DaeSolution& x) const {
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    vcap_cache_ = (ip >= 0 ? x.data()[ip] : 0.0)
                - (in >= 0 ? x.data()[in] : 0.0);
}

void Capacitor::loadF(const DaeSolution& x, double time,
                       std::vector<double>& F) const {
    (void)time;
    if (!isICOperatingPoint()) return;   // pure charge device otherwise

    refreshVcap(x);
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const int ib = branch_index_;

    // IC-OP: capacitor behaves as a voltage source of value IC
    if (ip >= 0) F[ip] += x.data()[ib];
    if (in >= 0) F[in] -= x.data()[ib];
    F[ib] += vcap_cache_ - initial_voltage_;
}

void Capacitor::loadQ(const DaeSolution& x, double time,
                      std::vector<double>& Q) const {
    (void)time;
    if (inICMode()) return;   // at the IC operating point: no charge load

    refreshVcap(x);
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const double q0 = capacitance_ * vcap_cache_;

    if (ip >= 0) Q[ip] += q0;
    if (in >= 0) Q[in] -= q0;
}

void Capacitor::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                              double time) const {
    (void)x; (void)time;
    if (inICMode()) {
        // IC-OP: 0-V source (branch equation v+ − v− = IC lives in F)
        emit(J, 2,  1.0);   // [p][bra]
        emit(J, 5, -1.0);   // [n][bra]
        emit(J, 6,  1.0);   // [bra][p]
        emit(J, 7, -1.0);   // [bra][n]
    } else if (ic_given_) {
        // Transient with IC: the branch unknown is deactivated; a 1 on its
        // diagonal keeps the (otherwise empty) row non-singular. (Xyce:
        // dFdx[li_Bra][ABraEquBraNodeOffset] += 1.0 in this situation.)
        emit(J, 8, 1.0);
    }
    // Normal operation: no static contribution at all.
}

void Capacitor::loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                              double time) const {
    (void)x; (void)time;
    if (inICMode()) return;   // no charge term at the IC-OP

    emit(J, 0,  capacitance_);   // [p][p]
    emit(J, 1, -capacitance_);   // [p][n]
    emit(J, 3, -capacitance_);   // [n][p]
    emit(J, 4,  capacitance_);   // [n][n]
}

void Capacitor::storeState(const DaeSolution& x, StateStore& store) const {
    refreshVcap(x);
    const double q0 = capacitance_ * vcap_cache_;
    store.setDeviceQ(state_index_, &q0, 1);
}

void Capacitor::loadChargeHistory(const StateStore& store, double q_coeff,
                                  double qdot_coeff,
                                  std::vector<double>& rhs) const {
    if (q_coeff == 0.0 || !state_index_.isValid()) return;

    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const double q0n  = store.deviceQ(state_index_, 0);      // Q(t_n)
    const double iqn  = qdot_coeff != 0.0
                        ? store.deviceQDot(state_index_, 0)  // dQ/dt(t_n)
                        : 0.0;
    const double qterm = q_coeff * q0n + iqn;

    if (ip >= 0) rhs[ip] += qterm;
    if (in >= 0) rhs[in] -= qterm;
}

void Capacitor::commitStates(const StateStore& store) {
    last_current_ = store.deviceQDot(state_index_, 0);
}

void Capacitor::seedInitialGuess(std::vector<double>& x) const {
    // IC-OP: the branch current is unknown but starts from zero current.
    if (isICOperatingPoint() && branch_index_ >= 0 &&
        branch_index_ < static_cast<int>(x.size())) {
        x[branch_index_] = 0.0;
    }
}

double Capacitor::getPower() const {
    return getVoltage() * last_current_;
}

} // namespace e_XSpice

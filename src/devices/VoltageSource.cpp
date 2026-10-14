/**
 * @author Gabriel Henrique Silva
 * @file VoltageSource.cpp
 * @brief Independent voltage source — MNA with branch-current unknown.
 */

#include "VoltageSource.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"

namespace e_XSpice {

VoltageSource::VoltageSource(const std::string& name,
                             Node* node_plus,
                             Node* node_minus,
                             SourceWaveform waveform)
    : Device(name, DeviceType::VOLTAGE_SOURCE)
    , waveform_(std::move(waveform))
{
    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    if (node_plus)  node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

// ---------------------------------------------------------------
// Setup
// ---------------------------------------------------------------

void VoltageSource::registerJacobianPattern(std::vector<int>& rows,
                                            std::vector<int>& cols) const {
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const int ib = branch_index_;

    rows.push_back(ip); cols.push_back(ib);  // 0: [p][bra]
    rows.push_back(in); cols.push_back(ib);  // 1: [n][bra]
    rows.push_back(ib); cols.push_back(ip);  // 2: [bra][p]
    rows.push_back(ib); cols.push_back(in);  // 3: [bra][n]
}

void VoltageSource::seedInitialGuess(std::vector<double>& x) const {
    if (branch_index_ >= 0 && branch_index_ < static_cast<int>(x.size())) {
        x[branch_index_] = 0.0;
    }
}

// ---------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------

void VoltageSource::loadF(const DaeSolution& x, double time,
                          std::vector<double>& F) const {
    (void)time;
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const int ib = branch_index_;

    const double i = x.data()[ib];
    last_current_ = i;
    last_vcap_ = (ip >= 0 ? x.data()[ip] : 0.0)
               - (in >= 0 ? x.data()[in] : 0.0);
    last_time_ = time;

    if (ip >= 0) F[ip] += i;
    if (in >= 0) F[in] -= i;
    F[ib] += last_vcap_;
}

void VoltageSource::loadB(double time, std::vector<double>& B) const {
    B[branch_index_] += waveform_.value(time);
}

void VoltageSource::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                                  double time) const {
    (void)x; (void)time;
    emit(J, 0,  1.0);   // [p][bra]
    emit(J, 1, -1.0);   // [n][bra]
    emit(J, 2,  1.0);   // [bra][p]
    emit(J, 3, -1.0);   // [bra][n]
}

void VoltageSource::requestBreakpoints(std::vector<double>& times,
                                       double tstop) const {
    waveform_.breakpoints(times, tstop);
}

double VoltageSource::maxTimeStepHint(double t) const {
    return waveform_.maxTimeStep(t);
}

double VoltageSource::getPower() const {
    // SPICE convention: the branch unknown is the current flowing from the
    // + node to the − node INSIDE the source ("into" the + terminal).
    // Power DELIVERED to the circuit is therefore  −V·i.
    return -getVoltage() * last_current_;
}

} // namespace e_XSpice

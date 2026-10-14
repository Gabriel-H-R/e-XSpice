/**
 * @author Gabriel Henique Silva
 * @file CurrentSource.cpp
 * @brief Independent current source — pure RHS injection (B vector).
 */

#include "CurrentSource.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"

namespace e_XSpice {

CurrentSource::CurrentSource(const std::string& name,
                             Node* node_plus,
                             Node* node_minus,
                             SourceWaveform waveform)
    : Device(name, DeviceType::CURRENT_SOURCE)
    , waveform_(std::move(waveform))
{
    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    if (node_plus)  node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

void CurrentSource::registerJacobianPattern(std::vector<int>& rows,
                                            std::vector<int>& cols) const {
    (void)rows; (void)cols;   // no matrix entries at all
}

void CurrentSource::loadB(double time, std::vector<double>& B) const {
    last_time_ = time;
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const double i = waveform_.value(time);

    if (ip >= 0) B[ip] -= i;
    if (in >= 0) B[in] += i;
}

void CurrentSource::requestBreakpoints(std::vector<double>& times,
                                       double tstop) const {
    waveform_.breakpoints(times, tstop);
}

double CurrentSource::maxTimeStepHint(double t) const {
    return waveform_.maxTimeStep(t);
}

double CurrentSource::getPower() const {
    const double vp = nodes_[0] ? nodes_[0]->getVoltage() : 0.0;
    const double vn = nodes_[1] ? nodes_[1]->getVoltage() : 0.0;
    // The current flows from node_plus to node_minus INSIDE the source, so
    // the power delivered to the circuit is −I·(v+ − v−).
    return -waveform_.value(last_time_) * (vp - vn);
}

} // namespace e_XSpice

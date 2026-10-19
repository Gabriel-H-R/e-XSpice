/**
 * @author Gabriel Henrique Silva
 * @file Resistor.cpp
 * @brief Linear resistor model (DAE formulation).
 */

#include "Resistor.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"
#include <stdexcept>

namespace e_XSpice {

Resistor::Resistor(const std::string& name,
                   Node* node_plus,
                   Node* node_minus,
                   double resistance,
                   RCP<const ResistorModel> model)
    : Device(name, DeviceType::RESISTOR)
    , resistance_(resistance)
    , conductance_(0.0)
    , model_(std::move(model))
{
    if (resistance <= 0.0) {
        throw std::invalid_argument("Resistor resistance must be positive");
    }
    // Effective value folds the shared card in ONCE (setup-time). The
    // conductance below is what the DAE assembly reads — hot paths are
    // untouched by the model layer.
    double effective = resistance;
    if (model_) {
        effective *= model_->factor(model_->TNOM);
    }
    conductance_ = 1.0 / effective;

    nodes_.push_back(node_plus);
    nodes_.push_back(node_minus);
    if (node_plus)  node_plus->addConnectedDevice(this);
    if (node_minus) node_minus->addConnectedDevice(this);
}

// ---------------------------------------------------------------
// Setup
// ---------------------------------------------------------------

void Resistor::registerJacobianPattern(std::vector<int>& rows,
                                       std::vector<int>& cols) const {
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);

    // 2x2 conductance stamp in a FIXED order:
    //   0: [i][i]   1: [i][j]   2: [j][i]   3: [j][j]
    rows.push_back(ip); cols.push_back(ip);   // 0
    rows.push_back(ip); cols.push_back(in);   // 1
    rows.push_back(in); cols.push_back(ip);   // 2
    rows.push_back(in); cols.push_back(in);   // 3
}

// ---------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------

void Resistor::loadF(const DaeSolution& x, double time,
                     std::vector<double>& F) const {
    (void)time;
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);

    const double v = (ip >= 0 ? x.data()[ip] : 0.0)
                   - (in >= 0 ? x.data()[in] : 0.0);
    const double i = conductance_ * v;

    if (ip >= 0) F[ip] += i;
    if (in >= 0) F[in] -= i;
}

void Resistor::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                             double time) const {
    (void)x; (void)time;
    const double g = conductance_;

    emit(J, 0,  g);   // [i][i]
    emit(J, 1, -g);   // [i][j]
    emit(J, 2, -g);   // [j][i]
    emit(J, 3,  g);   // [j][j]
}

// ---------------------------------------------------------------
// Sweep hook (.DC R1 ...) / circuit temperature (.TEMP)
// ---------------------------------------------------------------

double Resistor::getEffectiveResistance() const {
    return (conductance_ > 0.0) ? 1.0 / conductance_ : 0.0;
}

void Resistor::setEffectiveResistance(double ohms) {
    if (!(ohms > 0.0)) ohms = 1e-12;   // near short instead of division by zero
    resistance_ = ohms;
    conductance_ = 1.0 / ohms;
}

void Resistor::setCircuitTemperature(double kelvin) {
    if (!model_) return;
    const double effective = resistance_ * model_->factor(kelvin);
    conductance_ = (effective > 0.0) ? 1.0 / effective : 0.0;
}

// ---------------------------------------------------------------
// Outputs
// ---------------------------------------------------------------

double Resistor::getCurrent() const {
    const double vp = nodes_[0] ? nodes_[0]->getVoltage() : 0.0;
    const double vn = nodes_[1] ? nodes_[1]->getVoltage() : 0.0;
    return conductance_ * (vp - vn);
}

double Resistor::getPower() const {
    return getVoltage() * getCurrent();
}

double Resistor::getVoltage() const {
    const double vp = nodes_[0] ? nodes_[0]->getVoltage() : 0.0;
    const double vn = nodes_[1] ? nodes_[1]->getVoltage() : 0.0;
    return vp - vn;
}

// ---------------------------------------------------------------
// Noise (Johnson-Nyquist): a Norton pair of current-noise density
//      S_i = 4kT/R   [A^2/Hz]
// across the resistor. White (frequency independent), so the sweep
// only has to multiply by the constant computed here.
// The effective (temperature scaled) resistance is the one the circuit
// actually has, not the nominal card value.
// ---------------------------------------------------------------

void Resistor::loadNoiseSources(const NoiseContext& ctx,
                                std::vector<NoiseSource>& out) const {
    if (!(conductance_ > 0.0)) return;
    NoiseSource s;
    s.label = "thermal";
    s.pos = rowOf(nodes_[0]);
    s.neg = rowOf(nodes_[1]);
    s.psd_a2 = 4.0 * ctx.boltzmann * ctx.temperature * conductance_;
    s.f_exponent = 0.0;
    out.push_back(s);
}

} // namespace e_XSpice

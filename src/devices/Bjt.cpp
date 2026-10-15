/**
 * @author Gabriel Henrique Silva
 * @file Bjt.cpp
 * @brief Bipolar junction transistor (NPN and PNP) — Ebers-Moll transport
 *        model in the shared DAE formulation.
 */

#include "Bjt.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"
#include "../models/core/DeviceConstants.h"
#include "../models/core_models/Junction.h"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

namespace e_XSpice {

using std::exp;
using std::min;
using std::max;

bool Bjt::debug_ = false;

// ===============================================================
// Construction
// ===============================================================

Bjt::Bjt(const std::string& name,
         Node* collector,
         Node* base,
         Node* emitter,
         RCP<const BjtModel> model,
         const BjtInstanceParams& inst)
    : Device(name, DeviceType::BJT)
    , model_(std::move(model))
    , inst_(inst) {
    if (!model_) {
        throw std::invalid_argument("BJT requires a model card");
    }
    if (inst_.AREA <= 0.0) {
        throw std::invalid_argument("BJT AREA must be positive");
    }

    p_ = model_->isPnp() ? -1 : 1;

    // Bind the thermal snapshot ONCE (setup-time): the card's, or the
    // instance's own when TEMP is overridden. Evaluation then reads plain
    // doubles through one pointer.
    if (inst_.TEMPGiven) {
        if (inst_.TEMP <= 0.0) {
            throw std::invalid_argument("BJT TEMP must be positive (Kelvin)");
        }
        ownThermal_ = model_->computeThermal(inst_.TEMP + inst_.DTEMP);
        thermal_ = &ownThermal_;
    } else {
        thermal_ = &model_->thermal();
    }

    nodes_.push_back(collector);   // 0
    nodes_.push_back(base);        // 1
    nodes_.push_back(emitter);     // 2
    if (collector) collector->addConnectedDevice(this);
    if (base)     base->addConnectedDevice(this);
    if (emitter)  emitter->addConnectedDevice(this);
}

Bjt::Bjt(const std::string& name,
         Node* collector,
         Node* base,
         Node* emitter)
    : Bjt(name, collector, base, emitter,
          BjtModel::defaultModel(), BjtInstanceParams()) {}

// ===============================================================
// Setup
// ===============================================================

int Bjt::registerStates(StateStore& store) {
    setStateIndex(store.registerStates(2));   // Q_BE, Q_BC
    return 2;
}

void Bjt::registerJacobianPattern(std::vector<int>& rows,
                                  std::vector<int>& cols) const {
    // Unknowns: b (base), c (collector), e (emitter) — the 3x3 block.
    const int rb = rowOf(nodes_[1]);
    const int rc = rowOf(nodes_[0]);
    const int re = rowOf(nodes_[2]);

    // Registration order (shared by loadJacobianF and loadJacobianQ):
    //   0: [b][b]  1: [b][c]  2: [b][e]
    //   3: [c][b]  4: [c][c]  5: [c][e]
    //   6: [e][b]  7: [e][c]  8: [e][e]
    rows.push_back(rb); cols.push_back(rb);  // 0
    rows.push_back(rb); cols.push_back(rc);  // 1
    rows.push_back(rb); cols.push_back(re);  // 2
    rows.push_back(rc); cols.push_back(rb);  // 3
    rows.push_back(rc); cols.push_back(rc);  // 4
    rows.push_back(rc); cols.push_back(re);  // 5
    rows.push_back(re); cols.push_back(rb);  // 6
    rows.push_back(re); cols.push_back(rc);  // 7
    rows.push_back(re); cols.push_back(re);  // 8
}

void Bjt::beginSolve(bool initial_solve) {
    if (initial_solve) {
        // Start the first Newton solve from the most linear point of each
        // junction (Xyce initJctFlag_ behaviour).
        vbe_prev_ = thermal_->beVcrit;
        vbc_prev_ = 0.0;
        first_solve_ = true;
    }
}

// ===============================================================
// Junction evaluation (the transport/Ebers-Moll laws)
// ===============================================================

double Bjt::pnjlim(double vnew, double vold, double vt, double vcrit) {
    if (vnew > vcrit && std::abs(vnew - vold) > 2.0 * vt) {
        if (vold > 0.0) {
            const double arg = 1.0 + (vnew - vold) / vt;
            if (arg > 0.0) {
                vnew = vold + vt * log(arg);
            } else {
                vnew = vcrit;
            }
        } else {
            vnew = vt * log(vnew / vt);
        }
    }
    return vnew;
}

void Bjt::limitedJunctionVoltages(const DaeSolution& x,
                                 bool allow_limiting) const {
    const int rc = rowOf(nodes_[0]);
    const int rb = rowOf(nodes_[1]);
    const int re = rowOf(nodes_[2]);

    const double vc = (rc >= 0) ? x.data()[rc] : 0.0;
    const double vb = (rb >= 0) ? x.data()[rb] : 0.0;
    const double ve = (re >= 0) ? x.data()[re] : 0.0;

    // Raw junction voltages, NPN-oriented.
    // NPN-ORIENTED junction voltages: the polarity multiply is what makes
    // the limiters and the card laws polarity-agnostic (for a PNP the
    // physical forward direction is a NEGATIVE vb − ve).
    double vbe = p_ * (vb - ve);
    double vbc = p_ * (vb - vc);

    limiting_active_ = false;
    if (allow_limiting && !first_solve_) {
        const double nvbe = pnjlim(vbe, vbe_prev_, thermal_->vte,
                                   thermal_->beVcrit);
        if (nvbe != vbe) { vbe = nvbe; limiting_active_ = true; }
        const double nvbc = pnjlim(vbc, vbc_prev_, thermal_->vtc,
                                   thermal_->bcVcrit);
        if (nvbc != vbc) { vbc = nvbc; limiting_active_ = true; }
    }
    vbe_lim_ = vbe;   // NPN-oriented (already polarity-applied)
    vbc_lim_ = vbc;

    if (debug_) {
        std::fprintf(stderr,
                     "[BJ] raw vbe=%.8g vbc=%.8g -> lim vbe=%.8g vbc=%.8g "
                     "(prev=%.6g/%.6g, lim=%d)\n",
                     vb - ve, vb - vc, vbe, vbc, vbe_prev_, vbc_prev_,
                     static_cast<int>(limiting_active_));
    }
}

void Bjt::evaluate(double vbe, double vbc) const {
    const BjtThermal& t = *thermal_;
    const double area = inst_.AREA;

    const double isBE = t.beIsat * area;   // saturation currents, area-scaled
    const double isBC = t.bcIsat * area;
    const double betaF = t.betaf * area;
    const double betaR = t.betar * area;
    const double gmin = t.gmin;

    // ---- junction diode currents (nonideal base components) ----
    double arg = vbe / t.vte;
    arg = min(Const::MaxExpArg, arg);
    const double evbe = exp(arg);
    const double ibe = betaF * (evbe - 1.0) + gmin * vbe;
    const double gib = betaF * evbe / t.vte + gmin;

    arg = vbc / t.vtc;
    arg = min(Const::MaxExpArg, arg);
    const double evbc = exp(arg);
    const double ibc = betaR * (evbc - 1.0) + gmin * vbc;
    const double gibc = betaR * evbc / t.vtc + gmin;

    // ---- base-width modulation (Early effect) ----
    // K_E = 1 − vbc/VAF − vbe/VAR. The floor keeps ict and its derivative
    // finite deep in saturation (where the SPICE expression turns negative);
    // the derivative stays the linear one, so Newton still walks back out.
    const double dke_dvbe = -1.0 / model_->VAR;   // Vbe modulates with VAR
    const double dke_dvbc = -1.0 / model_->VAF;   // Vbc modulates with VAF
    double ke = 1.0 + dke_dvbe * vbe + dke_dvbc * vbc;
    const double ke_floor = 1e-6;
    if (ke < ke_floor) ke = ke_floor;

    // ---- transport current ----
    const double ict = (isBE * evbe - isBC * evbc) * ke;
    const double gm = isBE * evbe / t.vte * ke +
                      (isBE * evbe - isBC * evbc) * dke_dvbe;
    const double go = -isBC * evbc / t.vtc * ke +
                      (isBE * evbe - isBC * evbc) * dke_dvbc;

    // ---- terminal currents (positive INTO the transistor) ----
    ib_ = ibe + ibc;
    ic_ = ict - ibc;
    ie_ = -ib_ - ic_;

    gib_ = gib;
    gibc_ = gibc;
    gm_ = gm;
    go_ = go;

    // ---- charges ----
    // Depletion parts come from the shared junction model (tt = 0 — the
    // transit part is separate below); the capacitances scale with AREA
    // only, and the transit times with AREA as well (SPICE scaling).
    const JunctionGeom geomBE{t.beCap.cjo * area, 0.0, model_->MJE};
    const JunctionGeom geomBC{t.bcCap.cjo * area, 0.0, model_->MJC};
    const JunctionChargePoint depBE = junctionCharge(vbe, 0.0, 0.0, geomBE, t.beCap);
    const JunctionChargePoint depBC = junctionCharge(vbc, 0.0, 0.0, geomBC, t.bcCap);

    const double tf = t.tf * area;   // transit times, area-scaled
    const double tr = t.tr * area;

    qbe_ = depBE.q + tf * ict;
    qbc_ = depBC.q - tr * ict;

    cbe_v_ = depBE.c + tf * gm;
    cbe_w_ = tf * go;
    cbc_u_ = -tr * gm;
    cbc_v_ = depBC.c - tr * go;
}

// ===============================================================
// Loaders
// ===============================================================

void Bjt::loadF(const DaeSolution& x, double time,
                std::vector<double>& F) const {
    (void)time;
    const int rc = rowOf(nodes_[0]);
    const int rb = rowOf(nodes_[1]);
    const int re = rowOf(nodes_[2]);

    limitedJunctionVoltages(x, true);
    evaluate(vbe_lim_, vbc_lim_);
    // Report the TERMINAL voltages (what a scope would measure).
    vbe_cached_ = p_ * vbe_lim_;
    vbc_cached_ = p_ * vbc_lim_;

    if (rb >= 0) F[rb] += p_ * ib_;
    if (rc >= 0) F[rc] += p_ * ic_;
    if (re >= 0) F[re] += p_ * ie_;

    vbe_prev_ = vbe_lim_;
    vbc_prev_ = vbc_lim_;
    first_solve_ = false;
}

void Bjt::applySolutionLimit(std::vector<double>& x) const {
    if (!limiting_active_) return;

    // Rebuild the unknowns from the LIMITED junction voltages. The emitter
    // stays as-is (it is the anchor); p_*vje converts the NPN-oriented
    // limited voltage back to the terminal difference (1/p == p for ±1).
    const int rc = rowOf(nodes_[0]);
    const int rb = rowOf(nodes_[1]);
    const int re = rowOf(nodes_[2]);

    const double ve = (re >= 0) ? x[re] : 0.0;
    if (rb >= 0) x[rb] = ve + p_ * vbe_lim_;
    // A grounded base has no row; anchor the collector on it instead.
    const double vb = (rb >= 0) ? x[rb] : 0.0;
    if (rc >= 0) x[rc] = vb - p_ * vbc_lim_;
}

void Bjt::loadQ(const DaeSolution& x, double time,
                std::vector<double>& Q) const {
    (void)time;
    const int rc = rowOf(nodes_[0]);
    const int rb = rowOf(nodes_[1]);
    const int re = rowOf(nodes_[2]);

    limitedJunctionVoltages(x, true);
    evaluate(vbe_lim_, vbc_lim_);

    if (rb >= 0) Q[rb] += p_ * (qbe_ + qbc_);
    if (re >= 0) Q[re] -= p_ * qbe_;
    if (rc >= 0) Q[rc] -= p_ * qbc_;
}

void Bjt::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                        double time) const {
    (void)x; (void)time;
    emit(J, 0,  gib_ + gibc_);            // [b][b]
    emit(J, 1, -gibc_);                   // [b][c]
    emit(J, 2, -gib_);                    // [b][e]
    emit(J, 3,  gm_ + go_ - gibc_);       // [c][b]
    emit(J, 4, -go_ + gibc_);             // [c][c]
    emit(J, 5, -gm_);                     // [c][e]
    emit(J, 6, -(gib_ + gm_ + go_));      // [e][b]
    emit(J, 7,  go_);                     // [e][c]
    emit(J, 8,  gib_ + gm_);              // [e][e]
}

void Bjt::loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                        double time) const {
    (void)x; (void)time;
    emit(J, 0,  cbe_v_ + cbe_w_ + cbc_u_ + cbc_v_);  // [b][b]
    emit(J, 1, -cbe_w_ - cbc_v_);                    // [b][c]
    emit(J, 2, -cbe_v_ - cbc_u_);                    // [b][e]
    emit(J, 3, -cbc_u_ - cbc_v_);                    // [c][b]
    emit(J, 4,  cbc_v_);                             // [c][c]
    emit(J, 5,  cbc_u_);                             // [c][e]
    emit(J, 6, -cbe_v_ - cbe_w_);                    // [e][b]
    emit(J, 7,  cbe_w_);                             // [e][c]
    emit(J, 8,  cbe_v_);                             // [e][e]
}

void Bjt::storeState(const DaeSolution& x, StateStore& store) const {
    limitedJunctionVoltages(x, true);
    evaluate(vbe_lim_, vbc_lim_);
    const double q[2] = {qbe_, qbc_};
    store.setDeviceQ(state_index_, q, 2);
}

void Bjt::loadChargeHistory(const StateStore& store, double q_coeff,
                            double qdot_coeff,
                            std::vector<double>& rhs) const {
    if (q_coeff == 0.0) return;

    const int rc = rowOf(nodes_[0]);
    const int rb = rowOf(nodes_[1]);
    const int re = rowOf(nodes_[2]);

    const double qbe = q_coeff * store.deviceQ(state_index_, 0) +
                       qdot_coeff * store.deviceQDot(state_index_, 0);
    const double qbc = q_coeff * store.deviceQ(state_index_, 1) +
                       qdot_coeff * store.deviceQDot(state_index_, 1);

    if (rb >= 0) rhs[rb] += p_ * (qbe + qbc);
    if (re >= 0) rhs[re] -= p_ * qbe;
    if (rc >= 0) rhs[rc] -= p_ * qbc;
}

void Bjt::refreshOutputs(const DaeSolution& x, double time) const {
    (void)time;
    // Final cache refresh from the converged solution (no limiting — the
    // raw voltages are the answer). Needed by distributed simulation, where
    // only the owner rank evaluates the loaders.
    const int rc = rowOf(nodes_[0]);
    const int rb = rowOf(nodes_[1]);
    const int re = rowOf(nodes_[2]);

    const double vc = (rc >= 0) ? x.data()[rc] : 0.0;
    const double vb = (rb >= 0) ? x.data()[rb] : 0.0;
    const double ve = (re >= 0) ? x.data()[re] : 0.0;

    vbe_cached_ = vb - ve;   // terminal voltages (report values)
    vbc_cached_ = vb - vc;
    vbe_lim_ = p_ * vbe_cached_;   // NPN-oriented (what evaluate() wants)
    vbc_lim_ = p_ * vbc_cached_;
    evaluate(vbe_lim_, vbc_lim_);
}

bool Bjt::acceptSolution() const {
    // Newton iteration is only accepted when no limiter intervened
    // (Xyce: origFlag = false → isConverged() == false).
    return !limiting_active_;
}

void Bjt::confirmAccept() {
    limiting_active_ = false;
}

// ===============================================================
// Outputs
// ===============================================================

double Bjt::getPower() const {
    // Σ v_terminal·i_terminal, positive = the device absorbs power.
    const double vb = nodes_[1] ? nodes_[1]->getVoltage() : 0.0;
    const double vc = nodes_[0] ? nodes_[0]->getVoltage() : 0.0;
    const double ve = nodes_[2] ? nodes_[2]->getVoltage() : 0.0;
    return vb * (p_ * ib_) + vc * (p_ * ic_) + ve * (p_ * ie_);
}

} // namespace e_XSpice

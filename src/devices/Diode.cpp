/**
 * @author Gabriel Henrique Silva
 * @file Diode.cpp
 * @brief PN-junction diode — charge-based, robust formulation (Xyce level 1).
 */

#include "Diode.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"
#include "../models/core/DeviceConstants.h"
#include "../models/core_models/Junction.h"
#include <stdexcept>
#include <algorithm>

namespace e_XSpice {

using std::exp;
using std::log;
using std::min;
using std::max;

bool Diode::debug_ = false;

Diode::Diode(const std::string& name,
             Node* anode,
             Node* cathode,
             RCP<const DiodeModel> model,
             const DiodeInstanceParams& inst)
    : Device(name, DeviceType::DIODE)
    , model_(std::move(model))
    , inst_(inst)
{
    if (!model_) {
        throw std::invalid_argument("Diode requires a model card");
    }
    if (inst_.AREA <= 0.0) {
        throw std::invalid_argument("Diode AREA must be positive");
    }
    if (inst_.M <= 0.0) {
        throw std::invalid_argument("Diode multiplicity M must be positive");
    }

    // Bind the thermal snapshot once (setup-time): the model's, or the
    // instance's own when TEMP is overridden. Evaluation then reads plain
    // doubles through one pointer — no maps, no recomputation in loops.
    if (inst_.TEMPGiven) {
        if (inst_.TEMP <= 0.0) {
            throw std::invalid_argument("Diode TEMP must be positive (Kelvin)");
        }
        ownThermal_ = model_->computeThermal(inst_.TEMP + inst_.DTEMP);
        thermal_ = &ownThermal_;
    } else {
        thermal_ = &model_->thermal();
    }

    nodes_.push_back(anode);
    nodes_.push_back(cathode);
    if (anode)   anode->addConnectedDevice(this);
    if (cathode) cathode->addConnectedDevice(this);
}

Diode::Diode(const std::string& name,
             Node* anode,
             Node* cathode)
    : Diode(name, anode, cathode,
            DiodeModel::defaultModel(), DiodeInstanceParams()) {
}

// ===============================================================
// Setup
// ===============================================================

int Diode::registerStates(StateStore& store) {
    setStateIndex(store.registerStates(1));
    return 1;
}

void Diode::registerJacobianPattern(std::vector<int>& rows,
                                    std::vector<int>& cols) const {
    // Unknowns: 0 = anode(ext)  1 = cathode(ext)  2 = Pri (internal, or
    // aliased to anode when RS = 0).
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const int ir = priRow();

    //   0: [p][p]   1: [p][pri]
    //   2: [n][n]   3: [n][pri]
    //   4: [pri][p] 5: [pri][n] 6: [pri][pri]
    rows.push_back(ip); cols.push_back(ip);  // 0
    rows.push_back(ip); cols.push_back(ir);  // 1
    rows.push_back(in); cols.push_back(in);  // 2
    rows.push_back(in); cols.push_back(ir);  // 3
    rows.push_back(ir); cols.push_back(ip);  // 4
    rows.push_back(ir); cols.push_back(in);  // 5
    rows.push_back(ir); cols.push_back(ir);  // 6
}

void Diode::beginSolve(bool initial_solve) {
    if (initial_solve) {
        // Xyce initJctFlag_: start the first Newton solve from the critical
        // voltage (the most linear point of the exponential region).
        vd_prev_ = thermal_->tVcrit;
        first_solve_ = true;
    }
}

// ===============================================================
// Junction evaluation
// ===============================================================

void Diode::evaluateJunction(double vd) const {
    const DiodeThermal& t = *thermal_;
    const double Area = inst_.AREA;
    const double mult = inst_.M;
    const double Isat  = t.isat * Area;
    const double IsatR = t.isatR * Area;
    const double gmin  = model_->GMIN;

    const double vte    = t.vte;
    const double vteR   = t.vteR;
    const double vteBRK = t.vteBRK;

    double Id = 0.0, Gd = 0.0;

    if (vd >= -3.0 * vte) {
        // ---- forward / normal exponential region ----
        double arg1 = vd / vte;
        arg1 = min(Const::MaxExpArg, arg1);
        const double evd = exp(arg1);

        const double Inorm = Isat * (evd - 1.0) + gmin * vd;
        const double Gd1    = Isat * evd / vte + gmin;

        double arg2 = vd / vteR;
        arg2 = min(Const::MaxExpArg, arg2);
        const double evdR = exp(arg2);
        const double Irec = IsatR * (evdR - 1.0);
        const double Gd2  = IsatR * evdR / vteR;

        // high-injection knee: Khi = sqrt(IKF/(IKF+I_norm))
        double Khi = 1.0, DKhi = 0.0;
        if (t.ikf > 0.0 && (t.ikf + Inorm) > 0.0) {
            Khi = sqrt(t.ikf / (t.ikf + Inorm));
            DKhi = -0.5 * Khi * Gd1 / (t.ikf + Inorm);
        }

        // generation-depletion factor for the recombination current
        double Kgen = 0.0, DKgen = 0.0;
        if (Irec != 0.0 && t.tJctPot > 0.0) {
            const double u = (1.0 - vd / t.tJctPot) * (1.0 - vd / t.tJctPot) + 0.005;
            Kgen = exp(0.5 * model_->M * log(u));              // u^(M/2)
            const double du = 2.0 * (1.0 - vd / t.tJctPot);
            DKgen = -0.5 * model_->M * du * Kgen / u / t.tJctPot;
        }

        Id = Inorm * Khi + Irec * Kgen;
        Gd = Gd1 * Khi + Inorm * DKhi + Gd2 * Kgen + Irec * DKgen;
    } else if (!model_->BVGiven || vd >= -t.tBrkdwnV) {
        // ---- linear reverse-bias region (SPICE3 third-power law) ----
        double arg = 3.0 * vte / (vd * Const::E);
        arg = arg * arg * arg;
        Id = -Isat * (1.0 + arg) + gmin * vd;
        Gd = Isat * 3.0 * arg / vd + gmin;
    } else {
        // ---- reverse breakdown ----
        double arg1 = -(t.tBrkdwnV + vd) / vteBRK;
        arg1 = min(Const::MaxExpArg, arg1);
        const double evrev = exp(arg1);
        Id = -Isat * evrev + gmin * vd;
        Gd = Isat * evrev / vteBRK + gmin;
    }

    // ---------------- junction charge ----------------
    // Shared with the BJT junctions (Junction.h): the SPICE3f5 construction
    // where Q and C are both C0-continuous at V = FC·VJ. The TT transit term
    // no longer depends on CJO != 0 (SPICE semantics).
    const JunctionGeom geom{t.tJctCap * Area, model_->TT, model_->M};
    const JunctionCap  cap{t.tJctPot, t.tJctCap, t.tDepCap,
                           t.F2, t.F3, t.tF1};
    const JunctionChargePoint jc = junctionCharge(vd, Id, Gd, geom, cap);
    // Multiplicity scales the whole junction (parallel devices).
    vd_limited_ = vd;
    id_limited_ = Id * mult;
    gd_limited_ = Gd * mult;
    qd_limited_ = jc.q * mult;
    cd_limited_ = jc.c * mult;

    // Remember the limited voltage for the next Newton iteration
    // (Xyce keeps Vd_old in its store vector for exactly this purpose).
    vd_prev_ = vd;
}

double Diode::pnjlim(double vnew, double vold, double vt, double vcrit) {
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

double Diode::limitedJunctionVoltage(const DaeSolution& x) const {
    const int ir = priRow();
    const int in = rowOf(nodes_[1]);
    const int ip = rowOf(nodes_[0]);

    const double vpp = (ir >= 0) ? x.data()[ir] : 0.0;
    const double vn  = (in >= 0) ? x.data()[in] : 0.0;
    const double vp  = (ip >= 0) ? x.data()[ip] : 0.0;
    (void)vp;

    double vd = vpp - vn;

    const double vd_raw = vd;
    const DiodeThermal& t = *thermal_;
    limiting_active_ = false;
    if (!first_solve_) {
        if (model_->BVGiven && vd < min(0.0, -model_->BV + 10.0 * t.vteBRK)) {
            const double vdtmp = -(model_->BV + vd);
            const double vdtmp_old = -(vd_prev_ + model_->BV);
            const double vnew = pnjlim(vdtmp, vdtmp_old, t.vteBRK, t.tVcrit);
            if (vnew != vdtmp) { vd = -(vnew + model_->BV); limiting_active_ = true; }
        } else {
            const double vnew = pnjlim(vd, vd_prev_, t.vte, t.tVcrit);
            if (vnew != vd) { vd = vnew; limiting_active_ = true; }
        }
    }
    if (debug_) {
        std::fprintf(stderr, "[D1] vd raw=%.8g -> limited=%.8g (prev=%.8g, limited=%d)\n",
                     vd_raw, vd, vd_prev_, static_cast<int>(limiting_active_));
    }
    return vd;
}

double Diode::junctionVoltageFromNodes() const {
    const double va = nodes_[0] ? nodes_[0]->getVoltage() : 0.0;
    const double vc = nodes_[1] ? nodes_[1]->getVoltage() : 0.0;
    return va - vc;
}

// ===============================================================
// Evaluation
// ===============================================================

void Diode::loadF(const DaeSolution& x, double time,
                 std::vector<double>& F) const {
    (void)time;
    const int ip = rowOf(nodes_[0]);
    const int in = rowOf(nodes_[1]);
    const int ir = priRow();

    const double vd = limitedJunctionVoltage(x);
    evaluateJunction(vd);

    const double v_pri = (ir >= 0) ? x.data()[ir] : 0.0;
    const double v_an  = (ip >= 0) ? x.data()[ip] : 0.0;
    const double Ir = thermal_->gspr * (v_an - v_pri);
    const double Id = id_limited_;

    if (ip >= 0) F[ip] += Ir;      // current entering anode = resistor current
    if (in >= 0) F[in] -= Id;      // current leaving cathode
    // Internal node KCL (aliased to the anode row when RS = 0 — that is
    // where the junction current enters).
    if (ir >= 0) F[ir] += (Id - Ir);

    first_solve_ = false;
}

void Diode::applySolutionLimit(std::vector<double>& x) const {
    if (!limiting_active_) return;

    // SPICE-style limiting write-back: the next Newton iteration restarts
    // from the limited junction voltage, anchored on the cathode (the
    // cathode row is never the limited one). The anode side is shifted by
    // the SAME delta on both the Pri row and the external anode row, so a
    // series RS carries exactly its evaluated current into the next
    // iteration (no artificial RS residual — see Bjt::applySolutionLimit).
    // With RS = 0 the Pri row IS the anode row, so a single write covers
    // both cases; a grounded anode has no row at all (ir < 0).
    const int in = rowOf(nodes_[1]);
    const int ir = priRow();
    const int ip = rowOf(nodes_[0]);
    if (ir >= 0) {
        const double updated = (in >= 0) ? x[in] + vd_limited_ : vd_limited_;
        const double delta = updated - x[ir];
        x[ir] = updated;
        if (ip >= 0 && ip != ir) x[ip] += delta;
    }
}

void Diode::loadQ(const DaeSolution& x, double time,
                 std::vector<double>& Q) const {
    (void)time;
    const int in = rowOf(nodes_[1]);
    const int ir = priRow();

    const double vd = limitedJunctionVoltage(x);
    evaluateJunction(vd);

    if (in >= 0) Q[in] -= qd_limited_;
    // +qd goes to the Pri row; when RS=0 Pri aliases the anode row, and when
    // the anode is ground (e.g. "D2 0 5 ...") there is no row at all (ir<0).
    if (ir >= 0 && ir != in) Q[ir] += qd_limited_;
}

void Diode::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                         double time) const {
    (void)x; (void)time;
    const double Gd = gd_limited_;
    const double Gspr = thermal_->gspr;

    emit(J, 0,  Gspr);        // [p][p]
    emit(J, 1, -Gspr);        // [p][pri]
    emit(J, 2,  Gd);          // [n][n]
    emit(J, 3, -Gd);          // [n][pri]
    emit(J, 4, -Gspr);        // [pri][p]
    emit(J, 5, -Gd);          // [pri][n]
    emit(J, 6,  Gspr + Gd);   // [pri][pri]
}

void Diode::loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                         double time) const {
    (void)x; (void)time;
    const double Cd = cd_limited_;

    emit(J, 2,  Cd);   // [n][n]
    emit(J, 3, -Cd);   // [n][pri]
    emit(J, 5, -Cd);   // [pri][n]
    emit(J, 6,  Cd);   // [pri][pri]
}

void Diode::storeState(const DaeSolution& x, StateStore& store) const {
    const double vd = limitedJunctionVoltage(x);
    evaluateJunction(vd);
    store.setDeviceQ(state_index_, &qd_limited_, 1);
}

void Diode::loadChargeHistory(const StateStore& store, double q_coeff,
                             double qdot_coeff,
                             std::vector<double>& rhs) const {
    if (q_coeff == 0.0) return;

    const int in = rowOf(nodes_[1]);
    const int ir = priRow();
    const double qdn  = store.deviceQ(state_index_, 0);       // Q_j(t_n)
    const double qdn_dot = qdot_coeff != 0.0
                           ? store.deviceQDot(state_index_, 0)
                           : 0.0;
    const double qterm = q_coeff * qdn + qdn_dot;

    if (in >= 0) rhs[in] -= qterm;
    if (ir != in && ir >= 0) rhs[ir] += qterm;   // Pri==anode when RS=0
}

bool Diode::acceptSolution() const {
    // Newton iteration is only accepted when no limiter intervened
    // (Xyce: origFlag = false → isConverged() == false).
    return !limiting_active_;
}

void Diode::resetIterationState(const DaeSolution& x) {
    // Re-anchor the limiter history to the RAW junction voltage at the
    // restored point (no evaluation: the next assembly re-evaluates before
    // using any cache). pnjlim(raw, raw) cannot misfire.
    limiting_active_ = false;
    const int ir = priRow();
    const int in = rowOf(nodes_[1]);
    const double vpp = (ir >= 0) ? x.data()[ir] : 0.0;
    const double vn = (in >= 0) ? x.data()[in] : 0.0;
    vd_prev_ = vpp - vn;
}

void Diode::confirmAccept() {
    limiting_active_ = false;
}

// ===============================================================
// Outputs
// ===============================================================

double Diode::currentForVoltage(double vd) const {
    // Evaluates ONLY the conduction current (no cache mutation), used by
    // getCurrent()/output paths.
    const DiodeThermal& t = *thermal_;
    const double Area = inst_.AREA;
    const double mult = inst_.M;
    const double Isat = t.isat * Area;
    const double IsatR = t.isatR * Area;
    const double gmin = model_->GMIN;

    if (vd >= -3.0 * t.vte) {
        double arg1 = vd / t.vte;
        arg1 = std::min(100.0, arg1);
        const double evd = std::exp(arg1);
        const double Inorm = Isat * (evd - 1.0) + gmin * vd;

        double arg2 = vd / t.vteR;
        arg2 = std::min(100.0, arg2);
        const double evdR = std::exp(arg2);
        const double Irec = IsatR * (evdR - 1.0);

        double Khi = 1.0;
        if (t.ikf > 0.0 && (t.ikf + Inorm) > 0.0)
            Khi = std::sqrt(t.ikf / (t.ikf + Inorm));

        double Kgen = 0.0;
        if (Irec != 0.0 && t.tJctPot > 0.0) {
            const double u = (1.0 - vd / t.tJctPot) * (1.0 - vd / t.tJctPot) + 0.005;
            Kgen = std::exp(0.5 * model_->M * std::log(u));
        }
        return (Inorm * Khi + Irec * Kgen) * mult;
    }
    if (!model_->BVGiven || vd >= -t.tBrkdwnV) {
        double arg = 3.0 * t.vte / (vd * Const::E);
        arg = arg * arg * arg;
        return (-Isat * (1.0 + arg) + gmin * vd) * mult;
    }
    double arg1 = -(t.tBrkdwnV + vd) / t.vteBRK;
    arg1 = std::min(100.0, arg1);
    return (-Isat * std::exp(arg1) + gmin * vd) * mult;
}

double Diode::getCurrent() const {
    return currentForVoltage(junctionVoltageFromNodes());
}

double Diode::getPower() const {
    return getVoltage() * getCurrent();
}

} // namespace e_XSpice

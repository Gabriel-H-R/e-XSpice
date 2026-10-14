/**
 * @author Gabriel Henrique Silva
 * @file Diode.cpp
 * @brief PN-junction diode — charge-based, robust formulation (Xyce level 1).
 */

#include "Diode.h"
#include "../core/DaeSolution.h"
#include "../solver/DaeMatrix.h"
#include <stdexcept>
#include <algorithm>

namespace e_XSpice {

// Xyce constants (N_DEV_Const.h)
static constexpr double CONSTQ      = 1.6021918e-19;
static constexpr double CONSTboltz  = 1.3806226e-23;
static constexpr double CONSTKoverQ = CONSTboltz / CONSTQ;   // 8.6170344e-05 V/K
static constexpr double CONSTREFTEMP = 300.15;
static constexpr double CONSTEg0    = 1.16;
static constexpr double CONSTalphaEg = 7.02e-4;
static constexpr double CONSTbetaEg  = 1108.0;
static constexpr double CONSTEg300  = 1.1150877;
static constexpr double CONSTe      = 2.71828182845904523536028747135;  // e
static constexpr double CONSTroot2  = 1.41421356237309504880168872421;
static constexpr double CONSTMAX_EXP_ARG = 100.0;

using std::exp;
using std::log;
using std::min;
using std::max;

// ===============================================================
// Construction / temperature scaling
// ===============================================================

bool Diode::debug_ = false;

Diode::Diode(const std::string& name,
             Node* anode,
             Node* cathode,
             const ModelParams& params)
    : Device(name, DeviceType::DIODE)
    , params_(params)
{
    if (params_.IS <= 0.0)  throw std::invalid_argument("Diode IS must be positive");
    if (params_.N  <= 0.0)  throw std::invalid_argument("Diode N must be positive");
    if (params_.CJO < 0.0)  throw std::invalid_argument("Diode CJO must be non-negative");
    if (params_.VJ <= 0.0)  throw std::invalid_argument("Diode VJ must be positive");

    // History-dependent parameter guards (Xyce Model::processParams)
    params_.M  = min(params_.M, 0.9);
    params_.EG = max(params_.EG, 0.1);
    params_.FC = min(params_.FC, 0.95);

    nodes_.push_back(anode);
    nodes_.push_back(cathode);
    if (anode)   anode->addConnectedDevice(this);
    if (cathode) cathode->addConnectedDevice(this);

    updateTemperature();
    if (params_.BVGiven) fitBreakdownVoltage();
    // PNJLIM critical voltage: vcrit = vte·ln(vte/(√2·IS))
    tVcrit_ = vte_ * log(vte_ / (CONSTroot2 * isat_));
}

void Diode::updateTemperature() {
    const double Temp  = params_.TEMP;
    const double TNOM  = params_.TNOM;

    const double vtnom = CONSTKoverQ * TNOM;
    const double vt    = CONSTKoverQ * Temp;
    vt_   = vt;
    vte_  = params_.N  * vt;
    vteR_ = params_.NR * vt;
    vteBRK_ = (params_.NBV > 0.0 ? params_.NBV : 1.0) * vt;

    const double egfet  = CONSTEg0 - (CONSTalphaEg * Temp * Temp) / (Temp + CONSTbetaEg);
    const double egfet1 = CONSTEg0 - (CONSTalphaEg * TNOM * TNOM) / (TNOM + CONSTbetaEg);
    const double arg = -egfet / (2.0 * CONSTboltz * Temp)
                     + CONSTEg300 / (CONSTboltz * (CONSTREFTEMP + CONSTREFTEMP));
    const double arg1 = -egfet1 / (2.0 * CONSTboltz * TNOM)
                      + CONSTEg300 / (2.0 * CONSTboltz * CONSTREFTEMP);
    const double fact2 = Temp / CONSTREFTEMP;
    const double fact1 = TNOM / CONSTREFTEMP;
    const double pbfact  = -2.0 * vt  * (1.5 * log(fact2) + CONSTQ * arg);
    const double pbfact1 = -2.0 * vtnom * (1.5 * log(fact1) + CONSTQ * arg1);

    // junction potential at temperature
    const double pbo = (params_.VJ - pbfact1) / fact1;
    const double gmaold = (params_.VJ - pbo) / pbo;
    tJctPot_ = pbfact + fact2 * pbo;
    const double gmanew = (tJctPot_ - pbo) / pbo;

    // zero-bias capacitance at temperature
    tJctCap_ = params_.CJO / (1.0 + params_.M * (4.0e-4 * (TNOM - CONSTREFTEMP) - gmaold));
    tJctCap_ *= 1.0 + params_.M * (4.0e-4 * (Temp - CONSTREFTEMP) - gmanew);

    // saturation currents at temperature (SPICE3 formula)
    isat_  = params_.IS  * exp(((Temp / TNOM) - 1.0) * params_.EG / (params_.N  * vt)
                             + params_.XTI / params_.N  * log(Temp / TNOM));
    isatR_ = params_.ISR * exp(((Temp / TNOM) - 1.0) * params_.EG / (params_.NR * vt)
                             + params_.XTI / params_.NR * log(Temp / TNOM));

    tDepCap_ = params_.FC * tJctPot_;
    const double xfc = log(1.0 - params_.FC);
    F2_ = exp((1.0 + params_.M) * xfc);          // (1−FC)^(1+M)
    F3_ = 1.0 - params_.FC * (1.0 + params_.M);
    tF1_ = tJctPot_ * (1.0 - exp((1.0 - params_.M) * xfc)) / (1.0 - params_.M);

    gspr_ = (params_.RS > 0.0) ? 1.0 / params_.RS : 0.0;
    ikf_  = params_.IKF;
}

void Diode::fitBreakdownVoltage() {
    // Fixed-point fit of tBrkdwnV so that I_D = −IBV at Vd = −BV (Xyce)
    const double Temp = params_.TEMP, TNOM = params_.TNOM;
    const double tempBV = params_.BV * (1.0 + (Temp - TNOM) * (0.0 + 0.0 * (Temp - TNOM)));
    const double reltol = 1.0e-3;
    const double vt = vt_;

    double cbv = params_.IBV;
    double xbv = 0.0;

    if (cbv < isat_ * tempBV / vt) {
        xbv = tempBV;
    } else {
        const double tol = reltol * cbv;
        xbv = tempBV - params_.NBV * vt * log(1.0 + cbv / isat_);
        for (int i = 0; i < 25; ++i) {
            xbv = tempBV - params_.NBV * vt * log(cbv / isat_ + 1.0 - xbv / vt);
            const double xcbv = isat_ * (exp((tempBV - xbv) / (params_.NBV * vt)) - 1.0 + xbv / vt);
            if (std::abs(xcbv - cbv) <= tol) break;
        }
    }
    tBrkdwnV_ = xbv;
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
        vd_prev_ = tVcrit_;
        first_solve_ = true;
    }
}

// ===============================================================
// Junction evaluation
// ===============================================================

void Diode::evaluateJunction(double vd) const {
    const double Area = params_.AREA;
    const double Isat  = isat_ * Area;
    const double IsatR = isatR_ * Area;
    const double gmin  = params_.GMIN;

    const double vte   = vte_;
    const double vteR  = vteR_;
    const double vteBRK = vteBRK_;

    double Id = 0.0, Gd = 0.0;

    if (vd >= -3.0 * vte) {
        // ---- forward / normal exponential region ----
        double arg1 = vd / vte;
        arg1 = min(CONSTMAX_EXP_ARG, arg1);
        const double evd = exp(arg1);

        const double Inorm = Isat * (evd - 1.0) + gmin * vd;
        const double Gd1    = Isat * evd / vte + gmin;

        double arg2 = vd / vteR;
        arg2 = min(CONSTMAX_EXP_ARG, arg2);
        const double evdR = exp(arg2);
        const double Irec = IsatR * (evdR - 1.0);
        const double Gd2  = IsatR * evdR / vteR;

        // high-injection knee: Khi = sqrt(IKF/(IKF+I_norm))
        double Khi = 1.0, DKhi = 0.0;
        if (ikf_ > 0.0 && (ikf_ + Inorm) > 0.0) {
            Khi = sqrt(ikf_ / (ikf_ + Inorm));
            DKhi = -0.5 * Khi * Gd1 / (ikf_ + Inorm);
        }

        // generation-depletion factor for the recombination current
        double Kgen = 0.0, DKgen = 0.0;
        if (Irec != 0.0 && tJctPot_ > 0.0) {
            const double u = (1.0 - vd / tJctPot_) * (1.0 - vd / tJctPot_) + 0.005;
            Kgen = exp(0.5 * params_.M * log(u));              // u^(M/2)
            const double du = 2.0 * (1.0 - vd / tJctPot_);
            DKgen = -0.5 * params_.M * du * Kgen / u / tJctPot_;
        }

        Id = Inorm * Khi + Irec * Kgen;
        Gd = Gd1 * Khi + Inorm * DKhi + Gd2 * Kgen + Irec * DKgen;
    } else if (!params_.BVGiven || vd >= -tBrkdwnV_) {
        // ---- linear reverse-bias region (SPICE3 third-power law) ----
        double arg = 3.0 * vte / (vd * CONSTe);
        arg = arg * arg * arg;
        Id = -Isat * (1.0 + arg) + gmin * vd;
        Gd = Isat * 3.0 * arg / vd + gmin;
    } else {
        // ---- reverse breakdown ----
        double arg1 = -(tBrkdwnV_ + vd) / vteBRK;
        arg1 = min(CONSTMAX_EXP_ARG, arg1);
        const double evrev = exp(arg1);
        Id = -Isat * evrev + gmin * vd;
        Gd = Isat * evrev / vteBRK + gmin;
    }

    // ---------------- junction charge ----------------
    double Qd = 0.0, Cd = 0.0;
    if (params_.CJO != 0.0) {
        const double Czero = tJctCap_ * Area;
        if (vd < tDepCap_) {
            const double u  = 1.0 - vd / tJctPot_;
            const double sarg = exp(-params_.M * log(u));   // u^(−M)
            Qd = params_.TT * Id + tJctPot_ * Czero * (1.0 - u * sarg) / (1.0 - params_.M);
            Cd = params_.TT * Gd + Czero * sarg;
        } else {
            const double Czof2 = Czero / F2_;
            const double MotJctPot = params_.M / tJctPot_;
            Qd = params_.TT * Id + Czero * tF1_ +
                 Czof2 * (F3_ * (vd - tDepCap_) + 0.5 * MotJctPot *
                          (vd * vd - tDepCap_ * tDepCap_));
            Cd = params_.TT * Gd + Czof2 * (F3_ + MotJctPot * vd);
        }
    }

    vd_limited_ = vd;
    id_limited_ = Id;
    gd_limited_ = Gd;
    qd_limited_ = Qd;
    cd_limited_ = Cd;

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
    limiting_active_ = false;
    if (!first_solve_) {
        if (params_.BVGiven && vd < min(0.0, -params_.BV + 10.0 * vteBRK_)) {
            const double vdtmp = -(params_.BV + vd);
            const double vdtmp_old = -(vd_prev_ + params_.BV);
            const double vnew = pnjlim(vdtmp, vdtmp_old, vteBRK_, tVcrit_);
            if (vnew != vdtmp) { vd = -(vnew + params_.BV); limiting_active_ = true; }
        } else {
            const double vnew = pnjlim(vd, vd_prev_, vte_, tVcrit_);
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
    const double Ir = gspr_ * (v_an - v_pri);
    const double Id = id_limited_;

    if (ip >= 0) F[ip] += Ir;      // current entering anode = resistor current
    if (in >= 0) F[in] -= Id;      // current leaving cathode
    // Internal node KCL (aliased to the anode row when RS = 0 — that is
    // where the junction current enters).
    if (ir >= 0) F[ir] += (Id - Ir);

    first_solve_ = false;
}

void Diode::loadQ(const DaeSolution& x, double time,
                 std::vector<double>& Q) const {
    (void)time;
    const int in = rowOf(nodes_[1]);
    const int ir = priRow();

    const double vd = limitedJunctionVoltage(x);
    evaluateJunction(vd);

    if (in >= 0) Q[in] -= qd_limited_;
    if (ir != in) Q[ir] += qd_limited_;
}

void Diode::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                         double time) const {
    (void)x; (void)time;
    const double Gd = gd_limited_;
    const double Gspr = gspr_;

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

void Diode::confirmAccept() {
    limiting_active_ = false;
}

// ===============================================================
// Outputs
// ===============================================================

double Diode::currentForVoltage(double vd) const {
    // Evaluates ONLY the conduction current (no cache mutation), used by
    // getCurrent()/output paths.
    const double Area = params_.AREA;
    const double Isat = isat_ * Area;
    const double IsatR = isatR_ * Area;
    const double gmin = params_.GMIN;

    if (vd >= -3.0 * vte_) {
        double arg1 = vd / vte_;
        arg1 = std::min(100.0, arg1);
        const double evd = std::exp(arg1);
        const double Inorm = Isat * (evd - 1.0) + gmin * vd;

        double arg2 = vd / vteR_;
        arg2 = std::min(100.0, arg2);
        const double evdR = std::exp(arg2);
        const double Irec = IsatR * (evdR - 1.0);

        double Khi = 1.0;
        if (ikf_ > 0.0 && (ikf_ + Inorm) > 0.0)
            Khi = std::sqrt(ikf_ / (ikf_ + Inorm));

        double Kgen = 0.0;
        if (Irec != 0.0 && tJctPot_ > 0.0) {
            const double u = (1.0 - vd / tJctPot_) * (1.0 - vd / tJctPot_) + 0.005;
            Kgen = std::exp(0.5 * params_.M * std::log(u));
        }
        return Inorm * Khi + Irec * Kgen;
    }
    if (!params_.BVGiven || vd >= -tBrkdwnV_) {
        double arg = 3.0 * vte_ / (vd * CONSTe);
        arg = arg * arg * arg;
        return -Isat * (1.0 + arg) + gmin * vd;
    }
    double arg1 = -(tBrkdwnV_ + vd) / vteBRK_;
    arg1 = std::min(100.0, arg1);
    return -Isat * std::exp(arg1) + gmin * vd;
}

double Diode::getCurrent() const {
    return currentForVoltage(junctionVoltageFromNodes());
}

double Diode::getPower() const {
    return getVoltage() * getCurrent();
}

} // namespace e_XSpice

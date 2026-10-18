/**
 * @author Gabriel Henrique Silva
 * @file MosfetBSIM4Device.cpp
 * @brief BSIM4 level-2 MOSFET device — physics inline (no kernel adapter).
 *
 * Simple, scalable design. Physics computed directly here; the Mosfet
 * wrapper only applies polarity and stamps. If you later want the upstream
 * Berkeley C (b4ld.c/b4temp.c), replace evaluate() — nothing else changes.
 */

#include "MosfetBSIM4Device.h"
#include "../core/DeviceConstants.h"
#include "../../solver/DaeMatrix.h"
#include "../../core/DaeSolution.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace e_XSpice {

// EXSPICE_DEBUG_BSIM4=1 dumps one line per evaluation (diagnostics only:
// off by default, and the fprintf stays out of the Newton hot path).
bool MosfetBSIM4Device::debug_ =
    (std::getenv("EXSPICE_DEBUG_BSIM4") != nullptr);

// ε_SiO2 = 3.9·ε_0
static constexpr double EpsOx = 3.9 * 8.854214871e-12;

MosfetBSIM4Device::MosfetBSIM4Device(
    RCP<const MosfetBSIM4Model> model,
    double W, double L, double M, double AREA)
    : model_(std::move(model)), W_(W), L_(L), M_(M), AREA_(AREA) {
    if (!model_) {
        throw std::invalid_argument("BSIM4 device requires a model card");
    }
    if (W_ <= 0.0 || L_ <= 0.0) {
        throw std::invalid_argument("BSIM4 MOSFET W and L must be positive");
    }
}

double MosfetBSIM4Device::pnjlim(double vnew, double vold, double vt,
                                 double vcrit) {
    if (vnew > vcrit && std::abs(vnew - vold) > 2.0 * vt) {
        if (vold > 0.0) {
            const double arg = 1.0 + (vnew - vold) / vt;
            if (arg > 0.0) {
                vnew = vold + vt * std::log(arg);
            } else {
                vnew = vcrit;
            }
        } else {
            vnew = vt * std::log(vnew / vt);
        }
    }
    return vnew;
}

double MosfetBSIM4Device::clampStep(double vnew, double vold,
                                    double maxStep) {
    const double d = vnew - vold;
    if (d > maxStep) return vold + maxStep;
    if (d < -maxStep) return vold - maxStep;
    return vnew;
}

void MosfetBSIM4Device::limitedVoltages(const DaeSolution& x,
                                         int rd, int rg, int rs, int rb,
                                         bool allow_limiting) {
    const double vd = (rd >= 0) ? x.data()[rd] : 0.0;
    const double vg = (rg >= 0) ? x.data()[rg] : 0.0;
    const double vs = (rs >= 0) ? x.data()[rs] : 0.0;
    const double vb = (rb >= 0) ? x.data()[rb] : 0.0;

    // N-oriented controlling voltages (p applied by wrapper).
    double vgs = p_ * (vg - vs);
    double vds = p_ * (vd - vs);
    double vbs = p_ * (vb - vs);
    const double vbd = vbs - vds;

    limiting_active_ = false;
    if (allow_limiting && !first_solve_) {
        const double nvgs = clampStep(vgs, vgs_prev_, MaxStep_);
        if (nvgs != vgs) { vgs = nvgs; limiting_active_ = true; }
        const double nvds = clampStep(vds, vds_prev_, MaxStep_);
        if (nvds != vds) { vds = nvds; limiting_active_ = true; }
        const double nvbs = pnjlim(vbs, vbs_prev_, thermal_->vt,
                                   thermal_->vcrit);
        if (nvbs != vbs) { vbs = nvbs; limiting_active_ = true; }
        const double nvbd = pnjlim(vbd, vbd_prev_, thermal_->vt,
                                   thermal_->vcrit);
        if (nvbd != vbd) {
            vbs = nvbd + vds;
            limiting_active_ = true;
        }
    }
    vgs_lim_ = vgs;
    vds_lim_ = vds;
    vbs_lim_ = vbs;
}

void MosfetBSIM4Device::beginSolve(bool initial_solve) {
    if (initial_solve) {
        // Start from the most linear point of each junction (same Xyce
        // initJctFlag_ behaviour as Mosfet::beginSolve/Bjt::beginSolve).
        vbs_prev_ = thermal_ ? thermal_->vcrit : 0.5;
        vbd_prev_ = 0.0;
        vgs_prev_ = 0.0;
        vds_prev_ = 0.0;
        first_solve_ = true;
    }
}

void MosfetBSIM4Device::resetHistory(double vgs_raw, double vds_raw,
                                     double vbs_raw) {
    limiting_active_ = false;
    vbs_prev_ = vbs_raw;
    vbd_prev_ = vbs_raw - vds_raw;
    vgs_prev_ = vgs_raw;
    vds_prev_ = vds_raw;
}

void MosfetBSIM4Device::applySolutionLimit(std::vector<double>& x,
                                           int rd, int rg, int rs, int rb,
                                           int p) const {
    if (!limiting_active_) return;

    // Gate-anchored write-back (mirrors Mosfet::applySolutionLimit): the
    // gate is the control terminal, drain/source/bulk carry the excursion.
    const double vg = (rg >= 0) ? x[rg] : 0.0;
    if (rs >= 0) x[rs] = vg - p * vgs_lim_;
    const double vs = (rs >= 0) ? x[rs] : 0.0;
    if (rd >= 0) x[rd] = vs + p * vds_lim_;
    if (rb >= 0) x[rb] = vs + p * vbs_lim_;
}

void MosfetBSIM4Device::commitLimited() {
    vbs_prev_ = vbs_lim_;
    vbd_prev_ = vbs_lim_ - vds_lim_;
    vgs_prev_ = vgs_lim_;
    vds_prev_ = vds_lim_;
    first_solve_ = false;
}

void MosfetBSIM4Device::evaluate(double vgs, double vds, double vbs) {
    const MosfetBSIM4Thermal& t = *thermal_;
    const double mult = AREA_;
    const double gmin = t.gmin;

    // ---- source/drain swap for vds < 0 (device is symmetric) ----
    const bool sw = (vds < 0.0);
    swapped_ = sw;
    const double vgs_e = sw ? vgs - vds : vgs;
    const double vds_e = sw ? -vds : vds;
    const double vbs_e = sw ? vbs - vds : vbs;

    vgs_cached_ = vgs;
    vds_cached_ = vds;
    vbs_cached_ = vbs;

    // ---- threshold voltage with body effect ----
    // VTH0 carries the SPICE sign (negative for PMOS); the evaluation frame
    // is N-oriented, so the polarity is applied here exactly like the
    // level-1 wrapper does with VTO (without it a PMOS would see a NEGATIVE
    // threshold and get stronger as the gate approaches the source).
    const double vto_n = p_ * model_->VTH0;
    const double sqrtPhi = std::sqrt(std::max(model_->PHI, 1e-6));
    const double phib = model_->PHI - vbs_e;
    double dvth_dvbs = 0.0;
    double vth = vto_n;
    if (phib > 1e-3) {
        const double sq = std::sqrt(phib);
        vth += model_->GAMMA * (sq - sqrtPhi);
        dvth_dvbs = -model_->GAMMA / (2.0 * sq);
    } else {
        vth += model_->GAMMA * (std::sqrt(1e-3) - sqrtPhi);
    }
    vth_cached_ = vth;

    // ---- effective mobility ----
    const double ov = vgs_e - vth;
    double ueff = t.u0;
    const double leff = L_ - 2.0 * model_->LD;
    if (leff <= 0.0) {
        throw std::runtime_error("BSIM4: Leff = L - 2*LD must be positive");
    }
    const double denom_mob = 1.0 + (t.ua + 2.0 * t.ub * ov) * ov / leff;
    ueff = t.u0 / denom_mob;

    // ---- beta = KP * (W/Leff) * M ----
    const double beta = model_->KP * (W_ / leff) * M_;
    const double beta_over_u0 = beta * ueff / t.u0;

    // ---- drain current ----
    double Id = 0.0, gm = 0.0, gds = 0.0, gmbs = 0.0;
    saturated_ = false;

    if (ov > 0.0) {
        const double vdsat = ov;
        const double lame = 1.0 + model_->LAMBDA * vds_e;
        if (vds_e < vdsat) {
            // Linear (triode) region
            Id = beta_over_u0 * (ov * vds_e - 0.5 * vds_e * vds_e) * lame;
            gm = beta_over_u0 * vds_e * lame;
            gds = beta_over_u0 * (ov - vds_e) * lame +
                  beta_over_u0 * (ov * vds_e - 0.5 * vds_e * vds_e) *
                  model_->LAMBDA;
            gmbs = -gm * dvth_dvbs;
            saturated_ = false;
        } else {
            // Saturation region
            Id = 0.5 * beta_over_u0 * ov * ov * lame;
            gm = beta_over_u0 * ov * lame;
            gds = 0.5 * beta_over_u0 * ov * ov * model_->LAMBDA;
            gmbs = -gm * dvth_dvbs;
            saturated_ = true;
        }
    } else {
        Id = 0.0;
        gm = 0.0;
        gds = 0.0;
        gmbs = 0.0;
        saturated_ = false;
    }

    // ---- bulk diodes (PHYSICAL frame, symmetric) ----
    const double vbd_p = vbs - vds;
    const double isBulk = t.isat * mult;
    double arg = std::min(Const::MaxExpArg, vbs / t.vt);
    const double evbs = std::exp(arg);
    const double Ibs = isBulk * (evbs - 1.0) + gmin * vbs;
    arg = std::min(Const::MaxExpArg, vbd_p / t.vt);
    const double evbd = std::exp(arg);
    const double Ibd = isBulk * (evbd - 1.0) + gmin * vbd_p;

    // ---- map channel to physical frame ----
    const double sid = sw ? -1.0 : 1.0;
    id_ = sid * Id + Ibd;
    is_ = -sid * Id + Ibs;
    ib_ = -(Ibs + Ibd);

    gm_ = gm;
    gds_ = gds;
    gmbs_ = gmbs;

    if (debug_) {
        std::fprintf(stderr,
            "[B4 p=%+d] vgs=%+.6g vds=%+.6g vbs=%+.6g | vth=%+.6g ov=%+.6g"
            " u0=%.6g ua=%.6g ub=%.6g ueff=%.6g W=%.6g Leff=%.6g KP=%.6g"
            " beta=%.6g Id=%+.6g gm=%.6g gds=%.6g sat=%d\n",
            p_, vgs, vds, vbs, vth, ov, t.u0, t.ua, t.ub, ueff, W_, leff,
            model_->KP, beta_over_u0, id_, gm_, gds_, (int)saturated_);
    }

    // ---- charges (Meyer partition) ----
    const double cox = EpsOx / model_->TOX * W_ * leff * mult;
    const double cgso = model_->CGSO * W_ * mult;
    const double cgdo = model_->CGDO * W_ * mult;
    const double cgbo = model_->CGBO * leff * mult;

    double Cgs_p, Cgd_p, Cgb_p;
    if (!sw) {
        if (saturated_) {
            Cgs_p = (2.0 / 3.0) * cox + cgso;
            Cgd_p = cgdo;
        } else if (ov > 0.0) {
            Cgs_p = 0.5 * cox + cgso;
            Cgd_p = 0.5 * cox + cgdo;
        } else {
            Cgs_p = cgso;
            Cgd_p = cgdo;
        }
        Cgb_p = cox + cgbo;
    } else {
        if (saturated_) {
            Cgs_p = cgdo;
            Cgd_p = (2.0 / 3.0) * cox + cgso;
        } else if (ov > 0.0) {
            Cgs_p = cgdo;
            Cgd_p = 0.5 * cox + cgso;
        } else {
            Cgs_p = cgdo;
            Cgd_p = cgso;
        }
        Cgb_p = cox + cgbo;
    }

    // Bulk depletion (exact Q/C from shared junction model)
    const double factor = t.bulkCap.cjo;
    const JunctionGeom geomS{model_->CBS * mult * factor, 0.0, model_->MJ};
    const JunctionGeom geomD{model_->CBD * mult * factor, 0.0, model_->MJ};
    const JunctionChargePoint jbs =
        junctionCharge(vbs, 0.0, 0.0, geomS, t.bulkCap);
    const JunctionChargePoint jbd =
        junctionCharge(vbd_p, 0.0, 0.0, geomD, t.bulkCap);
    const double Cbs = jbs.c, Cbd = jbd.c;

    // Nodal charges (charge-conserving by construction)
    const double vgd_p = vgs - vds;
    const double vgb_p = vgs - vbs;
    qg_ = Cgs_p * vgs + Cgd_p * vgd_p + Cgb_p * vgb_p;
    qd_ = -Cgd_p * vgd_p - jbd.q;
    qb_ = -Cgb_p * vgb_p + jbs.q + jbd.q;

    // Capacitances (dQ/dV)
    cgs_ = Cgs_p;
    cgd_ = Cgd_p;
    cgb_ = Cgb_p;
    cdg_ = -Cgd_p;
    cdd_ = Cgd_p;
    cdb_ = -Cbd;
    csg_ = -Cgs_p;
    csd_ = 0.0;
    css_ = Cgs_p;
    csb_ = -Cbs;
    cbg_ = -Cgb_p;
    cbd_ = -Cbd;
    cbs_ = -Cbs;
    cbb_ = Cgb_p;
}

void MosfetBSIM4Device::loadF(std::vector<double>& F,
                              int rd, int rg, int rs, int rb) const {
    // Terminal order: D, G, S, B. Terminal currents are N-oriented
    // (positive INTO the device); p_ maps them to the physical frame, so
    // the stamp below is the current leaving each node of a PMOS.
    const int r[4] = {rd, rg, rs, rb};
    const double iv[4] = {id_, 0.0, is_, ib_};
    for (int i = 0; i < 4; ++i) {
        if (r[i] >= 0 && static_cast<size_t>(r[i]) < F.size()) {
            F[r[i]] += p_ * iv[i];
        }
    }
}

void MosfetBSIM4Device::loadQ(std::vector<double>& Q,
                              int rd, int rg, int rs, int rb) const {
    // Nodal charges (N-oriented): Qd, Qg, Qs derived by conservation, Qb.
    const int r[4] = {rd, rg, rs, rb};
    const double qv[4] = {qd_, qg_, -(qg_ + qd_ + qb_), qb_};
    for (int i = 0; i < 4; ++i) {
        if (r[i] >= 0 && static_cast<size_t>(r[i]) < Q.size()) {
            Q[r[i]] += p_ * qv[i];
        }
    }
}

void MosfetBSIM4Device::loadJacobianF(DaeMatrix& J,
                                      const std::vector<int>& slots) const {
    // slots: 16 entries in row-major order (D,G,S,B x D,G,S,B)
    double ef[16] = {0};

    const double gm = gm_;
    const double gds = gds_;
    const double gmbs = gmbs_;

    if (!swapped_) {
        ef[0*4+1] = gm;                              // d/dg
        ef[0*4+0] = gds;                             // d/dd
        ef[0*4+2] = -(gm + gds + gmbs);              // d/ds
        ef[0*4+3] = gmbs;                            // d/db
        ef[2*4+1] = -gm;                             // s/dg
        ef[2*4+0] = -gds;                            // s/dd
        ef[2*4+2] = (gm + gds + gmbs);               // s/ds
        ef[2*4+3] = -gmbs;                           // s/db
    } else {
        ef[0*4+1] = -gm;                             // d/dg
        ef[0*4+0] = (gm + gds + gmbs);               // d/dd
        ef[0*4+2] = -gds;                            // d/ds
        ef[0*4+3] = -gmbs;                           // d/db
        ef[2*4+1] = gm;                             // s/dg
        ef[2*4+0] = -(gm + gds + gmbs);              // s/dd
        ef[2*4+2] = gds;                             // s/ds
        ef[2*4+3] = gmbs;                            // s/db
    }
    // Bulk terminal: no separate bulk diode derivative in this subset.

    for (int k = 0; k < 16; ++k) {
        if (slots[k] >= 0) {
            J.accumulateSlot(slots[k], ef[k]);
        }
    }
}

void MosfetBSIM4Device::loadJacobianQ(DaeMatrix& J,
                                      const std::vector<int>& slots) const {
    double eq[16] = {0};

    eq[1*4+1] = cgs_ + cgd_ + cgb_;  // g/dg
    eq[1*4+0] = -cgd_;               // g/dd
    eq[1*4+2] = -cgs_;               // g/ds
    eq[1*4+3] = -cgb_;               // g/db
    eq[0*4+1] = -cgd_;               // d/dg
    eq[0*4+0] = cgd_ + cbd_;         // d/dd
    eq[0*4+3] = -cbd_;               // d/db
    eq[2*4+1] = -cgs_;               // s/dg
    eq[2*4+2] = cgs_ + cbs_;         // s/ds
    eq[2*4+3] = -cbs_;               // s/db
    eq[3*4+1] = -cgb_;               // b/dg
    eq[3*4+0] = -cbd_;               // b/dd
    eq[3*4+2] = -cbs_;               // b/ds
    eq[3*4+3] = cgb_ + cbs_ + cbd_;  // b/db

    for (int k = 0; k < 16; ++k) {
        if (slots[k] >= 0) {
            J.accumulateSlot(slots[k], eq[k]);
        }
    }
}

} // namespace e_XSpice

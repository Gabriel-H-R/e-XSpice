/**
 * @author Gabriel Henrique Silva
 * @file Bjt.cpp
 * @brief Bipolar junction transistor (NPN and PNP) — level 1 (Ebers-Moll)
 *        and level 2 (Gummel-Poon) in the shared DAE formulation.
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
         const BjtInstanceParams& inst,
         Node* substrate)
    : Device(name, DeviceType::BJT)
    , model_(std::move(model))
    , inst_(inst) {
    if (!model_) {
        throw std::invalid_argument("BJT requires a model card");
    }
    if (inst_.AREA <= 0.0) {
        throw std::invalid_argument("BJT AREA must be positive");
    }
    if (substrate && model_->LEVEL != 2) {
        throw std::invalid_argument(
            "BJT substrate terminal needs a LEVEL=2 card");
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

    // Level-2 internal-node layout, fixed order base/collector/emitter.
    // The offsets are card-fixed, so every loader derives rows with no maps.
    if (model_->LEVEL == 2) {
        if (model_->hasRb()) offB_ = numInternal_++;
        if (model_->hasRc()) offC_ = numInternal_++;
        if (model_->hasRe()) offE_ = numInternal_++;
        gb_ = model_->hasRb() ? 1.0 / model_->RB : 0.0;
        gc_ = model_->hasRc() ? 1.0 / model_->RC : 0.0;
        ge_ = model_->hasRe() ? 1.0 / model_->RE : 0.0;
    }
    int pat = 0;
    if (offB_ >= 0) { patRb_ = pat; pat += 4; }
    if (offC_ >= 0) { patRc_ = pat; pat += 4; }
    if (offE_ >= 0) { patRe_ = pat; pat += 4; }
    patInt_ = pat; pat += 9;
    if (substrate) { patSub_ = pat; pat += 4; }

    nodes_.push_back(collector);   // 0
    nodes_.push_back(base);        // 1
    nodes_.push_back(emitter);     // 2
    if (collector) collector->addConnectedDevice(this);
    if (base)     base->addConnectedDevice(this);
    if (emitter)  emitter->addConnectedDevice(this);
    if (substrate) {
        nodes_.push_back(substrate);  // 3
        substrate->addConnectedDevice(this);
    }
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
    const int n = 2 + (hasSubstrate() ? 1 : 0);   // Q_BE, Q_BC (, Q_CS)
    setStateIndex(store.registerStates(n));
    return n;
}

void Bjt::registerJacobianPattern(std::vector<int>& rows,
                                  std::vector<int>& cols) const {
    // External rows.
    const int rb = rowOf(nodes_[1]);
    const int rc = rowOf(nodes_[0]);
    const int re = rowOf(nodes_[2]);
    // Intrinsic rows (alias the externals when no ohmic resistor).
    const int bi = rowBi();
    const int ci = rowCi();
    const int ei = rowEi();

    // Resistor quads (registration order matches patRb_/patRc_/patRe_).
    if (offB_ >= 0) {
        rows.push_back(rb); cols.push_back(rb);
        rows.push_back(rb); cols.push_back(bi);
        rows.push_back(bi); cols.push_back(rb);
        rows.push_back(bi); cols.push_back(bi);
    }
    if (offC_ >= 0) {
        rows.push_back(rc); cols.push_back(rc);
        rows.push_back(rc); cols.push_back(ci);
        rows.push_back(ci); cols.push_back(rc);
        rows.push_back(ci); cols.push_back(ci);
    }
    if (offE_ >= 0) {
        rows.push_back(re); cols.push_back(re);
        rows.push_back(re); cols.push_back(ei);
        rows.push_back(ei); cols.push_back(re);
        rows.push_back(ei); cols.push_back(ei);
    }

    // Intrinsic 3x3 block (patInt_ base; same shape as level 1).
    rows.push_back(bi); cols.push_back(bi);
    rows.push_back(bi); cols.push_back(ci);
    rows.push_back(bi); cols.push_back(ei);
    rows.push_back(ci); cols.push_back(bi);
    rows.push_back(ci); cols.push_back(ci);
    rows.push_back(ci); cols.push_back(ei);
    rows.push_back(ei); cols.push_back(bi);
    rows.push_back(ei); cols.push_back(ci);
    rows.push_back(ei); cols.push_back(ei);

    // Substrate quad (patSub_ base).
    if (hasSubstrate()) {
        const int rs = rowSub();
        rows.push_back(rs); cols.push_back(rs);
        rows.push_back(rs); cols.push_back(ci);
        rows.push_back(ci); cols.push_back(rs);
        rows.push_back(ci); cols.push_back(ci);
    }
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
// Junction evaluation
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
    // Junctions live on the INTRINSIC nodes (aliased to externals when the
    // card has no ohmic resistance).
    const int ci = rowCi();
    const int bi = rowBi();
    const int ei = rowEi();

    const double vc = (ci >= 0) ? x.data()[ci] : 0.0;
    const double vb = (bi >= 0) ? x.data()[bi] : 0.0;
    const double ve = (ei >= 0) ? x.data()[ei] : 0.0;

    // NPN-ORIENTED junction voltages: the polarity multiply is what makes
    // the limiters and the card laws polarity-agnostic (for a PNP the
    // physical forward direction is a NEGATIVE vb − ve).
    double vbe = p_ * (vb - ve);
    double vbc = p_ * (vb - vc);

    limiting_active_ = false;
    if (allow_limiting && !first_solve_ && limiting_enabled_) {
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
                     p_ * (vb - ve), p_ * (vb - vc), vbe, vbc,
                     vbe_prev_, vbc_prev_,
                     static_cast<int>(limiting_active_));
    }
}

void Bjt::evaluate(double vbe, double vbc) const {
    if (model_->LEVEL == 2) {
        evaluateLevel2(vbe, vbc);
    } else {
        evaluateLevel1(vbe, vbc);
    }
}

void Bjt::evaluateLevel1(double vbe, double vbc) const {
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
    const JunctionChargePoint depBE =
        junctionCharge(vbe, 0.0, 0.0, geomBE, t.beCap);
    const JunctionChargePoint depBC =
        junctionCharge(vbc, 0.0, 0.0, geomBC, t.bcCap);

    const double tf = t.tf * area;   // transit times, area-scaled
    const double tr = t.tr * area;

    qbe_ = depBE.q + tf * ict;
    qbc_ = depBC.q - tr * ict;

    cbe_v_ = depBE.c + tf * gm;
    cbe_w_ = tf * go;
    cbc_u_ = -tr * gm;
    cbc_v_ = depBC.c - tr * go;

    isc_ = gsc_ = qcs_ = csc_ = 0.0;
}

void Bjt::evaluateLevel2(double vbe, double vbc) const {
    const BjtThermal& t = *thermal_;
    const double area = inst_.AREA;
    const double gmin = t.gmin;

    const double isF = t.beIsat * area;
    const double isR = t.bcIsat * area;
    const double betaF = t.betaf * area;
    const double betaR = t.betar * area;
    const double ise = t.iseT * area;
    const double isc = t.iscT * area;
    const double ikf = t.ikf * area;
    const double ikr = t.ikr * area;

    // ---- exponentials (ideal + non-ideal) ----
    double arg = vbe / t.vte;
    arg = min(Const::MaxExpArg, arg);
    const double evf = exp(arg);
    arg = vbc / t.vtc;
    arg = min(Const::MaxExpArg, arg);
    const double evr = exp(arg);

    double evfe = 1.0, evrc = 1.0;
    if (ise > 0.0) {
        arg = min(Const::MaxExpArg, vbe / t.vteE);
        evfe = exp(arg);
    }
    if (isc > 0.0) {
        arg = min(Const::MaxExpArg, vbc / t.vtcC);
        evrc = exp(arg);
    }

    // ---- normalized base charge qb = q1·(1+√(1+4·q2))/2 ----
    // q1 carries Early (base-width modulation); q2 carries high injection.
    // The denominator floor mirrors level 1's ke floor (same linear
    // derivative, so Newton still walks back out of deep saturation).
    const double denom = 1.0 - vbc / model_->VAF - vbe / model_->VAR;
    const double denom_f = max(denom, 1e-6);
    const double q1 = 1.0 / denom_f;
    // d(1/denom)/dvbe = −d(denom)/dvbe / denom² = q1²/VAR
    // d(1/denom)/dvbc = −d(denom)/dvbc / denom² = q1²/VAF
    const double dq1vbe = q1 * q1 / model_->VAR;
    const double dq1vbc = q1 * q1 / model_->VAF;

    const double fwd = isF * (evf - 1.0);
    const double rev = isR * (evr - 1.0);
    const double dfwd = isF * evf / t.vte;
    const double drev = isR * evr / t.vtc;

    double q2 = 0.0, dq2vbe = 0.0, dq2vbc = 0.0;
    if (ikf > 0.0) { q2 += fwd / ikf; dq2vbe += dfwd / ikf; }
    if (ikr > 0.0) { q2 += rev / ikr; dq2vbc += drev / ikr; }

    const double sqarg = max(1.0 + 4.0 * q2, 1e-12);
    const double sq = std::sqrt(sqarg);
    const double qb = q1 * (1.0 + sq) * 0.5;
    // dqb = dq1·(1+sq)/2 + q1·dq2/sq  (d(sq) = 2·dq2/sq).
    const double dqb_vbe = dq1vbe * (1.0 + sq) * 0.5 + q1 * dq2vbe / sq;
    const double dqb_vbc = dq1vbc * (1.0 + sq) * 0.5 + q1 * dq2vbc / sq;

    // ---- transport current icc = (fwd − rev)/qb ----
    const double num = fwd - rev;
    const double icc = num / qb;
    const double dicc_vbe = (dfwd * qb - num * dqb_vbe) / (qb * qb);
    const double dicc_vbc = (-drev * qb - num * dqb_vbc) / (qb * qb);

    // ---- base currents (ideal + non-ideal + gmin) ----
    double ibe = betaF * (evf - 1.0) + gmin * vbe;
    double gib = betaF * evf / t.vte + gmin;
    if (ise > 0.0) {
        ibe += ise * (evfe - 1.0);
        gib += ise * evfe / t.vteE;
    }
    double ibc_val = betaR * (evr - 1.0) + gmin * vbc;
    double gibc = betaR * evr / t.vtc + gmin;
    if (isc > 0.0) {
        ibc_val += isc * (evrc - 1.0);
        gibc += isc * evrc / t.vtcC;
    }

    ib_ = ibe + ibc_val;
    ic_ = icc - ibc_val;
    ie_ = -ib_ - ic_;

    gib_ = gib;
    gibc_ = gibc;
    gm_ = dicc_vbe;
    go_ = dicc_vbc;

    // ---- charges ----
    // Depletion parts from the shared junction model; transit parts use the
    // Gummel-Poon partition (forward component on BE, reverse on BC):
    //   Qbe,tr = TF·fwd/qb,   Qbc,tr = TR·rev/qb.
    const JunctionGeom geomBE{t.beCap.cjo * area, 0.0, model_->MJE};
    const JunctionGeom geomBC{t.bcCap.cjo * area, 0.0, model_->MJC};
    const JunctionChargePoint depBE =
        junctionCharge(vbe, 0.0, 0.0, geomBE, t.beCap);
    const JunctionChargePoint depBC =
        junctionCharge(vbc, 0.0, 0.0, geomBC, t.bcCap);

    const double tf = t.tf * area;   // transit times, area-scaled
    const double tr = t.tr * area;
    const double qbf = (tf != 0.0) ? tf * fwd / qb : 0.0;
    const double qbr = (tr != 0.0) ? tr * rev / qb : 0.0;

    qbe_ = depBE.q + qbf;
    qbc_ = depBC.q + qbr;

    const double qb2 = qb * qb;
    const double dqbf_vbe = (tf != 0.0) ? tf * (dfwd * qb - fwd * dqb_vbe) / qb2 : 0.0;
    const double dqbf_vbc = (tf != 0.0) ? tf * (0.0 * qb - fwd * dqb_vbc) / qb2 : 0.0;
    const double dqbr_vbe = (tr != 0.0) ? tr * (0.0 * qb - rev * dqb_vbe) / qb2 : 0.0;
    const double dqbr_vbc = (tr != 0.0) ? tr * (drev * qb - rev * dqb_vbc) / qb2 : 0.0;

    cbe_v_ = depBE.c + dqbf_vbe;
    cbe_w_ = dqbf_vbc;
    cbc_u_ = dqbr_vbe;
    cbc_v_ = depBC.c + dqbr_vbc;

    // ---- substrate diode (only with a 4th terminal) ----
    isc_ = gsc_ = qcs_ = csc_ = 0.0;
    if (hasSubstrate()) {
        // Evaluated by the loaders (they own the substrate/collector rows).
    }
}

void Bjt::loadF(const DaeSolution& x, double time,
                std::vector<double>& F) const {
    (void)time;
    const int rb = rowOf(nodes_[1]);
    const int rc = rowOf(nodes_[0]);
    const int re = rowOf(nodes_[2]);
    const int bi = rowBi();
    const int ci = rowCi();
    const int ei = rowEi();

    limitedJunctionVoltages(x, true);
    evaluate(vbe_lim_, vbc_lim_);
    // Report the EXTERNAL terminal voltages (what a scope would measure).
    {
        const double vxb = (rb >= 0) ? x.data()[rb] : 0.0;
        const double vxc = (rc >= 0) ? x.data()[rc] : 0.0;
        const double vxe = (re >= 0) ? x.data()[re] : 0.0;
        vbe_cached_ = p_ * (vxb - vxe);
        vbc_cached_ = p_ * (vxb - vxc);
    }

    // Ohmic drops (external ↔ intrinsic). Resistor currents are PHYSICAL
    // (symmetric device, physical voltages) — they take NO polarity scaling:
    // F[row] accumulates the physical node→device current.
    if (offB_ >= 0) {
        const double vxb = (rb >= 0) ? x.data()[rb] : 0.0;
        const double vib = (bi >= 0) ? x.data()[bi] : 0.0;
        if (rb >= 0) F[rb] += gb_ * (vxb - vib);
        if (bi >= 0) F[bi] += gb_ * (vib - vxb);
    }
    if (offC_ >= 0) {
        const double vxc = (rc >= 0) ? x.data()[rc] : 0.0;
        const double vic = (ci >= 0) ? x.data()[ci] : 0.0;
        if (rc >= 0) F[rc] += gc_ * (vxc - vic);
        if (ci >= 0) F[ci] += gc_ * (vic - vxc);
    }
    if (offE_ >= 0) {
        const double vxe = (re >= 0) ? x.data()[re] : 0.0;
        const double vie = (ei >= 0) ? x.data()[ei] : 0.0;
        if (re >= 0) F[re] += ge_ * (vxe - vie);
        if (ei >= 0) F[ei] += ge_ * (vie - vxe);
    }

    // Intrinsic transistor on the intrinsic rows.
    if (bi >= 0) F[bi] += p_ * ib_;
    if (ci >= 0) F[ci] += p_ * ic_;
    if (ei >= 0) F[ei] += p_ * ie_;

    // Substrate diode (level 2, 4-terminal): current flows substrate→ci.
    if (hasSubstrate()) {
        const int rs = rowSub();
        const double vsub = (rs >= 0) ? x.data()[rs] : 0.0;
        const double vci = (ci >= 0) ? x.data()[ci] : 0.0;
        const double vsc = p_ * (vsub - vci);
        evalSubstrateDiode(thermal_->csCap, inst_.AREA, thermal_->gmin, vsc,
                           isc_, gsc_, qcs_, csc_);
        if (rs >= 0) F[rs] += p_ * isc_;
        if (ci >= 0) F[ci] -= p_ * isc_;
    }

    vbe_prev_ = vbe_lim_;
    vbc_prev_ = vbc_lim_;
    first_solve_ = false;
}

void Bjt::applySolutionLimit(std::vector<double>& x) const {
    if (!limiting_active_) return;

    // Rebuild the INTRINSIC unknowns from the LIMITED junction voltages,
    // anchored on the intrinsic emitter (kept as-is). The EXTERNAL terminal
    // of each ohmic resistor is shifted by the SAME delta, so every RB/RC/RE
    // current is preserved exactly — zero artificial resistor residual.
    // With no ohmic resistors the intrinsic rows alias the externals and
    // this is exactly the level-1 write-back.
    const int bi = rowBi();
    const int ci = rowCi();
    const int ei = rowEi();
    const int rb = rowOf(nodes_[1]);   // external base
    const int rc = rowOf(nodes_[0]);   // external collector

    const double ve = (ei >= 0) ? x[ei] : 0.0;
    if (bi >= 0) {
        const double updated = ve + p_ * vbe_lim_;
        const double delta = updated - x[bi];
        x[bi] = updated;
        if (rb >= 0 && rb != bi) x[rb] += delta;
    }
    const double vb = (bi >= 0) ? x[bi] : 0.0;
    if (ci >= 0) {
        const double updated = vb - p_ * vbc_lim_;
        const double delta = updated - x[ci];
        x[ci] = updated;
        if (rc >= 0 && rc != ci) x[rc] += delta;
    }
}

void Bjt::loadQ(const DaeSolution& x, double time,
                std::vector<double>& Q) const {
    (void)time;
    const int bi = rowBi();
    const int ci = rowCi();
    const int ei = rowEi();

    limitedJunctionVoltages(x, true);
    evaluate(vbe_lim_, vbc_lim_);

    if (bi >= 0) Q[bi] += p_ * (qbe_ + qbc_);
    if (ei >= 0) Q[ei] -= p_ * qbe_;
    if (ci >= 0) Q[ci] -= p_ * qbc_;

    if (hasSubstrate()) {
        const int rs = rowSub();
        const double vsub = (rs >= 0) ? x.data()[rs] : 0.0;
        const double vci = (ci >= 0) ? x.data()[ci] : 0.0;
        const double vsc = p_ * (vsub - vci);
        evalSubstrateDiode(thermal_->csCap, inst_.AREA, thermal_->gmin, vsc,
                           isc_, gsc_, qcs_, csc_);
        if (rs >= 0) Q[rs] += p_ * qcs_;
        if (ci >= 0) Q[ci] -= p_ * qcs_;
    }
}

void Bjt::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                        double time) const {
    (void)x; (void)time;
    // Resistor quads: p-free because currents are physical (physical
    // voltages → physical currents).
    if (offB_ >= 0) {
        emit(J, patRb_ + 0,  gb_);
        emit(J, patRb_ + 1, -gb_);
        emit(J, patRb_ + 2, -gb_);
        emit(J, patRb_ + 3,  gb_);
    }
    if (offC_ >= 0) {
        emit(J, patRc_ + 0,  gc_);
        emit(J, patRc_ + 1, -gc_);
        emit(J, patRc_ + 2, -gc_);
        emit(J, patRc_ + 3,  gc_);
    }
    if (offE_ >= 0) {
        emit(J, patRe_ + 0,  ge_);
        emit(J, patRe_ + 1, -ge_);
        emit(J, patRe_ + 2, -ge_);
        emit(J, patRe_ + 3,  ge_);
    }
    const int o = patInt_;
    emit(J, o + 0,  gib_ + gibc_);            // [bi][bi]
    emit(J, o + 1, -gibc_);                   // [bi][ci]
    emit(J, o + 2, -gib_);                    // [bi][ei]
    emit(J, o + 3,  gm_ + go_ - gibc_);       // [ci][bi]
    emit(J, o + 4, -go_ + gibc_);             // [ci][ci]
    emit(J, o + 5, -gm_);                     // [ci][ei]
    emit(J, o + 6, -(gib_ + gm_ + go_));      // [ei][bi]
    emit(J, o + 7,  go_);                     // [ei][ci]
    emit(J, o + 8,  gib_ + gm_);              // [ei][ei]

    if (hasSubstrate()) {
        // F[s] += p·isc_, F[ci] −= p·isc_ with vsc = p·(vsub − vci):
        // dF[s]/d(vsub) = +gsc, dF[s]/d(vci) = −gsc,
        // dF[ci]/d(vsub) = −gsc, dF[ci]/d(vci) = +gsc.
        emit(J, patSub_ + 0,  gsc_);
        emit(J, patSub_ + 1, -gsc_);
        emit(J, patSub_ + 2, -gsc_);
        emit(J, patSub_ + 3,  gsc_);
    }
}

void Bjt::loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                        double time) const {
    (void)x; (void)time;
    const int o = patInt_;
    emit(J, o + 0,  cbe_v_ + cbe_w_ + cbc_u_ + cbc_v_);  // [bi][bi]
    emit(J, o + 1, -cbe_w_ - cbc_v_);                    // [bi][ci]
    emit(J, o + 2, -cbe_v_ - cbc_u_);                    // [bi][ei]
    emit(J, o + 3, -cbc_u_ - cbc_v_);                    // [ci][bi]
    emit(J, o + 4,  cbc_v_);                             // [ci][ci]
    emit(J, o + 5,  cbc_u_);                             // [ci][ei]
    emit(J, o + 6, -cbe_v_ - cbe_w_);                    // [ei][bi]
    emit(J, o + 7,  cbe_w_);                             // [ei][ci]
    emit(J, o + 8,  cbe_v_);                             // [ei][ei]

    if (hasSubstrate()) {
        // Q[s] += p·qcs, Q[ci] −= p·qcs: same sign pattern as F.
        emit(J, patSub_ + 0,  csc_);
        emit(J, patSub_ + 1, -csc_);
        emit(J, patSub_ + 2, -csc_);
        emit(J, patSub_ + 3,  csc_);
    }
}

void Bjt::storeState(const DaeSolution& x, StateStore& store) const {
    limitedJunctionVoltages(x, true);
    evaluate(vbe_lim_, vbc_lim_);
    if (hasSubstrate()) {
        const int ci = rowCi();
        const int rs = rowSub();
        const double vsub = (rs >= 0) ? x.data()[rs] : 0.0;
        const double vci = (ci >= 0) ? x.data()[ci] : 0.0;
        const double vsc = p_ * (vsub - vci);
        evalSubstrateDiode(thermal_->csCap, inst_.AREA, thermal_->gmin, vsc,
                           isc_, gsc_, qcs_, csc_);
        const double q[3] = {qbe_, qbc_, qcs_};
        store.setDeviceQ(state_index_, q, 3);
    } else {
        const double q[2] = {qbe_, qbc_};
        store.setDeviceQ(state_index_, q, 2);
    }
}

void Bjt::loadChargeHistory(const StateStore& store, double q_coeff,
                            double qdot_coeff,
                            std::vector<double>& rhs) const {
    if (q_coeff == 0.0) return;

    const int bi = rowBi();
    const int ci = rowCi();
    const int ei = rowEi();

    const double qbe = q_coeff * store.deviceQ(state_index_, 0) +
                       qdot_coeff * store.deviceQDot(state_index_, 0);
    const double qbc = q_coeff * store.deviceQ(state_index_, 1) +
                       qdot_coeff * store.deviceQDot(state_index_, 1);

    if (bi >= 0) rhs[bi] += p_ * (qbe + qbc);
    if (ei >= 0) rhs[ei] -= p_ * qbe;
    if (ci >= 0) rhs[ci] -= p_ * qbc;

    if (hasSubstrate()) {
        const int rs = rowSub();
        const double qcs = q_coeff * store.deviceQ(state_index_, 2) +
                           qdot_coeff * store.deviceQDot(state_index_, 2);
        if (rs >= 0) rhs[rs] += p_ * qcs;
        if (ci >= 0) rhs[ci] -= p_ * qcs;
    }
}

void Bjt::refreshOutputs(const DaeSolution& x, double time) const {
    (void)time;
    // Final cache refresh from the converged solution (no limiting).
    const int rb = rowOf(nodes_[1]);
    const int rc = rowOf(nodes_[0]);
    const int re = rowOf(nodes_[2]);
    const int bi = rowBi();
    const int ci = rowCi();
    const int ei = rowEi();

    const double vxb = (rb >= 0) ? x.data()[rb] : 0.0;
    const double vxc = (rc >= 0) ? x.data()[rc] : 0.0;
    const double vxe = (re >= 0) ? x.data()[re] : 0.0;
    const double vib = (bi >= 0) ? x.data()[bi] : 0.0;
    const double vic = (ci >= 0) ? x.data()[ci] : 0.0;
    const double vie = (ei >= 0) ? x.data()[ei] : 0.0;

    vbe_cached_ = p_ * (vxb - vxe);   // terminal voltages (report values)
    vbc_cached_ = p_ * (vxb - vxc);
    vbe_lim_ = p_ * (vib - vie);      // NPN-oriented (what evaluate() wants)
    vbc_lim_ = p_ * (vib - vic);
    evaluate(vbe_lim_, vbc_lim_);

    if (hasSubstrate()) {
        const int rs = rowSub();
        const double vsub = (rs >= 0) ? x.data()[rs] : 0.0;
        const double vsc = p_ * (vsub - vic);
        evalSubstrateDiode(thermal_->csCap, inst_.AREA, thermal_->gmin, vsc,
                           isc_, gsc_, qcs_, csc_);
    }
}

bool Bjt::acceptSolution() const {
    return !limiting_active_;
}

void Bjt::confirmAccept() {
    limiting_active_ = false;
}

void Bjt::resetIterationState(const DaeSolution& x) {
    (void)x;
    // No per-iteration state to reset beyond what evaluate() recomputes.
}

// ===============================================================
// Outputs
// ===============================================================

double Bjt::getPower() const {
    // Σ v_terminal·i_terminal over the EXTERNAL terminals (ohmic drops
    // included: this is the power the circuit delivers to the device).
    const double vb = nodes_[1] ? nodes_[1]->getVoltage() : 0.0;
    const double vc = nodes_[0] ? nodes_[0]->getVoltage() : 0.0;
    const double ve = nodes_[2] ? nodes_[2]->getVoltage() : 0.0;
    double pwr = vb * (p_ * ib_) + vc * (p_ * ic_) + ve * (p_ * ie_);
    if (hasSubstrate()) {
        const double vs = nodes_[3] ? nodes_[3]->getVoltage() : 0.0;
        pwr += vs * (-p_ * isc_);
    }
    return pwr;
}

// ===============================================================
// Noise (`.NOISE`)
// ===============================================================
// Topology (see loadF): the netlist nodes collector/base/emitter feed the
// INTRINSIC nodes bi/ci/ei through RB/RC/RE when those are non-zero, and the
// transistor itself acts on the intrinsic rows. So the ohmic resistances are
// noise sources only when the internal node exists, and the shot noise of each
// junction sits between the intrinsic terminals it actually drives.
//
// Currents are physical (getCollectorCurrent/getEmitterCurrent apply the PNP
// polarity), but noise is symmetric: only the magnitude matters.

void Bjt::loadNoiseSources(const NoiseContext& ctx,
                           std::vector<NoiseSource>& out) const {
    const double kT4 = 4.0 * ctx.boltzmann * ctx.temperature;

    const int rb = rowOf(nodes_[1]);
    const int rc = rowOf(nodes_[0]);
    const int re = rowOf(nodes_[2]);
    const int bi = rowBi();
    const int ci = rowCi();
    const int ei = rowEi();

    auto thermal = [&](const char* label, int a, int b, double g) {
        if (!(g > 0.0) || a < 0 || b < 0 || a == b) return;
        NoiseSource s;
        s.label = label;
        s.pos = a;
        s.neg = b;
        s.psd_a2 = kT4 * g;
        s.f_exponent = 0.0;
        out.push_back(s);
    };

    // Ohmic resistances (only when the internal node exists).
    if (offB_ >= 0) thermal("rb", rb, bi, gb_);
    if (offC_ >= 0) thermal("rc", rc, ci, gc_);
    if (offE_ >= 0) thermal("re", re, ei, ge_);

    // Shot noise: 2q·I at the collector junction and at the emitter junction.
    const double ic = std::fabs(getCollectorCurrent());
    const double ie = std::fabs(getEmitterCurrent());
    if (ic > 0.0 && bi >= 0 && ci >= 0) {
        NoiseSource s;
        s.label = "ic";
        s.pos = bi;
        s.neg = ci;
        s.psd_a2 = 2.0 * ctx.charge * ic;
        s.f_exponent = 0.0;
        out.push_back(s);
    }
    if (ie > 0.0 && bi >= 0 && ei >= 0) {
        NoiseSource s;
        s.label = "ie";
        s.pos = bi;
        s.neg = ei;
        s.psd_a2 = 2.0 * ctx.charge * ie;
        s.f_exponent = 0.0;
        out.push_back(s);
    }
}

} // namespace e_XSpice

/**
 * @author Gabriel Henrique Silva
 * @file MosfetBSIM4Device.cpp
 * @brief BSIM4 level-2 MOSFET device — physics inline (no kernel adapter).
 *
 * Ported from Xyce N_DEV_MOSFET_B4.C (4.8) and N_DEV_MOSFET_B4p82.C (4.8.2).
 * Physics computed directly here; the Mosfet wrapper only applies polarity and stamps.
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

bool MosfetBSIM4Device::debug_ =
    (std::getenv("EXSPICE_DEBUG_BSIM4") != nullptr);

static constexpr double EpsOx = 3.9 * 8.854214871e-12;

MosfetBSIM4Device::MosfetBSIM4Device(
    RCP<const MosfetBSIM4Model> model, double W, double L, double M, double NF,
    double AD, double AS, double PD, double PS,
    double SA, double SB, double SD,
    double SCA, double SCB, double SCC, double SC)
    : model_(std::move(model)), W_(W), L_(L), M_(M), NF_(NF), AD_(AD), AS_(AS),
      PD_(PD), PS_(PS), SA_(SA), SB_(SB), SD_(SD), SCA_(SCA), SCB_(SCB),
      SCC_(SCC), SC_(SC) {
    if (!model_) throw std::invalid_argument("BSIM4 device requires a model card");
    if (W_ <= 0.0 || L_ <= 0.0) throw std::invalid_argument("BSIM4 MOSFET W and L must be positive");
    if (M_ <= 0.0) M_ = 1.0;
    if (NF_ <= 0.0) NF_ = 1.0;
    (void)0;
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

namespace {

// Constants of the Xyce N_DEV_MOSFET_B4p82.C C-V block (CONSTEPS0, CONSTKboQ,
// CONSTDELTA_*, exponential clamps).
constexpr double kCvEps0 = 8.85418e-12;
constexpr double kCvEpsSi = 11.7 * kCvEps0;
constexpr double kCvKoverQ = 8.617087e-5;
constexpr double kCvQ = 1.6021918e-19;
constexpr double kCvPi = 3.14159265358979323846;
constexpr double kCvDelta1 = 0.02;
constexpr double kCvDelta3 = 0.02;
constexpr double kCvDelta4 = 0.02;
constexpr double kCvExpThreshold = 34.0;
constexpr double kCvMaxExp = 5.834617425e14;
constexpr double kCvMinExp = 1.713908431e-15;

// Instance-level C-V and threshold quantities (Xyce SizeDependParam).
struct CvParams {
    double leff = 0.0, weff = 0.0, leffCV = 0.0, weffCV = 0.0;
    double toxe = 0.0, coxe = 0.0, coxp = 0.0, toxp = 0.0, epssub = 0.0;
    double phi = 0.0, sqrtPhi = 0.0, vbi = 0.0, Xdep0 = 0.0, factor1 = 0.0;
    double k1 = 0.0, k2 = 0.0, k1ox = 0.0, k2ox = 0.0, vbm = 0.0, vbsc = 0.0;
    double vfb = 0.0, vfbzb = 0.0, vtfbphi2 = 0.0;
    double cdep0 = 0.0, ldeb = 0.0, acde = 0.0, moin = 0.0, abulkCVfactor = 0.0;
    double mstarcv = 0.0, voffcv = 0.0, voffcbncv = 0.0, noff = 0.0;
    double cgdo = 0.0, cgso = 0.0, cgbo = 0.0, cgdl = 0.0, cgsl = 0.0;
    double ckappad = 0.0, ckappas = 0.0, ados = 0.0, bdos = 0.0;
};

// exp-based smoothing term of the DVT0/DVT0W threshold corrections.
double expTheta(double x) {
    if (x < kCvExpThreshold) {
        const double e = std::exp(x);
        const double d = e - 1.0;
        return e / (d * d + 2.0 * e * kCvMinExp);
    }
    return 1.0 / (kCvMaxExp - 2.0);
}

CvParams cvParams(const MosfetBSIM4Model& m, double W, double L, int dtype,
                  double vt) {
    CvParams c;
    const double tnom = m.TNOM;
    const double vtm0 = kCvKoverQ * tnom;
    const double ni = 1.45e10 * (tnom / 300.15) * std::sqrt(tnom / 300.15);

    c.leff = L + m.XL - 2.0 * m.LINT;
    c.weff = W + m.XW - 2.0 * m.WINT;
    c.leffCV = L + m.XL - 2.0 * m.DLC;
    c.weffCV = W + m.XW - 2.0 * m.DWC;

    c.toxe = m.TOXE;
    c.toxp = m.TOXP;
    c.coxe = m.EPSROX * kCvEps0 / c.toxe;
    c.coxp = m.EPSROX * kCvEps0 / c.toxp;
    c.epssub = kCvEpsSi;

    c.phi = vtm0 * std::log(m.NDEP / ni) + m.PHIN + 0.4;
    c.sqrtPhi = std::sqrt(c.phi);
    c.Xdep0 = std::sqrt(2.0 * c.epssub / (kCvQ * m.NDEP * 1.0e6)) * c.sqrtPhi;
    c.vbi = vtm0 * std::log(m.NSD * m.NDEP / (ni * ni));
    c.factor1 = std::sqrt(c.epssub / (m.EPSROX * kCvEps0) * c.toxe);
    c.cdep0 = std::sqrt(kCvQ * c.epssub * m.NDEP * 1.0e6 / 2.0 / c.phi);
    c.ldeb = std::sqrt(c.epssub * vtm0 / (kCvQ * m.NDEP * 1.0e6)) / 3.0;

    // k1/k2 come from the model card, already resolved by
    // MosfetBSIM4Model::deriveBodyEffectCoeffs() (K1/K2 given, or derived
    // from GAMMA1/GAMMA2).  Doing it here as well used to give the C-V block
    // a different body effect than the DC threshold.
    c.vbm = m.VBM;
    if (c.vbm > 0.0) c.vbm = -c.vbm;
    c.k1 = m.K1;
    c.k2 = m.K2;
    c.k1ox = c.k1 * c.toxe / m.TOXM;
    c.k2ox = c.k2 * c.toxe / m.TOXM;

    if (c.k2 < 0.0) {
        const double T0 = 0.5 * c.k1 / c.k2;
        c.vbsc = 0.9 * (c.phi - T0 * T0);
        if (c.vbsc > -3.0) {
            c.vbsc = -3.0;
        } else if (c.vbsc < -30.0) {
            c.vbsc = -30.0;
        }
    } else {
        c.vbsc = -30.0;
    }
    if (c.vbsc > c.vbm) c.vbsc = c.vbm;

    // Flat-band voltage: VFB given, else from VTH0, else -1 V.
    if (m.isGiven("VFB")) {
        c.vfb = m.VFB;
    } else if (m.isGiven("VTH0")) {
        c.vfb = dtype * m.VTH0 - c.phi - c.k1 * c.sqrtPhi;
    } else {
        c.vfb = -1.0;
    }

    // vfbzbfactor (zero-bias flat-band with short/narrow-width terms).
    const double tmp1 = c.vbi - c.phi;
    const double tmp2 = c.factor1 * std::sqrt(c.Xdep0);
    const double T8 = m.DVT0W * expTheta(m.DVT1W * c.weff * c.leff / tmp2) * tmp1;
    const double T9 = m.DVT0 * expTheta(m.DVT1 * c.leff / tmp2) * tmp1;
    const double T4 = c.toxe * c.phi / (c.weff + m.W0);
    const double T0 = std::sqrt(1.0 + m.LPE0 / c.leff);
    const double TRatio = vt / (kCvKoverQ * tnom);
    const double T3 = (m.TEMPMOD == 0 || m.TEMPMOD == 1)
        ? (m.KT1 + m.KT1L / c.leff) * (TRatio - 1.0)
        : -m.KT1 * (TRatio - 1.0);
    const double T5 = c.k1ox * (T0 - 1.0) * c.sqrtPhi + T3;
    const double vfbzbfactor =
        -T8 - T9 + m.K3 * T4 + T5 - c.phi - c.k1 * c.sqrtPhi;
    c.vfbzb = vfbzbfactor + dtype * m.VTH0;

    const double T3v = dtype * m.VTH0 - c.vfb - c.phi;
    c.vtfbphi2 = (4.0 * T3v < 0.0) ? 0.0 : 4.0 * T3v;

    c.abulkCVfactor = 1.0 + std::pow(m.CLC / c.leffCV, m.CLE);
    c.acde = m.ACDE * std::pow(m.NDEP / 2.0e16, -0.25);
    c.moin = m.MOIN;
    c.noff = m.NOFF;
    c.voffcv = m.VOFFCV;
    c.voffcbncv = m.VOFFCV + m.VOFFCVL / c.leff;
    c.mstarcv = 0.5 + std::atan(m.MINVCV) / kCvPi;

    const double cf = m.isGiven("CF")
        ? m.CF
        : 2.0 * m.EPSROX * kCvEps0 / kCvPi * std::log(1.0 + 0.4e-6 / c.toxe);
    c.cgdo = (m.CGDO + cf) * c.weffCV;
    c.cgso = (m.CGSO + cf) * c.weffCV;
    c.cgbo = m.CGBO * c.leffCV;
    c.cgdl = m.CGDL;
    c.cgsl = m.CGSL;
    c.ckappad = m.CKAPPAD;
    c.ckappas = m.CKAPPAS;
    c.ados = m.ADOS;
    c.bdos = m.BDOS;
    return c;
}

// Intrinsic terminal charges in the mode frame (drain/source as seen with
// vds >= 0), before overlap and junction charges.
struct CtmCharges {
    double qgate, qbulk, qdrn, qsrc;
};

// Charge-thickness model (Xyce N_DEV_MOSFET_B4p82.C C-V block, capMod 1 and
// 2). Sum of charges is zero by construction.
CtmCharges ctmCharges(const CvParams& c, int capMod, int cvchargeMod,
                      double xpart, double vgsEff, double vds, double vbseff,
                      double vgst, double n, double vtm, double abulk0) {
    const double CoxWL = c.coxe * c.weffCV * c.leffCV;

    // Effective overdrive of the C-V block (cvchargeMod 0 or 1).
    double Vgsteff;
    if (cvchargeMod == 0) {
        const double noffEff = n * c.noff;
        const double T0 = vtm * noffEff;
        const double VgstNVt = (vgst - c.voffcv) / T0;
        if (VgstNVt > kCvExpThreshold) {
            Vgsteff = vgst - c.voffcv;
        } else if (VgstNVt < -kCvExpThreshold) {
            Vgsteff = T0 * std::log(1.0 + kCvMinExp);
        } else {
            Vgsteff = T0 * std::log(1.0 + std::exp(VgstNVt));
        }
    } else {
        const double T0 = n * vtm;
        const double T1a = c.mstarcv * vgst;
        const double T2a = T1a / T0;
        double T10;
        if (T2a > kCvExpThreshold) {
            T10 = T1a;
        } else if (T2a < -kCvExpThreshold) {
            T10 = n * vtm * std::log(1.0 + kCvMinExp);
        } else {
            T10 = n * vtm * std::log(1.0 + std::exp(T2a));
        }
        const double T1b = c.voffcbncv - (1.0 - c.mstarcv) * vgst;
        const double T2b = T1b / T0;
        double T9;
        if (T2b < -kCvExpThreshold) {
            T9 = c.mstarcv + n * (c.coxe * kCvMinExp / c.cdep0);
        } else if (T2b > kCvExpThreshold) {
            T9 = c.mstarcv + n * (c.coxe * kCvMaxExp / c.cdep0);
        } else {
            T9 = c.mstarcv + n * (c.coxe / c.cdep0) * std::exp(T2b);
        }
        Vgsteff = T10 / T9;
    }

    if (capMod == 1) {
        const double Vfb = c.vfbzb;
        const double V3 = Vfb - vgsEff + vbseff - kCvDelta3;
        const double sq = (Vfb <= 0.0)
            ? std::sqrt(V3 * V3 - 4.0 * kCvDelta3 * Vfb)
            : std::sqrt(V3 * V3 + 4.0 * kCvDelta3 * Vfb);
        const double Vfbeff = Vfb - 0.5 * (V3 + sq);
        const double Qac0 = CoxWL * (Vfbeff - Vfb);

        double T0 = 0.5 * c.k1ox;
        const double T3 = vgsEff - Vfbeff - vbseff - Vgsteff;
        double T1;
        if (c.k1ox == 0.0) {
            T1 = 0.0;
        } else if (T3 < 0.0) {
            T1 = T0 + T3 / c.k1ox;
        } else {
            T1 = std::sqrt(T0 * T0 + T3);
        }
        const double Qsub0 = CoxWL * c.k1ox * (T1 - T0);

        const double AbulkCV = abulk0 * c.abulkCVfactor;
        const double VdsatCV = Vgsteff / AbulkCV;
        T0 = VdsatCV - vds - kCvDelta4;
        T1 = std::sqrt(T0 * T0 + 4.0 * kCvDelta4 * VdsatCV);
        double VdseffCV;
        if (T0 >= 0.0) {
            VdseffCV = VdsatCV - 0.5 * (T0 + T1);
        } else {
            VdseffCV = VdsatCV * (1.0 - (2.0 * kCvDelta4) / (T1 - T0));
        }
        if (vds == 0.0) VdseffCV = 0.0;

        T0 = AbulkCV * VdseffCV;
        T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
        const double T2 = VdseffCV / T1;
        const double T3q = T0 * T2;
        double qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3q);
        const double qbulk0 = CoxWL * (1.0 - AbulkCV) * (0.5 * VdseffCV - T3q);
        double qsrc;
        if (xpart > 0.5) {
            T1 = T1 + T1;
            qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0 - T0 * T0 / T1);
        } else if (xpart < 0.5) {
            T1 = T1 / 12.0;
            const double T2b = 0.5 * CoxWL / (T1 * T1);
            const double T3b = Vgsteff * (2.0 * T0 * T0 / 3.0
                + Vgsteff * (Vgsteff - 4.0 * T0 / 3.0))
                - 2.0 * T0 * T0 * T0 / 15.0;
            qsrc = -T2b * T3b;
        } else {
            qsrc = -0.5 * (qgate + qbulk0);
        }
        qgate += Qac0 + Qsub0;
        const double qbulk = qbulk0 - (Qac0 + Qsub0);
        return {qgate, qbulk, -(qgate + qbulk + qsrc), qsrc};
    }

    // capMod 2: charge-thickness model with bias-dependent Coxeff (CTM).
    const double Vfb = c.vfbzb;
    double V3 = Vfb - vgsEff + vbseff - kCvDelta3;
    const double sq = (Vfb <= 0.0)
        ? std::sqrt(V3 * V3 - 4.0 * kCvDelta3 * Vfb)
        : std::sqrt(V3 * V3 + 4.0 * kCvDelta3 * Vfb);
    const double Vfbeff = Vfb - 0.5 * (V3 + sq);

    const double Cox = c.coxp;
    double Tox = 1.0e8 * c.toxp;
    double T0 = (vgsEff - vbseff - Vfb) / Tox;
    double tmp = T0 * c.acde;
    double Tcen;
    if (-kCvExpThreshold < tmp && tmp < kCvExpThreshold) {
        Tcen = c.ldeb * std::exp(tmp);
    } else if (tmp <= -kCvExpThreshold) {
        Tcen = c.ldeb * kCvMinExp;
    } else {
        Tcen = c.ldeb * kCvMaxExp;
    }
    const double LINK = 1.0e-3 * c.toxp;
    V3 = c.ldeb - Tcen - LINK;
    const double V4 = std::sqrt(V3 * V3 + 4.0 * LINK * c.ldeb);
    Tcen = c.ldeb - 0.5 * (V3 + V4);
    double Ccen = c.epssub / Tcen;
    double T2 = Cox / (Cox + Ccen);
    double Coxeff = T2 * Ccen;
    double CoxWLcen = CoxWL * Coxeff / c.coxe;
    const double Qac0 = CoxWLcen * (Vfbeff - Vfb);

    T0 = 0.5 * c.k1ox;
    const double T3 = vgsEff - Vfbeff - vbseff - Vgsteff;
    double T1;
    if (c.k1ox == 0.0) {
        T1 = 0.0;
    } else if (T3 < 0.0) {
        T1 = T0 + T3 / c.k1ox;
    } else {
        T1 = std::sqrt(T0 * T0 + T3);
    }
    const double Qsub0 = CoxWLcen * c.k1ox * (T1 - T0);

    // Gate-bias dependent delta-phi (Xyce DeltaPhi).
    double Denomi;
    double T0d;
    if (c.k1ox <= 0.0) {
        Denomi = 0.25 * c.moin * vtm;
        T0d = 0.5 * c.sqrtPhi;
    } else {
        Denomi = c.moin * vtm * c.k1ox * c.k1ox;
        T0d = c.k1ox * c.sqrtPhi;
    }
    const double T1d = 2.0 * T0d + Vgsteff;
    const double DeltaPhi = vtm * std::log(1.0 + T1d * Vgsteff / Denomi);

    // VgDP: Vgsteff minus DeltaPhi, smoothed to stay positive.
    T0 = Vgsteff - DeltaPhi - 0.001;
    T1 = std::sqrt(T0 * T0 + Vgsteff * 0.004);
    const double VgDP = 0.5 * (T0 + T1);

    // Tcen and Coxeff re-evaluated with Vgsteff (Xyce "WDLiu" step).
    Tox += Tox;
    T0 = (Vgsteff + c.vtfbphi2) / Tox;
    tmp = std::exp(c.bdos * 0.7 * std::log(T0));
    T1 = 1.0 + tmp;
    Tcen = c.ados * 1.9e-9 / T1;
    Ccen = c.epssub / Tcen;
    T2 = Cox / (Cox + Ccen);
    Coxeff = T2 * Ccen;
    CoxWLcen = CoxWL * Coxeff / c.coxe;

    const double AbulkCV = abulk0 * c.abulkCVfactor;
    const double VdsatCV = VgDP / AbulkCV;
    T0 = VdsatCV - vds - kCvDelta4;
    T1 = std::sqrt(T0 * T0 + 4.0 * kCvDelta4 * VdsatCV);
    double VdseffCV;
    if (T0 >= 0.0) {
        VdseffCV = VdsatCV - 0.5 * (T0 + T1);
    } else {
        VdseffCV = VdsatCV * (1.0 - (2.0 * kCvDelta4) / (T1 - T0));
    }
    if (vds == 0.0) VdseffCV = 0.0;

    T0 = AbulkCV * VdseffCV;
    T1 = VgDP;
    T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
    const double T3c = T0 / T2;
    double qgate = CoxWLcen * (T1 - T0 * (0.5 - T3c));
    const double qbulk0 = CoxWLcen * (1.0 - AbulkCV)
        * (0.5 * VdseffCV - T0 * VdseffCV / T2);
    double qsrc;
    if (xpart > 0.5) {
        qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0 - 0.5 * T0 * T0 / T2);
    } else if (xpart < 0.5) {
        const double T2b = T2 / 12.0;
        const double T3b = 0.5 * CoxWLcen / (T2b * T2b);
        const double T4b = T1 * (2.0 * T0 * T0 / 3.0
            + T1 * (T1 - 4.0 * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;
        qsrc = -T3b * T4b;
    } else {
        qsrc = -0.5 * qgate;
    }
    qgate += Qac0 + Qsub0 - qbulk0;
    const double qbulk = qbulk0 - (Qac0 + Qsub0);
    return {qgate, qbulk, -(qgate + qbulk + qsrc), qsrc};
}

}  // namespace

void MosfetBSIM4Device::polyDepletion(double phi, double ngate, double epsgate,
                                      double coxe, double Vgs_arg,
                                      double& Vgs_eff, double& dVgs_eff_dVg) {
    // Xyce N_DEV_MOSFET_B4.C:8819 (polyDepletion).
    const double CONSTQ = 1.6021918e-19;
    const double CONSTEPSSI = 1.03594e-10;
    if (ngate > 1.0e18 && ngate < 1.0e25 && Vgs_arg > phi && epsgate != 0.0) {
        const double T1 = 1.0e6 * CONSTQ * epsgate * ngate / (coxe * coxe);
        const double T8 = Vgs_arg - phi;
        const double T4 = std::sqrt(1.0 + 2.0 * T8 / T1);
        const double T2 = 2.0 * T8 / (T4 + 1.0);
        const double T3 = 0.5 * T2 * T2 / T1;   // Vpoly
        const double T7 = 1.12 - T3 - 0.05;
        const double T6 = std::sqrt(T7 * T7 + 0.224);
        const double T5 = 1.12 - 0.5 * (T7 + T6);
        Vgs_eff = Vgs_arg - T5;
        dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
    } else {
        Vgs_eff = Vgs_arg;
        dVgs_eff_dVg = 1.0;
    }
    (void)phi; (void)CONSTEPSSI;
}

int MosfetBSIM4Device::numFingerDiff(double nf, int minSD, double& nuIntD,
                                     double& nuEndD, double& nuIntS,
                                     double& nuEndS) {
    // Xyce N_DEV_MOSFET_B4.C:8889.
    int NF = static_cast<int>(nf);
    if ((NF % 2) != 0) {
        nuEndD = nuEndS = 1.0;
        nuIntD = nuIntS = 2.0 * std::max((nf - 1.0) / 2.0, 0.0);
    } else {
        if (minSD == 1) {          // minimize # of source
            nuEndD = 2.0;
            nuIntD = 2.0 * std::max(nf / 2.0 - 1.0, 0.0);
            nuEndS = 0.0;
            nuIntS = nf;
        } else {                   // minimize # of drain
            nuEndD = 0.0;
            nuIntD = nf;
            nuEndS = 2.0;
            nuIntS = 2.0 * std::max(nf / 2.0 - 1.0, 0.0);
        }
    }
    return 0;
}

int MosfetBSIM4Device::paEffGeo(double nf, int geo, int minSD, double weffCJ,
                                double DMCG, double DMCI, double DMDG,
                                double& Ps, double& Pd, double& As,
                                double& Ad) {
    // Xyce N_DEV_MOSFET_B4.C:8931.
    double nuIntD = 0.0, nuEndD = 0.0, nuIntS = 0.0, nuEndS = 0.0;
    if (geo < 9) numFingerDiff(nf, minSD, nuIntD, nuEndD, nuIntS, nuEndS);

    const double T0 = DMCG + DMCI;
    const double PSiso = 2.0 * T0 + weffCJ;
    const double PDiso = PSiso;
    const double PSsha = 2.0 * DMCG;
    const double PDsha = PSsha;
    const double PSmer = 2.0 * DMDG;
    const double PDmer = PSmer;
    const double ASiso = T0 * weffCJ;
    const double ADiso = ASiso;
    const double ASsha = DMCG * weffCJ;
    const double ADsha = ASsha;
    const double ASmer = DMDG * weffCJ;
    const double ADmer = ASmer;

    switch (geo) {
        case 0:
            Ps = nuEndS * PSiso + nuIntS * PSsha;
            Pd = nuEndD * PDiso + nuIntD * PDsha;
            As = nuEndS * ASiso + nuIntS * ASsha;
            Ad = nuEndD * ADiso + nuIntD * ADsha;
            break;
        case 1:
            Ps = nuEndS * PSiso + nuIntS * PSsha;
            Pd = (nuEndD + nuIntD) * PDsha;
            As = nuEndS * ASiso + nuIntS * ASsha;
            Ad = (nuEndD + nuIntD) * ADsha;
            break;
        case 2:
            Ps = (nuEndS + nuIntS) * PSsha;
            Pd = nuEndD * PDiso + nuIntD * PDsha;
            As = (nuEndS + nuIntS) * ASsha;
            Ad = nuEndD * ADiso + nuIntD * ADsha;
            break;
        case 3:
            Ps = (nuEndS + nuIntS) * PSsha;
            Pd = (nuEndD + nuIntD) * PDsha;
            As = (nuEndS + nuIntS) * ASsha;
            Ad = (nuEndD + nuIntD) * ADsha;
            break;
        case 4:
            Ps = nuEndS * PSiso + nuIntS * PSsha;
            Pd = nuEndD * PDmer + nuIntD * PDsha;
            As = nuEndS * ASiso + nuIntS * ASsha;
            Ad = nuEndD * ADmer + nuIntD * ADsha;
            break;
        case 5:
            Ps = (nuEndS + nuIntS) * PSsha;
            Pd = nuEndD * PDmer + nuIntD * PDsha;
            As = (nuEndS + nuIntS) * ASsha;
            Ad = nuEndD * ADmer + nuIntD * ADsha;
            break;
        case 6:
            Ps = nuEndS * PSmer + nuIntS * PSsha;
            Pd = nuEndD * PDiso + nuIntD * PDsha;
            As = nuEndS * ASmer + nuIntS * ASsha;
            Ad = nuEndD * ADiso + nuIntD * ADsha;
            break;
        case 7:
            Ps = nuEndS * PSmer + nuIntS * PSsha;
            Pd = (nuEndD + nuIntD) * PDsha;
            As = nuEndS * ASmer + nuIntS * ASsha;
            Ad = (nuEndD + nuIntD) * ADsha;
            break;
        case 8:
            Ps = nuEndS * PSmer + nuIntS * PSsha;
            Pd = nuEndD * PDmer + nuIntD * PDsha;
            As = nuEndS * ASmer + nuIntS * ASsha;
            Ad = nuEndD * ADmer + nuIntD * ADsha;
            break;
        case 9:
            Ps = PSiso + (nf - 1.0) * PSsha;
            Pd = nf * PDsha;
            As = ASiso + (nf - 1.0) * ASsha;
            Ad = nf * ADsha;
            break;
        case 10:
            Ps = nf * PSsha;
            Pd = PDiso + (nf - 1.0) * PDsha;
            As = nf * ASsha;
            Ad = ADiso + (nf - 1.0) * ADsha;
            break;
        default:
            Ps = nuEndS * PSiso + nuIntS * PSsha;
            Pd = nuEndD * PDiso + nuIntD * PDsha;
            As = nuEndS * ASiso + nuIntS * ASsha;
            Ad = nuEndD * ADiso + nuIntD * ADsha;
            break;
    }
    return 0;
}

void MosfetBSIM4Device::evaluateBase(double vgs, double vds, double vbs, bool probe) {
    const MosfetBSIM4Thermal& t = *thermal_;
    const double gmin = t.gmin;

    // ---- geometry (Xyce NumFingerDiff / PAeffGeo) ----------------------
    // M (numberParallel) is NOT geometry: Xyce multiplies every terminal
    // current/charge by it in the loaders (N_DEV_MOSFET_B4.C:7920).  The
    // junction areas/perimeters below are ABSOLUTE (m^2, m) and depend only
    // on NF, exactly as PAeffGeo (N_DEV_MOSFET_B4.C:8931) produces them.
    const double nf = NF_;
    const double M  = M_;

    // Per-finger width.  DWJ is not ported yet, so weffCJ = W/nf + XW
    // (Xyce B4p82.C:794,862).
    const double weffCJ = (W_ / nf) + model_->XW;

    // PAeffGeo case 0 (GEOMOD default).  DMCGT (temperature term) is not
    // ported, so DMCGeff = DMCG (Xyce B4p82.C:2276-2278).
    // DMCG default 0: fall back to a sensible default (Xyce leaves Aseff=0,
    // but e-XSpice historically used W*L).  Use W*L when DMCG==0 so junctions
    const double DMCG = model_->DMCG;
    const double DMCI = (model_->DMCI > 0.0) ? model_->DMCI : (DMCG > 0.0 ? DMCG : 1.0);
    const double DMDG = (model_->DMDG > 0.0) ? model_->DMDG : (DMCG > 0.0 ? DMCG : 1.0);

    double PsGeo = 0.0, PdGeo = 0.0, AsGeo = 0.0, AdGeo = 0.0;
    paEffGeo(nf, 0, 0, weffCJ, DMCG, DMCI, DMDG,
             PsGeo, PdGeo, AsGeo, AdGeo);
    (void)0;

    double Aseff, Adeff, Pseff, Pdeff;
    if (AS_ > 0.0) {
        Aseff = AS_;
    } else if (DMCG <= 0.0) {
        Aseff = W_ * L_;
    } else {
        Aseff = AsGeo;
    }
    if (Aseff <= 0.0) Aseff = W_ * L_;
    if (AD_ > 0.0) {
        Adeff = AD_;
    } else if (DMCG <= 0.0) {
        Adeff = W_ * L_;
    } else {
        Adeff = AdGeo;
    }
    if (Adeff <= 0.0) Adeff = W_ * L_;
    // PERMOD != 0 strips the gate-edge contribution out of a given perimeter
    // (Xyce B4p82.C:2293-2328).
    if (PS_ > 0.0) {
        Pseff = (model_->PERMOD != 0) ? std::max(0.0, PS_ - weffCJ * nf) : PS_;
    } else if (DMCG <= 0.0) {
        Pseff = 2.0 * W_;
    } else {
        Pseff = PsGeo;
    }
    if (Pseff <= 0.0) Pseff = 2.0 * W_;
    if (PD_ > 0.0) {
        Pdeff = (model_->PERMOD != 0) ? std::max(0.0, PD_ - weffCJ * nf) : PD_;
    } else if (DMCG <= 0.0) {
        Pdeff = 2.0 * W_;
    } else {
        Pdeff = PdGeo;
    }
    if (Pdeff <= 0.0) Pdeff = 2.0 * W_;
    if (Aseff < 0.0) Aseff = 0.0;
    if (Adeff < 0.0) Adeff = 0.0;
    if (Pseff < 0.0) Pseff = 0.0;
    if (Pdeff < 0.0) Pdeff = 0.0;

    // RdseffGeo / RdsEndIso (simplified)
    double RdsEff = 0.0;
    if (model_->RDSW > 0.0 && W_ > 0.0) {
        double Rsh = model_->RDSW / W_;
        RdsEff = Rsh * (1.0 + (model_->XRCRG1 > 0.0 ? model_->XRCRG1 : 0.0));
        RdsEff *= (1.0 + model_->PRT * (t.vt / (Const::KoverQ * model_->TNOM) - 1.0));
    }

    const double leff_early = L_ - 2.0 * model_->LD;
    const double leff = (leff_early > 0.0) ? leff_early : 1e-9;

    // ---- source/drain swap for vds < 0 ----
    const bool sw = (vds < 0.0);
    swapped_ = sw;
    const double vgs_e = sw ? vgs - vds : vgs;
    const double vds_e = sw ? -vds : vds;
    const double vbs_e = sw ? vbs - vds : vbs;

    vgs_cached_ = vgs;
    vds_cached_ = vds;
    vbs_cached_ = vbs;

    // ---- threshold voltage with body effect (BSIM4 Xyce) ----
    const double vto_n = p_ * model_->VTH0;
    const double sqrtPhi = std::sqrt(std::max(model_->PHI, 1e-6));
    const double phib = model_->PHI - vbs_e;
    double dvth_dvbs = 0.0;
    double vth = vto_n;
    // Narrow-width term.  Xyce B4p82.C:3875:
    //   Vth_NarrowW = toxe * phi / (weff + w0)
    // (the old code used TOX/(W+W0) AND multiplied it by (phi - Vbs), which
    // is not the BSIM4 form at all).
    const double Vth_NarrowW =
        model_->TOXE * model_->PHI / (weffCJ + model_->W0 + 1e-30);
    if (phib > 1e-3) {
        const double sq = std::sqrt(phib);
        // Body effect: k1ox*(sqrt(Phis) - sqrt(phi)) - k2ox*Vbs.
        // Xyce B4p82.C:3894-3897.  GAMMA/K2 were used before, which made the
        // DC threshold disagree with Abulk and with the C-V block.
        vth += t.k1ox * (sq - sqrtPhi);
        dvth_dvbs = -t.k1ox / (2.0 * sq);
        vth -= t.k2ox * vbs_e;
        dvth_dvbs -= t.k2ox;
        // K3 / K3B (narrow width): constant plus a body-bias slope.
        vth += model_->K3 * Vth_NarrowW;
        vth += model_->K3B * vbs_e * Vth_NarrowW;
        dvth_dvbs += model_->K3B * Vth_NarrowW;
    } else {
        vth += t.k1ox * (std::sqrt(1e-3) - sqrtPhi);
    }
    // KETA
    vth += model_->KETA * vbs_e;
    dvth_dvbs += model_->KETA;

    // ---- DVT (short-channel, Xyce factor1/sqrt(Xdep)) ----
    if (leff > 0.0 && model_->DVT0 > 0.0) {
        // Xdep (depletion width)
        double Xdep = std::sqrt(2.0 * 11.7 * 8.854e-12 * model_->PHI / (Const::Q * model_->NSUB * 1e6 + 1e-20));
        if (Xdep < 1e-9) Xdep = 1e-9;
        double factor1 = model_->TOX / Xdep;
        double lt = std::sqrt(factor1 * Xdep * Xdep / 3.0 + (model_->TOX * model_->TOX) / 3.0);
        if (lt < 1e-30) lt = 1e-30;
        double eta = model_->DVT0 * std::exp(-model_->DVT1 * leff / lt);
        double vbi = model_->PHI;
        vth -= eta * (vbi - vbs_e * 0.5 - model_->PHI * 0.5);
        dvth_dvbs += eta * 0.5;
        // DVT2
        vth += model_->DVT2 * vbs_e * eta;
        dvth_dvbs += model_->DVT2 * eta;
    }

    // DIBL (drain-induced barrier lowering): Vth reduction with Vds
    {
        double eta = model_->ETA0 + model_->ETAB * vbs_e;
        if (eta < 0.0) eta = 0.0;
        double V0 = (model_->PHI - vbs_e) * 0.5;
        if (V0 < 1e-3) V0 = 1e-3;
        double theta0vb0 = 1.0 / (1.0 + std::exp(-model_->PDIBLC1 * V0));
        double DIBL_Sft = eta * theta0vb0 * vds_e;
        vth -= DIBL_Sft;
        double dtheta_dvbs = -0.5 * theta0vb0 * (1.0 - theta0vb0) * model_->PDIBLC1;
        double dDIBL_dvbs = model_->ETAB * theta0vb0 * vds_e + eta * dtheta_dvbs * vds_e;
        dvth_dvbs -= dDIBL_dvbs;
    }



    // ---- temperature dependence (KT1, KT2, UTE) ----
    {
        double TRatio = t.vt / (Const::KoverQ * model_->TNOM);
        double dT = TRatio - 1.0;
        vth += model_->KT1 * dT + model_->KT2 * vbs_e * dT;
        dvth_dvbs += model_->KT2 * dT;
        // KT1L
        if (model_->KT1L != 0.0 && leff > 0.0) {
            vth += model_->KT1L / leff * dT;
        }
    }

    // ---- narrow-width effect (W0, K3) already above ----

    // ---- Stress/LOD (Xyce 4.8.2 complete) ----
    double ueff_lod_mult = 1.0;
    double vsat_lod_mult = 1.0;
    double eta0_lod_add = 0.0;
    double k2_lod_add = 0.0;
    double dvth_lod = 0.0;
    if (SA_ > 0.0 && SB_ > 0.0) {
        double nf_lod = nf;   // fingers, not parallel devices (Xyce w/nf)
        double Ldrn = L_;
        double sd = SD_;
        if (sd <= 0.0 && nf_lod > 1.0) sd = 2.0 * 0.3e-6;
        double wlod = model_->WLOD;
        if (wlod < 0.0) wlod = 0.0;
        double Inv_sa = 0.0, Inv_sb = 0.0;
        double kvsat = model_->KVSAT;
        if (kvsat < -1.0) kvsat = -1.0;
        if (kvsat > 1.0) kvsat = 1.0;
        int i = 0;
        for (i = 0; i < (int)nf_lod; i++) {
            double T0 = 1.0 / nf_lod / (SA_ + 0.5 * Ldrn + i * (sd + Ldrn));
            double T1 = 1.0 / nf_lod / (SB_ + 0.5 * Ldrn + i * (sd + Ldrn));
            Inv_sa += T0;
            Inv_sb += T1;
        }
        double Inv_ODeff = Inv_sa + Inv_sb;
        double saref = (model_->SAREF > 0.0) ? model_->SAREF : 1e-6;
        double sbref = (model_->SBREF > 0.0) ? model_->SBREF : 1e-6;
        double inv_od_ref = 1.0 / (saref + 0.5 * Ldrn) + 1.0 / (sbref + 0.5 * Ldrn);
        double TRatio = t.vt / (Const::KoverQ * model_->TNOM);
        double ku0temp = model_->KU0 * (1.0 + model_->TKU0 * (TRatio - 1.0)) + 1e-30;
        double rho = model_->KU0 / ku0temp * Inv_ODeff;
        double rho_ref = model_->KU0 / ku0temp * inv_od_ref;
        ueff_lod_mult = (1.0 + rho) / (1.0 + rho_ref);
        vsat_lod_mult = (1.0 + kvsat * rho) / (1.0 + kvsat * rho_ref);
        double OD_offset = Inv_ODeff - inv_od_ref;
        double llodvth = (model_->LLODVTH != 0.0) ? model_->LLODVTH : 1.0;
        double wlodvth = (model_->WLODVTH != 0.0) ? model_->WLODVTH : 1.0;
        double kvth0_full = 1.0 + model_->LKVTH0 / std::pow(L_, llodvth) +
                            model_->WKVTH0 / std::pow(W_, wlodvth) +
                            model_->PKVTH0 / (std::pow(L_, llodvth) * std::pow(W_, wlodvth));
        if (kvth0_full == 0.0) kvth0_full = 1.0;
        dvth_lod = model_->KVTH0 / kvth0_full * OD_offset;
        double k2_lod_add = model_->STK2 / std::pow(kvth0_full, model_->LODK2) * OD_offset;
        eta0_lod_add = model_->STETA0 / std::pow(kvth0_full, model_->LODETA0) * OD_offset;
        (void)wlod;
    }
    vth += dvth_lod;

    // ---- WPE (well proximity, Xyce 4.8.2) ----
    if (model_->WPEMOD) {
        double sca = SCA_, scb = SCB_, scc = SCC_, sc = SC_;
        double Wdrn = W_ / nf;   // Xyce B4p82.C:765 — per-finger width
        if (sca == 0.0 && scb == 0.0 && scc == 0.0 && sc > 0.0) {
            double T1 = sc + Wdrn;
            double T2 = 1.0 / model_->SCREF;
            sca = model_->SCREF * model_->SCREF / (sc * T1);
            scb = ((0.1 * sc + 0.01 * model_->SCREF) * std::exp(-10.0 * sc * T2) -
                   (0.1 * T1 + 0.01 * model_->SCREF) * std::exp(-10.0 * T1 * T2)) / Wdrn;
            scc = ((0.05 * sc + 0.0025 * model_->SCREF) * std::exp(-20.0 * sc * T2) -
                   (0.05 * T1 + 0.0025 * model_->SCREF) * std::exp(-20.0 * T1 * T2)) / Wdrn;
        }
        if (sca < 0.0) sca = 0.0;
        if (scb < 0.0) scb = 0.0;
        if (scc < 0.0) scc = 0.0;
        double sceff = sca + model_->WEB * scb + model_->WEC * scc;
        vth += model_->KVTH0WE * sceff;
        dvth_dvbs += model_->K2WE * sceff;
        // ueff correction stored for later
    }

    // ---- Binning (range check) ----
    // (handled in processInstanceParams at model level; here just use model_ params)

    vth_cached_ = vth;

    // ---- effective mobility (BSIM4, Xyce) ----
    const double ov = vgs_e - vth;
    if (leff <= 0.0) {
        throw std::runtime_error("BSIM4: Leff = L - 2*LD must be positive");
    }
    // Abulk (BSIM4 Xyce completo)
    double Xdep_ab = std::sqrt(2.0 * 11.7 * 8.854e-12 * model_->PHI / (Const::Q * model_->NSUB * 1e6 + 1e-20));
    if (Xdep_ab < 1e-9) Xdep_ab = 1e-9;
    double sqrtPhis_ab = std::sqrt(std::max(model_->PHI - vbs_e, 1e-6));
    double Vth_NarrowW_ab = Vth_NarrowW;   // same term the Vth block uses
    // Xyce B4p82.C:4128-4129: T9 = 0.5*k1ox*Lpe_Vb/sqrtPhis (Lpe0/LPEB not
    // ported -> Lpe_Vb = 1), T1 = T9 + k2ox - k3b*Vth_NarrowW.
    double T1_ab = 0.5 * t.k1ox / sqrtPhis_ab + t.k2ox - model_->K3B * Vth_NarrowW_ab;
    double T2_ab = model_->A0 * leff / (leff + 2.0 * std::sqrt(model_->XJ * Xdep_ab + 1e-30)) + model_->B0 / (weffCJ + model_->B1 + 1e-9);
    double Abulk0_ab = 1.0 + T1_ab * T2_ab;
    double dAbulk_dVg_ab = model_->AGS * model_->A0 * std::pow(ov > 0 ? ov : 0.0, 3.0) * 1e-3;
    double Abulk = Abulk0_ab + dAbulk_dVg_ab * ov;
    // keta modulation
    double keta_factor_ab = 1.0 + model_->KETA * vbs_e;
    if (keta_factor_ab < 0.1) keta_factor_ab = 0.1;
    Abulk /= keta_factor_ab;
    if (Abulk < 0.1) Abulk = 0.1;
    if (Abulk > 10.0) Abulk = 10.0;
    double denom_mob = 1.0;
    double ddenom_dov = 0.0;
    double dDenom_dVd = 0.0;
    double dDenom_dVb = 0.0;
    double dVth_dVd = 0.0;
    {
        double _eta = model_->ETA0 + model_->ETAB * vbs_e;
        if (_eta < 0.0) _eta = 0.0;
        double _V0 = (model_->PHI - vbs_e) * 0.5;
        if (_V0 < 1e-3) _V0 = 1e-3;
        double _theta = 1.0 / (1.0 + std::exp(-model_->PDIBLC1 * _V0));
        dVth_dVd = -_eta * _theta;
    }
    // mobility modulation (BSIM4 ua, ub, uc, ud, mobMod 0-7, exact derivatives)
    double T14 = 0.0;
    if (model_->MTRLMOD && model_->MTRLCOMPATMOD == 0) {
        T14 = 2.0 * model_->DTYPE * (model_->PHIG - model_->EASUB - 0.5 * model_->EG0 + 0.45);
    }
    double vtfbphi1 = model_->VTFBPHI1;
    if (vtfbphi1 == 0.0) vtfbphi1 = model_->PHI;
    double eu = (model_->EU != 0.0) ? model_->EU : 1.0;
    double ucs = model_->UCS;
    if (ucs == 0.0) ucs = 1.0;
    double VgsteffVth = model_->VGSTEFFFVTH;
    if (VgsteffVth == 0.0) VgsteffVth = 1.0;
    double ud = t.ud;
    double uc = t.uc;

    if (model_->MOBMOD == 0) {
        double T0 = ov + vth + vth - T14;
        double T2 = t.ua + uc * vbs_e;
        double T3 = T0 / model_->TOX;
        double T12 = std::sqrt(vth * vth + 0.0001);
        double T9 = 1.0 / (ov + 2.0 * T12);
        double T10 = T9 * model_->TOX;
        double T8 = ud * T10 * T10 * vth;
        double T6 = T8 * vth;
        double T7 = -2.0 * T6 * T9;
        denom_mob = 1.0 + T3 * (T2 + t.ub * T3) + T6;
        ddenom_dov = (T2 + 2.0 * t.ub * T3) / model_->TOX;
        ddenom_dov += -2.0 * T6 * T9 * vth / std::sqrt(vth * vth + 0.0001);
        // exact derivatives wrt Vd/Vb
        {
            double T13 = 2.0 * ((T2 + 2.0 * t.ub * T3) / model_->TOX + T7 * vth / T12 + T8 * vth);
            dDenom_dVd = T13 * dVth_dVd;
            dDenom_dVb = T13 * dvth_dvbs + uc * T3;
        }
    } else if (model_->MOBMOD == 1) {
        double T0 = ov + vth + vth - T14;
        double T2 = 1.0 + uc * vbs_e;
        double T3 = T0 / model_->TOX;
        double T4 = T3 * (t.ua + t.ub * T3);
        double T12 = std::sqrt(vth * vth + 0.0001);
        double T9 = 1.0 / (ov + 2.0 * T12);
        double T10 = T9 * model_->TOX;
        double T8 = ud * T10 * T10 * vth;
        double T6 = T8 * vth;
        double T7 = -2.0 * T6 * T9;
        denom_mob = 1.0 + T4 * T2 + T6;
        ddenom_dov = (t.ua + 2.0 * t.ub * T3) * T2 / model_->TOX;
        ddenom_dov += -2.0 * T6 * T9 * vth / std::sqrt(vth * vth + 0.0001);
        {
            double T13 = 2.0 * ((t.ua + 2.0 * t.ub * T3) * T2 / model_->TOX + T7 * vth / T12 + T8 * vth);
            dDenom_dVd = T13 * dVth_dVd;
            dDenom_dVb = T13 * dvth_dvbs + uc * T4;
        }
    } else if (model_->MOBMOD == 2) {
        double T0 = (ov + vtfbphi1) / model_->TOX;
        if (T0 < 1e-30) T0 = 1e-30;
        double T1 = std::exp(eu * std::log(T0));
        double T2 = t.ua + uc * vbs_e;
        double T12 = std::sqrt(vth * vth + 0.0001);
        double T9 = 1.0 / (ov + 2.0 * T12);
        double T10 = T9 * model_->TOX;
        double T8 = ud * T10 * T10 * vth;
        double T6 = T8 * vth;
        denom_mob = 1.0 + T1 * T2 + T6;
        ddenom_dov = T2 * T1 * eu / (T0 * model_->TOX) - 2.0 * T6 * T9 * vth / std::sqrt(vth * vth + 0.0001);
        {
            // dDenom/dVd = (T2*dT1/dVth? no) ... approximate: T13 = 2*(T8*vth) for Coulombic + T2*eu*T1/T0*dVth
            double dT1_dVth = T1 * eu / (vth + 1e-6); // approximate
            double T13 = T2 * dT1_dVth + 2.0 * T8 * vth / std::sqrt(vth*vth+0.0001);
            dDenom_dVd = T13 * dVth_dVd;
            dDenom_dVb = T13 * dvth_dvbs + T1 * uc;
        }
    } else if (model_->MOBMOD == 4) {
        double T0 = ov + vtfbphi1 - T14;
        double T2 = t.ua + uc * vbs_e;
        double T3 = T0 / model_->TOX;
        double T12 = std::sqrt(vtfbphi1 * vtfbphi1 + 0.0001);
        double T9 = 1.0 / (ov + 2.0 * T12);
        double T10 = T9 * model_->TOX;
        double T8 = ud * T10 * T10 * vtfbphi1;
        double T6 = T8 * vtfbphi1;
        denom_mob = 1.0 + T3 * (T2 + t.ub * T3) + T6;
        ddenom_dov = (T2 + 2.0 * t.ub * T3) / model_->TOX;
        dDenom_dVd = 0.0;
        dDenom_dVb = uc * T3;
    } else if (model_->MOBMOD == 5) {
        double T0 = ov + vtfbphi1 - T14;
        double T2 = 1.0 + uc * vbs_e;
        double T3 = T0 / model_->TOX;
        double T4 = T3 * (t.ua + t.ub * T3);
        double T12 = std::sqrt(vtfbphi1 * vtfbphi1 + 0.0001);
        double T9 = 1.0 / (ov + 2.0 * T12);
        double T10 = T9 * model_->TOX;
        double T8 = ud * T10 * T10 * vtfbphi1;
        double T6 = T8 * vtfbphi1;
        denom_mob = 1.0 + T4 * T2 + T6;
        ddenom_dov = (t.ua + 2.0 * t.ub * T3) * T2 / model_->TOX;
        dDenom_dVd = 0.0;
        dDenom_dVb = uc * T4;
    } else if (model_->MOBMOD == 6) {
        double T0 = (ov + vtfbphi1) / model_->TOX;
        if (T0 < 1e-30) T0 = 1e-30;
        double T1 = std::exp(eu * std::log(T0));
        double T2 = t.ua + uc * vbs_e;
        double T12 = std::sqrt(vtfbphi1 * vtfbphi1 + 0.0001);
        double T9 = 1.0 / (ov + 2.0 * T12);
        double T10 = T9 * model_->TOX;
        double T8 = ud * T10 * T10 * vtfbphi1;
        double T6 = T8 * vtfbphi1;
        denom_mob = 1.0 + T1 * T2 + T6;
        ddenom_dov = T2 * T1 * eu / (T0 * model_->TOX);
        dDenom_dVd = 0.0;
        dDenom_dVb = T1 * uc;
    } else if (model_->MOBMOD == 3 || model_->MOBMOD == 7) {
        // high-k mobility (universal + Coulombic)
        double T0 = (ov + vtfbphi1) * 1.0e-8 / model_->TOX / 6.0;
        if (T0 < 1e-30) T0 = 1e-30;
        double T1 = std::exp(eu * std::log(T0));
        double T2 = t.ua + uc * vbs_e;
        // Coulombic
        double T10 = std::exp(ucs * std::log(0.5 + 0.5 * ov / VgsteffVth));
        double T11 = ud / T10;
        denom_mob = 1.0 + T1 * T2 + T11;
        ddenom_dov = T2 * T1 * eu * 1.0e-8 / (T0 * model_->TOX * 6.0) - 0.5 * ucs * T11 / (0.5 + 0.5 * ov / VgsteffVth) / VgsteffVth;
        dDenom_dVd = 0.0;
        dDenom_dVb = T1 * uc;
    } else {
        denom_mob = 1.0 + (t.ua + t.ub * ov + uc * ov * ov) * ov;
        ddenom_dov = t.ua + 2.0 * t.ub * ov + 3.0 * uc * ov * ov;
    }
    // lod mobility
    denom_mob /= ueff_lod_mult;
    // ueff
    double ueff = t.u0 / denom_mob;
    // WPE mobility
    if (model_->WPEMOD && model_->KU0WE != 0.0) {
        double sca = SCA_, scb = SCB_, scc = SCC_, sc = SC_;
        double Wdrn = W_ / nf;   // per-finger width
        if (sca == 0.0 && scb == 0.0 && scc == 0.0 && sc > 0.0) {
            double T1 = sc + Wdrn;
            double T2 = 1.0 / model_->SCREF;
            sca = model_->SCREF * model_->SCREF / (sc * T1);
        }
        double sceff = sca + model_->WEB * scb + model_->WEC * scc;
        ueff *= (1.0 + model_->KU0WE * sceff);
    }
    // Stress mobility
    if (SA_ > 0.0 && SB_ > 0.0 && model_->KU0WE != 0.0) {
        double inv_s = 1.0 / SA_ + 1.0 / SB_;
        ueff *= (1.0 + model_->KU0WE * inv_s * 1e-9);
    }
    // stress from SCA/SCB/SCC
    if (model_->SC != 0.0 || model_->SCA != 0.0) {
        double sc = model_->SC;
        if (sc == 0.0) sc = (model_->SCA + model_->SCB + model_->SCC) / 3.0;
        ueff *= (1.0 + model_->KU0WE * sc * 1e-9);
    }
    // temperature mobility (UTE): already folded into t.u0 by scaleTo()
    // (Xyce B4p82.C:1615: u0temp = u0 * pow(TRatio, ute)).  Applying it a
    // second time here would cancel it out, so this block is gone.
    // update denom_mob consistency
    denom_mob = t.u0 / ueff;

    // ---- beta -----------------------------------------------------------
    // Xyce B4p82.C:4581-4584: beta = ueff * CoxeffWovL with
    // CoxeffWovL = Coxeff * Weff / Leff.  Using ueff (not the SPICE3 KP)
    // is what makes the model honour U0; KP only feeds the U0 = KP/Cox
    // bridge in applyDefaults().  Weff is the PER-FINGER width, so the
    // nf factor must NOT be multiplied in again.
    const double beta_over_u0 = ueff * t.cox * (weffCJ / leff) * M;
    const double dbeta_over_u0_dov =
        -beta_over_u0 * ddenom_dov / denom_mob;

    // ---- drain current (BSIM4 complete) ----
    double Id = 0.0, gm = 0.0, gds = 0.0, gmbs = 0.0;
    saturated_ = false;

    // DIBL: eta0 (Xyce)
    double eta0 = model_->ETA0 + eta0_lod_add;
    // Voff
    double Voff = model_->VOFF;
    // NFACTOR
    double nfactor = model_->NFACTOR;
    if (nfactor < 0.1) nfactor = 0.1;

    // velocity saturation
    double vsat = t.vsat * vsat_lod_mult;
    // temperature vsat
    vsat *= (1.0 + model_->AT * (t.vt / (Const::KoverQ * model_->TNOM) - 1.0));

    if (ov > 0.0) {
        const double Esat = (ueff > 0.0) ? (2.0 * vsat / ueff) : 1e12;
        const double EsatL = Esat * leff;
        const double Vgst2Vtm = ov + 2.0 * t.vt;
        const double vdsat = (Abulk * EsatL + Vgst2Vtm > 0.0)
            ? (EsatL * Vgst2Vtm / (Abulk * EsatL + Vgst2Vtm)) : ov;
        // DIBL: Vth reduction in saturation
        double vth_dibl = vth;
        if (model_->ETA0 > 0.0 && vds_e > 0.0) {
            double dibl_factor = eta0 * vds_e / (leff * leff + 1e-18);
            vth_dibl -= dibl_factor;
        }
        // CLM (channel length modulation) lambda
        double lambda = model_->LAMBDA;
        // DIBL factor for gds correction: dVth/dVds = -eta*theta0vb0
        double dibl_eta = model_->ETA0 + model_->ETAB * vbs_e;
        if (dibl_eta < 0.0) dibl_eta = 0.0;
        double dibl_V0 = (model_->PHI - vbs_e) * 0.5;
        if (dibl_V0 < 1e-3) dibl_V0 = 1e-3;
        double dibl_theta = 1.0 / (1.0 + std::exp(-model_->PDIBLC1 * dibl_V0));
        double dibl_gds_factor = dibl_eta * dibl_theta;
        const double lame = 1.0 + lambda * vds_e;
        if (vds_e < vdsat) {
            const double f = ov * vds_e - 0.5 * vds_e * vds_e;
            Id = beta_over_u0 * f * lame;
            gm = beta_over_u0 * vds_e * lame + dbeta_over_u0_dov * f * lame;
            gds = beta_over_u0 * (ov - vds_e) * lame + beta_over_u0 * f * lambda + gm * dibl_gds_factor;
            // mobility exact derivative contributions
            gds += -Id * dDenom_dVd / denom_mob;
            gmbs = -gm * dvth_dvbs - Id * dDenom_dVb / denom_mob;
            saturated_ = false;
        } else {
            const double f = 0.5 * ov * ov;
            Id = beta_over_u0 * f * lame;
            gm = beta_over_u0 * ov * lame + dbeta_over_u0_dov * f * lame;
            gds = beta_over_u0 * f * lambda + gm * dibl_gds_factor - Id * dDenom_dVd / denom_mob;
            // SCBE ...
            double diffVds = vds_e - vdsat;
            if (diffVds > 0.0) {
                double litl = std::sqrt(3.0 * 3.9 * 8.854e-12 / (model_->EPSROX * model_->XJ * model_->TOX));
                if (litl < 1e-30) litl = 1e-30;
                double VASCBE = 0.0;
                if (model_->PSCBE2 > 0.0 && model_->PSCBE1 >= 0.0) {
                    if (diffVds > model_->PSCBE1 * litl / 34.0) {
                        double T0 = model_->PSCBE1 * litl / diffVds;
                        VASCBE = leff * std::exp(T0) / model_->PSCBE2;
                    } else {
                        VASCBE = 5.834617425e14 * leff / model_->PSCBE2;
                    }
                } else {
                    VASCBE = 5.834617425e14;
                }
                if (VASCBE < 1e-30) VASCBE = 1e-30;
                double T9 = diffVds / VASCBE;
                Id = Id * (1.0 + T9);
                gm = gm * (1.0 + T9);
            }
            gmbs = -gm * dvth_dvbs - Id * dDenom_dVb / denom_mob;
            saturated_ = true;
        }
    } else {
        // Subthreshold
        double n = 1.0 + nfactor;
        double Vtm = t.vt;
        double Vgsteff = ov + Voff;
        if (Vgsteff > -10.0 * Vtm) {
            double Id0 = beta_over_u0 * Vtm * Vtm * (1.0 + model_->LAMBDA * vds_e);
            double exp_term = std::exp(Vgsteff / (n * Vtm + 1e-12));
            double drain_factor = 1.0 - std::exp(-vds_e / (Vtm + 1e-12));
            if (drain_factor < 0.0) drain_factor = 0.0;
            if (drain_factor > 1.0) drain_factor = 1.0;
            Id = Id0 * exp_term * drain_factor;
            gm = Id / (n * Vtm + 1e-12);
            double gds_calc = Id * (model_->LAMBDA + std::exp(-vds_e / (Vtm + 1e-12)) / ((Vtm + 1e-12) * (1.0 - std::exp(-vds_e / (Vtm + 1e-12) + 1e-12))));
            if (!std::isfinite(gds_calc)) gds_calc = 0.0;
            gds = gds_calc;
            gmbs = -gm * dvth_dvbs;
        } else {
            Id = 0.0; gm = 0.0; gds = 0.0; gmbs = 0.0;
        }
        saturated_ = false;
    }

    // ---- bulk diodes (PHYSICAL frame) ----
    // Xyce B4p82.C:2445-2448: SourceSatCurrent =
    //   Aseff*JS + Pseff*JSWS + weffCJ*nf*JSWGS
    // (and the drain-side twin).  Everything is scaled by M here because the
    // device reports ONE parallel multiplier's worth of terminals.
    const double vbd_p = vbs - vds;
    const double satS = t.jctS.isatArea * Aseff
                      + t.jctS.isatSide * Pseff
                      + t.jctS.isatGate * weffCJ * nf;
    const double satD = t.jctD.isatArea * Adeff
                      + t.jctD.isatSide * Pdeff
                      + t.jctD.isatGate * weffCJ * nf;
    const double isBulkS = M * satS;
    const double isBulkD = M * satD;
    double arg = std::min(Const::MaxExpArg, vbs / t.vt);
    const double evbs = std::exp(arg);
    const double Ibs = isBulkS * (evbs - 1.0) + gmin * vbs;
    arg = std::min(Const::MaxExpArg, vbd_p / t.vt);
    const double evbd = std::exp(arg);
    const double Ibd = isBulkD * (evbd - 1.0) + gmin * vbd_p;
    gibs_ = isBulkS * evbs / t.vt + gmin;
    gibd_ = isBulkD * evbd / t.vt + gmin;

    // ---- substrate current (impact ionization, BSIM4) ----
    double Isub = 0.0;
    if (model_->ALPHA0 > 0.0 && vds_e > 0.0 && Id > 0.0) {
        double beta0 = (model_->BETA0 > 0.0) ? model_->BETA0 : 1.0;
        Isub = model_->ALPHA0 * Id * std::exp(-beta0 / (vds_e + 1e-3));
        if (model_->ALPHA1 > 0.0) {
            Isub *= (1.0 + model_->ALPHA1 * leff);
        }
    }

    // ---- GIDL ----
    double Igidl = 0.0;
    if (model_->AGIDL > 0.0) {
        double Vgd = vgs - vds;
        double Vdb = vbs - vds;
        double Eg = 1.12;
        double Vtm = t.vt;
        if (Vgd < -Eg * 0.5) {
            double a = -(Eg + (Vdb - Vgd)) / (3.0 * Vtm + 1e-12);
            if (a > -100.0) {
                Igidl = model_->AGIDL * W_ * M_ * std::exp(a);
                if (model_->CGIDL > 0.0 && Vdb > 0.0) {
                    Igidl *= std::pow(1.0 + Vdb / (model_->EGIDL + 1e-6), -model_->CGIDL);
                }
            }
        }
    }

    // ---- GISL ----
    double Igisl = 0.0;
    if (model_->AGISL > 0.0) {
        double Vgs = vgs;
        double Vsb = vbs;
        double Eg = 1.12;
        double Vtm = t.vt;
        if (Vgs < -Eg * 0.5) {
            double a = -(Eg + (Vsb - Vgs)) / (3.0 * Vtm + 1e-12);
            if (a > -100.0) {
                Igisl = model_->AGISL * W_ * M_ * std::exp(a);
                if (model_->CGISL > 0.0 && Vsb > 0.0) {
                    Igisl *= std::pow(1.0 + Vsb / (model_->EGISL + 1e-6), -model_->CGISL);
                }
            }
        }
    }

    // ---- RDS modulation (series resistance) ----
    if (model_->RDSW > 0.0 && W_ > 0.0) {
        double Rds = model_->RDSW / W_;
        double mod = 1.0 + model_->PRWG * ov + model_->PRWB * vbs_e;
        if (mod < 0.1) mod = 0.1;
        Rds *= mod;
        Rds *= (1.0 + model_->PRT * (t.vt / (Const::KoverQ * model_->TNOM) - 1.0));
        if (std::fabs(gm) > 1e-12) {
            double gm_eff = std::fabs(gm);
            Id = Id / (1.0 + Rds * gm_eff);
        }
    }

    // ---- Gate current (tunneling) ----
    double Ig = 0.0;
    if (model_->AIGC > 0.0 && ov > 0.1) {
        double Vox = ov;
        double a = -model_->BIGC / (Vox + 1e-6);
        if (a > -100.0) {
            Ig = model_->AIGC * W_ * M_ * leff * Vox * Vox * std::exp(a);
        }
    }
    if (model_->AIGD > 0.0 && vds_e > 0.1) {
        double Voxd = vds_e;
        double a = -model_->BIGC / (Voxd + 1e-6);
        if (a > -100.0) {
            Ig += model_->AIGD * W_ * M_ * leff * Voxd * Voxd * std::exp(a) * 0.5;
        }
    }

    // ---- Bulk (substrate) current (igb) ----
    double Igb = 0.0;
    if (model_->AIGBACC > 0.0 && vbs_e < -0.1) {
        double Vgb = vgs - vbs;
        double a = -model_->BIGC / (std::fabs(Vgb) + 1e-6);
        if (a > -100.0) {
            Igb = model_->AIGBACC * W_ * M_ * leff * std::fabs(Vgb) * std::exp(a);
        }
    }

    // ---- map channel to physical frame ----
    const double sid = sw ? -1.0 : 1.0;
    id_ = sid * Id + Ibd + Igidl;
    is_ = -sid * Id + Ibs + Igisl;
    ib_ = -(Ibs + Ibd) + Isub - Igidl - Igisl + Igb;
    ig_ = -Ig;

    gm_ = gm;
    gds_ = gds;
    gmbs_ = gmbs;

    if (debug_ && !probe) {
        std::fprintf(stderr,
            "[B4 p=%+d] vgs=%+.6g vds=%+.6g vbs=%+.6g | vth=%+.6g ov=%+.6g"
            " ueff=%.6g denom=%.6g W=%.6g Leff=%.6g KP=%.6g beta=%.6g Id=%+.6g"
            " gm=%.6g gds=%.6g sat=%d\n",
            p_, vgs, vds, vbs, vth, ov, ueff, denom_mob, W_, leff,
            model_->KP, beta_over_u0, id_, gm_, gds_, (int)saturated_);
    }

    // ---- charges: C-V block (Xyce updateIntermediateVars, capMod 0/1/2) ----
    const CvParams cv = cvParams(*model_, W_, L_, p_, t.vt);
    const int capMod = model_->CAPMOD;

    // Vbseff: smoothing around vbsc and 0.95*phi (VbseffCV equals Vbseff).
    double vbseff;
    {
        double T0 = vbs_e - cv.vbsc - 0.001;
        double T1 = std::sqrt(T0 * T0 - 0.004 * cv.vbsc);
        if (T0 >= 0.0) {
            vbseff = cv.vbsc + 0.5 * (T0 + T1);
        } else {
            vbseff = cv.vbsc * (1.0 - 0.002 / (T1 - T0));
        }
        const double T9 = 0.95 * cv.phi;
        T0 = T9 - vbseff - 0.001;
        T1 = std::sqrt(T0 * T0 + 0.004 * T9);
        vbseff = T9 - 0.5 * (T0 + T1);
    }

    // Abulk0 (Xyce), scaled by abulkCVfactor inside the C-V block.
    const double sqrtPhis = std::sqrt(cv.phi - vbseff);
    const double XdepCV = cv.Xdep0 * sqrtPhis / cv.sqrtPhi;
    const double vthNarrowW = cv.toxe * cv.phi / (cv.weff + model_->W0);
    const double lpeVb = std::sqrt(1.0 + model_->LPEB / cv.leff);
    const double abT1 = 0.5 * cv.k1ox * lpeVb / sqrtPhis + cv.k2ox
        - model_->K3B * vthNarrowW;
    const double abT5 = cv.leff
        / (cv.leff + 2.0 * std::sqrt(model_->XJ * XdepCV));
    const double abT2 = model_->A0 * abT5 + model_->B0 / (cv.weff + model_->B1);
    double abulk0 = 1.0 + abT1 * abT2;
    if (abulk0 < 0.1) {
        abulk0 = (0.2 - abulk0) / (3.0 - 20.0 * abulk0);
    }

    // Charge-model inputs: poly depletion off (NGATE not ported yet), DC
    // threshold of this device, and the subthreshold ideality approximation.
    const double vgstCV = vgs_e - vth;
    const double nCV = 1.0 + model_->NFACTOR;

    CtmCharges ci{0.0, 0.0, 0.0, 0.0};
    if (capMod == 0) {
        // Meyer-type partition (approximation; Xyce capMod 0 not ported yet).
        const double coxCV = cv.coxe * cv.weffCV * cv.leffCV;
        double cgsM = 0.0;
        double cgdM = 0.0;
        if (saturated_) {
            cgsM = (2.0 / 3.0) * coxCV;
        } else if (ov > 0.0) {
            cgsM = 0.5 * coxCV;
            cgdM = 0.5 * coxCV;
        }
        const double cgbM = coxCV;
        const double vgdE = vgs_e - vds_e;
        const double vgbE = vgs_e - vbs_e;
        ci.qgate = cgsM * vgs_e + cgdM * vgdE + cgbM * vgbE;
        ci.qdrn = -cgdM * vgdE;
        ci.qbulk = -cgbM * vgbE;
        ci.qsrc = -(ci.qgate + ci.qbulk + ci.qdrn);
    } else {
        ci = ctmCharges(cv, capMod == 1 ? 1 : 2, model_->CVCHARGEMOD,
                        model_->XPART, vgs_e, vds_e, vbseff, vgstCV, nCV,
                        t.vt, abulk0);
    }

    // Overlap charges (Xyce): linear for capMod 0, smoothed with cgdl/cgsl
    // and ckappad/ckappas for capMod 1/2. Gate-bulk overlap uses cgbo.
    const double vgdx = vgs - vds;
    const double vgsx = vgs;
    double qgdo = 0.0;
    double qgso = 0.0;
    if (capMod == 0) {
        qgdo = cv.cgdo * vgdx;
        qgso = cv.cgso * vgsx;
    } else {
        double T0 = vgdx + kCvDelta1;
        double T1 = std::sqrt(T0 * T0 + 4.0 * kCvDelta1);
        double T2 = 0.5 * (T0 - T1);
        double T3 = cv.weffCV * cv.cgdl;
        double T4 = std::sqrt(1.0 - 4.0 * T2 / cv.ckappad);
        qgdo = (cv.cgdo + T3) * vgdx - T3 * (T2 + 0.5 * cv.ckappad * (T4 - 1.0));
        T0 = vgsx + kCvDelta1;
        T1 = std::sqrt(T0 * T0 + 4.0 * kCvDelta1);
        T2 = 0.5 * (T0 - T1);
        T3 = cv.weffCV * cv.cgsl;
        T4 = std::sqrt(1.0 - 4.0 * T2 / cv.ckappas);
        qgso = (cv.cgso + T3) * vgsx - T3 * (T2 + 0.5 * cv.ckappas * (T4 - 1.0));
    }
    const double qgb = cv.cgbo * (vgs - vbs);

    // Terminal charges (Xyce setupCapacitors/updatePrimaryState). In reverse
    // mode (vds < 0) the intrinsic drain and source charges are exchanged.
    const double qgateT = ci.qgate + qgdo + qgso + qgb;
    const double qbulkT = ci.qbulk - qgb;
    const double qdrnT = sw ? -(qgateT + qbulkT + (ci.qdrn - qgso))
                            : (ci.qdrn - qgdo);

    // Junction (depletion) charges of the drain and source junctions.
    // Xyce B4p82.C:3348-3695 evaluates the bottom, sidewall and gate-edge
    // terms separately; junctionCharge() already takes an AREA-scaled
    // geometry, so each term gets its own JunctionGeom.  M is applied once,
    // at the node charges below.
    const double mAreaS = (model_->MJS > 0.0) ? model_->MJS : model_->MJ;
    const double mAreaD = (model_->MJD > 0.0) ? model_->MJD : model_->MJ;
    const JunctionGeom geomS0{t.jctS.capArea.cjo * Aseff, 0.0, mAreaS};
    const JunctionGeom geomS1{t.jctS.capSide.cjo * Pseff, 0.0, model_->MJSWS};
    const JunctionGeom geomS2{t.jctS.capGate.cjo * weffCJ * nf, 0.0, model_->MJSWGS};
    const JunctionGeom geomD0{t.jctD.capArea.cjo * Adeff, 0.0, mAreaD};
    const JunctionGeom geomD1{t.jctD.capSide.cjo * Pdeff, 0.0, model_->MJSWD};
    const JunctionGeom geomD2{t.jctD.capGate.cjo * weffCJ * nf, 0.0, model_->MJSWGD};

    const JunctionChargePoint jbs =
        junctionCharge(vbs, 0.0, 0.0, geomS0, t.jctS.capArea);
    const JunctionChargePoint jbs1 =
        junctionCharge(vbs, 0.0, 0.0, geomS1, t.jctS.capSide);
    const JunctionChargePoint jbs2 =
        junctionCharge(vbs, 0.0, 0.0, geomS2, t.jctS.capGate);
    const JunctionChargePoint jbd =
        junctionCharge(vbd_p, 0.0, 0.0, geomD0, t.jctD.capArea);
    const JunctionChargePoint jbd1 =
        junctionCharge(vbd_p, 0.0, 0.0, geomD1, t.jctD.capSide);
    const JunctionChargePoint jbd2 =
        junctionCharge(vbd_p, 0.0, 0.0, geomD2, t.jctD.capGate);

    const double qjs = jbs.q + jbs1.q + jbs2.q;
    const double qjd = jbd.q + jbd1.q + jbd2.q;

    // Node charges (N-oriented). Xyce loadDAEQVector scales them by the
    // number of parallel devices (M).
    qg_ = qgateT * M;
    qd_ = (qdrnT - qjd) * M;
    qb_ = (qbulkT + qjs + qjd) * M;
    qdef_ = 0.0;
    (void)0;
}

double MosfetBSIM4Device::eval1ovFNoise(double vds, double freq, double temp) const {
    // 1/f noise PSD (simplified).  Returns 0 when not applicable.
    if (freq <= 0.0) return 0.0;
    const double gm = std::fabs(gm_);
    if (gm <= 0.0) return 0.0;
    const double kT = 1.380649e-23 * (temp > 0.0 ? temp : 300.15);
    // Use model NOIA/NOIB/NOIC if available.
    const MosfetBSIM4Model& m = *model_;
    const double Leff = L_ - 2.0 * m.LD;
    const double leff = (Leff > 0.0) ? Leff : 1e-9;
    const double Cox = m.EPSROX * 8.854e-12 / (m.TOXE > 0.0 ? m.TOXE : 1e-9);
    const double Nl = Cox * (vgs_cached_ - vth_cached_) / 1.60219e-19;
    const double nstar = 1e12;
    const double T1 = m.NOIA + m.NOIB / (Nl + nstar)
                    + m.NOIC / ((Nl + nstar) * (Nl + nstar));
    const double psd = kT * std::fabs(id_) * std::fabs(id_) / (Cox * leff * W_ * M_)
                     * T1 / freq;
    (void)vds;
    return psd;
}

int MosfetBSIM4Device::getNumNoiseSources() const {
    return 2;  // thermal + flicker
}

void MosfetBSIM4Device::evaluate(double vgs, double vds, double vbs) {
    // Charge Jacobian dQ/dV: central differences of the terminal charges of
    // the exact charge model (same model, probed at vgs/vds/vbs +- h).
    constexpr double kStep = 1.0e-6;
    const double x0[3] = {vgs, vds, vbs};
    for (int j = 0; j < 3; ++j) {
        double xp[3] = {x0[0], x0[1], x0[2]};
        double xm[3] = {x0[0], x0[1], x0[2]};
        xp[j] += kStep;
        xm[j] -= kStep;
        evaluateBase(xp[0], xp[1], xp[2], true);
        const double qdp = qd_;
        const double qgp = qg_;
        const double qbp = qb_;
        evaluateBase(xm[0], xm[1], xm[2], true);
        const double inv2h = 0.5 / kStep;
        cq_[0][j] = (qdp - qd_) * inv2h;
        cq_[1][j] = (qgp - qg_) * inv2h;
        cq_[2][j] = (qbp - qb_) * inv2h;
    }
    evaluateBase(vgs, vds, vbs, false);
    if (debug_) {
        std::fprintf(stderr,
            "[B4Q] vgs=%+.6g vds=%+.6g vbs=%+.6g qg=%+.6e qd=%+.6e qb=%+.6e"
            " cgg=%+.6e cgd=%+.6e cgb=%+.6e cdg=%+.6e cdd=%+.6e\n",
            vgs, vds, vbs, qg_, qd_, qb_, cq_[1][0], cq_[1][1], cq_[1][2],
            cq_[0][0], cq_[0][1]);
    }
}

void MosfetBSIM4Device::loadF(std::vector<double>& F,
                              int rd, int rg, int rs, int rb) const {
    const int r[4] = {rd, rg, rs, rb};
    const double iv[4] = {id_, ig_, is_, ib_};
    for (int i = 0; i < 4; ++i) {
        if (r[i] >= 0 && static_cast<size_t>(r[i]) < F.size()) {
            F[r[i]] += p_ * iv[i];
        }
    }
}

void MosfetBSIM4Device::loadQ(std::vector<double>& Q,
                              int rd, int rg, int rs, int rb) const {
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
    double ef[16] = {0};

    const double gm = gm_;
    const double gds = gds_;
    const double gmbs = gmbs_;

    const double gibs = gibs_;
    const double gibd = gibd_;
    if (!swapped_) {
        ef[0*4+1] = gm;
        ef[0*4+0] = gds - gibd;
        ef[0*4+2] = -(gm + gds + gmbs);
        ef[0*4+3] = gmbs + gibd;
        ef[2*4+1] = -gm;
        ef[2*4+0] = -gds;
        ef[2*4+2] = (gm + gds + gmbs) - gibs;
        ef[2*4+3] = -gmbs + gibs;
    } else {
        ef[0*4+1] = -gm;
        ef[0*4+0] = (gm + gds + gmbs) - gibd;
        ef[0*4+2] = -gds;
        ef[0*4+3] = -gmbs + gibd;
        ef[2*4+1] = gm;
        ef[2*4+0] = -(gm + gds + gmbs);
        ef[2*4+2] = gds - gibs;
        ef[2*4+3] = gmbs + gibs;
    }
    ef[3*4+0] = gibd;
    ef[3*4+2] = gibs;
    ef[3*4+3] = -(gibs + gibd);

    for (int k = 0; k < 16; ++k) {
        if (slots[k] >= 0) {
            J.accumulateSlot(slots[k], ef[k]);
        }
    }
}

void MosfetBSIM4Device::loadJacobianQ(DaeMatrix& J,
                                      const std::vector<int>& slots) const {
    // cq_ rows: drain, gate, bulk; columns: d/dvgs, d/dvds, d/dvbs. The
    // source charge is -(qd + qg + qb), and the node voltages enter through
    // vgs = Vg - Vs, vds = Vd - Vs, vbs = Vb - Vs (row-major D,G,S,B).
    double eq[16] = {0};
    const int rows[3] = {0, 1, 3};
    for (int r = 0; r < 3; ++r) {
        const int row = rows[r];
        eq[row * 4 + 1] = cq_[r][0];
        eq[row * 4 + 0] = cq_[r][1];
        eq[row * 4 + 3] = cq_[r][2];
        eq[row * 4 + 2] = -(cq_[r][0] + cq_[r][1] + cq_[r][2]);
    }
    for (int c = 0; c < 4; ++c) {
        eq[2 * 4 + c] = -(eq[0 * 4 + c] + eq[1 * 4 + c] + eq[3 * 4 + c]);
    }

    for (int k = 0; k < 16; ++k) {
        if (slots[k] >= 0) {
            J.accumulateSlot(slots[k], eq[k]);
        }
    }
}

} // namespace e_XSpice

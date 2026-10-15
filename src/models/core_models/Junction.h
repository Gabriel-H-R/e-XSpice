/**
 * @author Gabriel Henrique Silva
 * @file Junction.h
 * @brief Shared PN-junction physics used by EVERY model with junctions
 *        (diode today, BJT base-emitter/base-collector tomorrow, JFET/MOS
 *        junctions later).
 *
 * Two quantities live here, both header-only INLINED on purpose — they are
 * evaluated once per junction per Newton iteration, and the compiler folds
 * them into the device's loadF/loadQ with zero call overhead:
 *
 *  1. junctionCharge():  Q(V) and C(V) = dQ/dV — the SPICE3f5/Xyce
 *     construction where the depletion charge and the capacitance are BOTH
 *     continuous (C0) at V = FC·VJ:
 *
 *        V < FC·VJ:  Q = TT·I + VJ·C0·(1 − u·sarg)/(1 − M),   C = TT·G + C0·sarg
 *        V ≥ FC·VJ:  Q = TT·I + C0·F1 + (C0/F2)·[F3·(V − FC·VJ)
 *                                          + ½·(M/VJ)·(V² − (FC·VJ)²)]
 *                    C = TT·G + (C0/F2)·[F3 + (M/VJ)·V]
 *        with u = 1 − V/VJ, sarg = u^(−M), F2 = (1−FC)^(1+M),
 *        F3 = 1−FC·(1+M), F1 = VJ·(1 − (1−FC)^(1−M))/(1−M).
 *
 *     The TT·I term (transit charge, SPICE's CJ-limited minority-carrier
 *     storage) is applied whenever TT ≠ 0, INDEPENDENTLY of CJO — SPICE
 *     semantics: a junction with TT>0 and CJO=0 still stores charge.
 *
 *  2. scaleJunctionCap(): VJ(T), CJO(T) and the F2/F3/F1 constants — the
 *     bandgap-based temperature scaling, computed ONCE per card (setup),
 *     never inside an iteration.
 *
 * Everything is plain data (no virtual calls, no maps) so the hot paths
 * read it through a pointer.
 */

#ifndef JUNCTION_H
#define JUNCTION_H

#include "../core/DeviceConstants.h"
#include <cmath>

namespace e_XSpice {

/**
 * @brief Per-instance geometry of one junction (area ALREADY applied).
 */
struct JunctionGeom {
    double cjo = 0.0;   // CJO(T)·AREA
    double tt  = 0.0;   // transit time coefficient (SPICE TT); 0 disables
    double m   = 0.33;  // grading coefficient
};

/**
 * @brief Temperature-scaled capacitance constants of one junction.
 *        cjo is per unit area (the instance multiplies by AREA).
 */
struct JunctionCap {
    double vj   = 0.0;  // VJ(T)      — junction potential
    double cjo  = 0.0;  // CJO(T)      — zero-bias capacitance (per area)
    double tdep = 0.0;  // FC·VJ(T)   — transition voltage
    double f2   = 0.0;  // (1−FC)^(1+M)
    double f3   = 0.0;  // 1−FC·(1+M)
    double f1   = 0.0;  // depletion charge at the transition point
};

/**
 * @brief Evaluated junction charge/capacitance at a given voltage.
 */
struct JunctionChargePoint {
    double q = 0.0;   // charge [C]
    double c = 0.0;   // capacitance dQ/dV [F]
};

/**
 * @brief Q(V) and C(V) = dQ/dV of a junction (C0-continuous at FC·VJ).
 *
 * @param vd   junction voltage (forward positive)
 * @param id   junction current at vd (for the TT transit term)
 * @param gd   dI/dV at vd (matching derivative)
 * @param geom instance geometry (AREA-scaled)
 * @param cap  card constants at the evaluation temperature
 */
inline JunctionChargePoint junctionCharge(double vd, double id, double gd,
                                         const JunctionGeom& geom,
                                         const JunctionCap& cap) {
    double q = 0.0;
    double c = 0.0;

    // Transit (minority-carrier) charge: TT·I. Independent of CJO — this is
    // what limits f_T in a device with negligible depletion capacitance.
    if (geom.tt != 0.0) {
        q += geom.tt * id;
        c += geom.tt * gd;
    }

    if (geom.cjo != 0.0 && cap.vj > 0.0) {
        if (vd < cap.tdep) {
            const double u = 1.0 - vd / cap.vj;
            const double sarg = std::exp(-geom.m * std::log(u));   // u^(−M)
            q += cap.vj * geom.cjo * (1.0 - u * sarg) / (1.0 - geom.m);
            c += geom.cjo * sarg;
        } else {
            const double czof2 = geom.cjo / cap.f2;
            const double motjctpot = geom.m / cap.vj;
            q += geom.cjo * cap.f1 +
                 czof2 * (cap.f3 * (vd - cap.tdep) +
                          0.5 * motjctpot *
                              (vd * vd - cap.tdep * cap.tdep));
            c += czof2 * (cap.f3 + motjctpot * vd);
        }
    }

    return {q, c};
}

/**
 * @brief Temperature scaling of the junction potential/capacitance plus the
 *        continuity constants F2/F3/F1 (computed once per card).
 *
 * Formula (SPICE3 / Xyce N_DEV_Diode+):
 *   VJ(T)  = VJ + 2·vt·[1.5·ln(T/Tref) + Q·(Eg(0)·...)]  (bandgap shift)
 *   CJO(T) = CJO·(1 + M·(4e-4·(T−Tref) − gmanew)) /
 *                (1 + M·(4e-4·(TNOM−Tref) − gmaold))
 */
inline JunctionCap scaleJunctionCap(double vj, double m, double fc,
                                    double cjo, double tempKelvin,
                                    double tnomKelvin) {
    using namespace Const;

    JunctionCap cap;
    const double vt    = KoverQ * tempKelvin;
    const double vtnom = KoverQ * tnomKelvin;

    const double egfet  = Eg0 - (AlphaEg * tempKelvin * tempKelvin) /
                                   (tempKelvin + BetaEg);
    const double egfet1 = Eg0 - (AlphaEg * tnomKelvin * tnomKelvin) /
                                   (tnomKelvin + BetaEg);
    const double arg  = -egfet  / (2.0 * Boltz * tempKelvin) +
                         Eg300 / (2.0 * Boltz * RefTemp);
    const double arg1 = -egfet1 / (2.0 * Boltz * tnomKelvin) +
                        Eg300 / (2.0 * Boltz * RefTemp);
    const double fact2 = tempKelvin / RefTemp;
    const double fact1 = tnomKelvin / RefTemp;
    const double pbfact  = -2.0 * vt    * (1.5 * std::log(fact2) + Q * arg);
    const double pbfact1 = -2.0 * vtnom * (1.5 * std::log(fact1) + Q * arg1);

    const double pbo = (vj - pbfact1) / fact1;
    const double gmaold = (vj - pbo) / pbo;
    cap.vj   = pbfact + fact2 * pbo;
    const double gmanew = (cap.vj - pbo) / pbo;

    cap.cjo = cjo / (1.0 + m * (4.0e-4 * (tnomKelvin - RefTemp) - gmaold));
    cap.cjo *= 1.0 + m * (4.0e-4 * (tempKelvin - RefTemp) - gmanew);

    cap.tdep = fc * cap.vj;
    const double xfc = std::log(1.0 - fc);
    cap.f2 = std::exp((1.0 + m) * xfc);        // (1−FC)^(1+M)
    cap.f3 = 1.0 - fc * (1.0 + m);
    cap.f1 = cap.vj * (1.0 - std::exp((1.0 - m) * xfc)) / (1.0 - m);

    return cap;
}

/**
 * @brief SPICE saturation-current temperature scaling
 *        IS(T) = IS·exp((T/TNOM − 1)·EG/(N·vt) + (XTI/N)·ln(T/TNOM)).
 */
inline double scaleSatCurrent(double is, double n, double eg, double xti,
                              double tempKelvin, double tnomKelvin) {
    const double vt = Const::KoverQ * tempKelvin;
    return is * std::exp(((tempKelvin / tnomKelvin) - 1.0) * eg / (n * vt) +
                         (xti / n) * std::log(tempKelvin / tnomKelvin));
}

} // namespace e_XSpice

#endif // JUNCTION_H

/**
 * @author Gabriel Henrique Silva
 * @file BjtModel.h
 * @brief Shared BJT model card — level 1 (Ebers-Moll transport) and level 2
 *        (Gummel-Poon) for NPN and PNP.
 *
 * ================================================================
 * LEVEL 1 — transport form of Ebers-Moll (SPICE3 / Xyce level 1)
 * ================================================================
 * With  vbe = v_b − v_e,  vbc = v_b − v_c  (PNP: see the polarity note in
 * Bjt.h — the device remaps the voltages, the card is polarity-agnostic),
 *
 *   vte = NF·Vt,  vtc = NR·Vt
 *
 *   I_BE = β_F·(exp(vbe/vte) − 1) + gmin·vbe      (nonideal base current)
 *   I_BC = β_R·(exp(vbc/vtc) − 1) + gmin·vbc
 *   I_CT = IS·(exp(vbe/vte) − exp(vbc/vtc))·K_E   (transport current)
 *
 *   β_F = IS/BF,  β_R = IS/BR
 *   K_E = 1 − vbc/VAF − vbe/VAR                    (base-width modulation)
 *
 * ================================================================
 * LEVEL 2 — Gummel-Poon (SPICE3 Gummel-Poon, Xyce N_DEV_BJT level 2)
 * ================================================================
 * Selected with LEVEL=2 on the card. Adds, on top of level 1:
 *
 *   Base charge (high injection + Early, one normalized charge):
 *     q1 = 1 / (1 − vbc/VAF − vbe/VAR)
 *     q2 = IS·(exp(vbe/NF·Vt) − 1)/IKF + IS·(exp(vbc/NR·Vt) − 1)/IKR
 *     qb = q1·(1 + sqrt(1 + 4·q2)) / 2
 *   The transport current is divided by qb:
 *     I_CT = IS·(exp(vbe/vte) − exp(vbc/vtc)) / qb
 *   so that β falls at high current (β_eff ≈ BF/qb) while the low-current
 *   β stays exactly BF — level 1 is the qb → 1 limit of level 2.
 *
 *   Non-ideal base currents (low-current slope degradation):
 *     I_BE += ISE·(exp(vbe/NE·Vt) − 1),  I_BC += ISC·(exp(vbc/NC·Vt) − 1)
 *
 *   Ohmic resistances (internal nodes, see Bjt.h):
 *     RB between external base and the intrinsic base,
 *     RC between external collector and the intrinsic collector,
 *     RE between external emitter and the intrinsic emitter.
 *   With RB=RC=RE=0 the intrinsic nodes alias the externals and level 2
 *   costs exactly the qb evaluation — no extra unknowns, no extra stamp.
 *
 *   Substrate (optional 4th terminal Q C B E S model):
 *     collector-substrate diode ISS·(exp(vsc/vt) − 1) + depletion charge
 *     CJS/VJS/MJS. Without the 4th terminal the substrate current/charge
 *     are identically zero and no state is allocated.
 *
 * Temperature (computed ONCE per card, never per instance):
 *   IS(T) per junction via the SPICE formula (see Junction.h),
 *   ISE(T)/ISC(T) with the NE/NC emission coefficients,
 *   VJ(T)/CJO(T) per junction with the bandgap shift, vcrit for PNJLIM.
 *   IKF/IKR and RB/RC/RE are temperature-independent (documented
 *   simplification: SPICE scales IKF with a power of T/TNOM whose effect
 *   is second-order next to the exponential IS(T)).
 *
 * Deliberately NOT in this card (queued as future work, documented in the
 * roadmap): base-resistance modulation (RBM/IRB), excess phase (PTF is
 * parsed and stored but only affects AC), avalanche (BV/IBV), tristate
 * substrate PNP handling.
 */

#ifndef BJT_MODEL_H
#define BJT_MODEL_H

#include "../core/DeviceModel.h"
#include "../core/ModelCard.h"
#include "Junction.h"
#include <string>

namespace e_XSpice {

/**
 * @brief Temperature-scaled derived quantities of the card (plain data).
 */
struct BjtThermal {
    double vt      = 0.0;   // kT/q
    double vte     = 0.0;   // NF·vt   (BE emission)
    double vtc     = 0.0;   // NR·vt   (BC emission)
    double vteE    = 0.0;   // NE·vt   (BE non-ideal, level 2)
    double vtcC    = 0.0;   // NC·vt   (BC non-ideal, level 2)
    double beIsat  = 0.0;   // IS(T) for the base-emitter junction
    double bcIsat  = 0.0;   // IS(T) for the base-collector junction
    double betaf   = 0.0;   // IS(T)/BF
    double betar   = 0.0;   // IS(T)/BR
    double iseT    = 0.0;   // ISE(T)  (level 2)
    double iscT    = 0.0;   // ISC(T)  (level 2)
    double issT    = 0.0;   // ISS(T)  (substrate, level 2)
    double ikf     = 0.0;   // forward high-injection knee (level 2)
    double ikr     = 0.0;   // reverse high-injection knee (level 2)
    double beVcrit = 0.0;   // PNJLIM critical voltage, BE junction
    double bcVcrit = 0.0;   // PNJLIM critical voltage, BC junction
    double tf      = 0.0;   // forward transit time (per unit area)
    double tr      = 0.0;   // reverse transit time (per unit area)
    double gmin    = 1e-12; // junction conductance in parallel
    JunctionCap beCap{};    // BE junction capacitance constants
    JunctionCap bcCap{};    // BC junction capacitance constants
    JunctionCap csCap{};    // collector-substrate capacitance constants
};

/**
 * @brief Per-instance overrides (one netlist line).
 *
 * Own option table (same machinery as card parameters): the parser fills
 * it generically. Range validation stays in the Bjt constructor.
 */
struct BjtInstanceParams {
    double AREA = 1.0;       // area multiplier (currents AND capacitances)
    double TEMP = 0.0;       // instance temperature (Kelvin, when TEMPGiven)
    double DTEMP = 0.0;      // delta temperature (K)
    bool   TEMPGiven = false;

    using Table = std::vector<TypedParamDesc<BjtInstanceParams>>;
    static const Table& paramTable() {
        static const Table t = {
            makeParam("AREA", &BjtInstanceParams::AREA,
                      "area multiplier"),
            makeCelsiusParam("TEMP", &BjtInstanceParams::TEMP,
                             "instance temperature",
                             [](BjtInstanceParams& o) { o.TEMPGiven = true; }),
            makeParam("DTEMP", &BjtInstanceParams::DTEMP,
                      "delta temperature (K)"),
        };
        return t;
    }
};

class BjtModel : public DeviceModel {
public:
    /**
     * @param pnp true → PNP card, false → NPN card
     */
    BjtModel(std::string name, bool pnp);

    std::string typeName() const override { return pnp_ ? "PNP" : "NPN"; }

    /**
     * @brief true for a PNP card (the device flips the terminal signs).
     */
    bool isPnp() const { return pnp_; }

    // ---- card parameters (SPICE3 Gummel-Poon defaults) ----
    int    LEVEL = 1;        // 1 = Ebers-Moll transport, 2 = Gummel-Poon
    double IS  = 1e-16;      // transport saturation current
    double BF  = 100.0;      // ideal forward beta
    double BR  = 1.0;        // ideal reverse beta
    double NF  = 1.0;        // forward current emission coefficient
    double NR  = 1.0;        // reverse current emission coefficient
    double VAF = 1e99;       // forward Early voltage
    double VAR = 1e99;       // reverse Early voltage
    double TF  = 0.0;        // forward transit time
    double TR  = 0.0;        // reverse transit time
    // base-emitter junction capacitance
    double CJE = 0.0;        // BE zero-bias depletion capacitance (per area)
    double VJE = 0.75;       // BE built-in potential
    double MJE = 0.33;       // BE grading coefficient
    double FE  = 0.5;        // BE forward-bias capacitance coefficient
    // base-collector junction capacitance
    double CJC = 0.0;        // BC zero-bias depletion capacitance (per area)
    double VJC = 0.75;       // BC built-in potential
    double MJC = 0.33;       // BC grading coefficient
    double FC  = 0.5;        // BC forward-bias capacitance coefficient
    // ---- level-2 (Gummel-Poon) parameters ----
    double IKF = 0.0;        // forward high-injection knee current (0 = off)
    double IKR = 0.0;        // reverse high-injection knee current (0 = off)
    double ISE = 0.0;        // BE non-ideal saturation current (0 = off)
    double NE  = 1.5;        // BE non-ideal emission coefficient
    double ISC = 0.0;        // BC non-ideal saturation current (0 = off)
    double NC  = 2.0;        // BC non-ideal emission coefficient
    double RB  = 0.0;        // zero-bias base resistance (0 = no int. node)
    double RC  = 0.0;        // collector ohmic resistance (0 = no int. node)
    double RE  = 0.0;        // emitter ohmic resistance (0 = no int. node)
    double ISS = 0.0;        // substrate saturation current (0 = off)
    double CJS = 0.0;        // zero-bias C-S depletion capacitance (per area)
    double VJS = 0.75;       // C-S built-in potential
    double MJS = 0.0;        // C-S grading coefficient (0 = SPICE default .33)
    double PTF = 0.0;        // excess phase (parsed; AC-only, see header)
    // temperature
    double EG   = 1.11;      // bandgap voltage
    double XTI  = 3.0;       // saturation-current temperature exponent
    double TNOM = 300.15;    // Kelvin
    double TEMP = 300.15;    // Kelvin
    double GMIN = 1e-12;     // parallel junction conductance

    /**
     * @brief True when any ohmic resistance needs an internal node.
     */
    bool hasRb() const { return RB > 0.0; }
    bool hasRc() const { return RC > 0.0; }
    bool hasRe() const { return RE > 0.0; }

    /**
     * @brief Number of internal unknowns (0..3) for this card.
     */
    int numInternalNodes() const {
        return (hasRb() ? 1 : 0) + (hasRc() ? 1 : 0) + (hasRe() ? 1 : 0);
    }

    /**
     * @brief Validates, applies guards and scales to TEMP.
     */
    bool processParams(std::string& message) override;

    /**
     * @brief Recomputes thermal_ for T (Kelvin).
     */
    void updateTemperature(double tempKelvin) override;

    /**
     * @brief Snapshot for an arbitrary T (instance TEMP overrides).
     */
    BjtThermal computeThermal(double tempKelvin) const;

    /**
     * @brief Active thermal snapshot (devices read it through a pointer).
     */
    const BjtThermal& thermal() const { return thermal_; }

    /**
     * @brief Process-wide default NPN card (SPICE defaults).
     */
    static RCP<const BjtModel> defaultModel();

private:
    bool pnp_;
    BjtThermal thermal_;

    void applyGuards();
    BjtThermal scaleTo(double tempKelvin) const;
};

} // namespace e_XSpice

#endif // BJT_MODEL_H

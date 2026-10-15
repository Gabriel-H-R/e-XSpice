/**
 * @author Gabriel Henrique Silva
 * @file BjtModel.h
 * @brief Shared BJT model card — transport (Ebers-Moll) model, Xyce level 1.
 *
 * ================================================================
 * CURRENT MODEL — transport form of Ebers-Moll (SPICE3 / Xyce level 1)
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
 * Terminal currents (positive INTO the transistor, NPN):
 *   i_b = I_BE + I_BC
 *   i_c = I_CT − I_BC
 *   i_e = −i_b − i_c
 *
 * Sanity checks of the formulation:
 *   active (vbc<0, vbe>0):  i_c ≈ IS·exp(vbe/vte), i_b ≈ β_F·exp(vbe/vte)
 *                           ⇒ β = i_c/i_b = BF  exactly ✓
 *   reverse (vbc>0, vbe<0): β = BR  exactly ✓
 *   saturation, vbe→vbc:    i_c → 0 and reverses ✓
 *   Early: Vce↑ ⇒ vbc↓ ⇒ K_E↑ ⇒ i_c↑, g_o = ∂i_c/∂v_c > 0 ✓
 *
 * ================================================================
 * CHARGE MODEL
 * ================================================================
 *   Q_BE = Q_dep(BE, vbe) + τF·I_CT        Q_BC = Q_dep(BC, vbc) − τR·I_CT
 * with Q_dep the shared SPICE3f5 junction charge (Junction.h). The transit
 * charge is attributed to the junction it is transported FROM (forward on
 * BE, reverse on BC — SPICE level-1 partition), so dQ/dx carries the exact
 * cross derivatives (∂Q_BE/∂vbc = τF·∂I_CT/∂vbc).
 *
 * ================================================================
 * TEMPERATURE (computed ONCE per card, never per instance)
 * ================================================================
 *   IS(T) per junction via the SPICE formula (see Junction.h),
 *   VJ(T)/CJO(T) per junction with the bandgap shift, vcrit for PNJLIM.
 * Ohmic resistances (RB/RC/RE), high-injection knees (IKF/IKR), the
 * substrate terminal and XPARTF are deliberately NOT in this card — they
 * need internal nodes and are queued as the next model level (roadmap).
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
    double beIsat  = 0.0;   // IS(T) for the base-emitter junction
    double bcIsat  = 0.0;   // IS(T) for the base-collector junction
    double betaf   = 0.0;   // IS(T)/BF
    double betar   = 0.0;   // IS(T)/BR
    double beVcrit = 0.0;   // PNJLIM critical voltage, BE junction
    double bcVcrit = 0.0;   // PNJLIM critical voltage, BC junction
    double tf      = 0.0;   // forward transit time (per unit area)
    double tr      = 0.0;   // reverse transit time (per unit area)
    double gmin    = 1e-12; // junction conductance in parallel
    JunctionCap beCap{};    // BE junction capacitance constants
    JunctionCap bcCap{};    // BC junction capacitance constants
};

/**
 * @brief Per-instance overrides (one netlist line).
 */
struct BjtInstanceParams {
    double AREA = 1.0;       // area multiplier (currents AND capacitances)
    double TEMP = 0.0;       // instance temperature (Kelvin, when TEMPGiven)
    double DTEMP = 0.0;      // delta temperature (K)
    bool   TEMPGiven = false;
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

    // ---- card parameters (SPICE3 level-1 defaults) ----
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
    // temperature
    double EG   = 1.11;      // bandgap voltage
    double XTI  = 3.0;       // saturation-current temperature exponent
    double TNOM = 300.15;    // Kelvin
    double TEMP = 300.15;    // Kelvin
    double GMIN = 1e-12;     // parallel junction conductance

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
     * @brief Builds a card from netlist key=value pairs (TNOM/TEMP arrive
     *        in CELSIUS, SPICE convention, and are stored Kelvin).
     *        Aliases: VA→VAF, VB→VAR, CJ→CJE, PE→VJE, ME→MJE,
     *                 PC→VJC, MC→MJC, CJS→CJC.
     */
    static RCP<DeviceModel> fromCard(const std::string& name,
                                     const ParamMap& params,
                                     bool pnp);

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

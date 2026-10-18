/**
 * @author Gabriel Henrique Silva
 * @file DiodeModel.h
 * @brief Shared PN-junction model card (Xyce level-1 physics).
 *
 * The card holds everything SHARED by all diodes referencing it: the raw
 * SPICE parameters plus the temperature-scaled derived quantities, computed
 * ONCE in processParams()/updateTemperature() instead of once per instance.
 *
 * Temperature model (Xyce N_DEV_Diode, SPICE3):
 *   Vt = kT/q,  Vte = N·Vt,  VteR = NR·Vt,  VteBRK = NBV·Vt
 *   IS(T)  = IS·exp(((T/TNOM)−1)·EG/(N·Vt) + (XTI/N)·ln(T/TNOM))   (+ ISR twin)
 *   VJ(T), CJO(T): bandgap-based junction potential/capacitance scaling
 *   tVcrit = Vte·ln(Vte/(√2·IS))                                   (PNJLIM)
 *   tBrkdwnV: fixed-point fit so that I_D = −IBV at Vd = −BV
 *   F2 = (1−FC)^(1+M), F3 = 1−FC·(1+M): make Q and C BOTH continuous at
 *   Vd = FC·VJ (SPICE3f5 construction)
 *
 * Per-instance overrides (AREA, multiplicity M, TEMP) live in
 * DiodeInstanceParams and are owned by the Diode, not the card. An instance
 * that sets TEMP gets its own Thermal snapshot (computed once in its
 * constructor); otherwise it reads the model's thermal directly — zero
 * per-iteration cost either way.
 */

#ifndef DIODE_MODEL_H
#define DIODE_MODEL_H

#include "../core/DeviceModel.h"
#include "../core/ModelCard.h"
#include <string>

namespace e_XSpice {

/**
 * @brief Temperature-scaled derived quantities (computed once per model).
 *        Plain data: devices read it through a pointer, no virtual calls,
 *        no maps, no recomputation inside Newton loops.
 */
struct DiodeThermal {
    double vt      = 0.0;   // kT/q
    double vte     = 0.0;   // N·vt
    double vteR    = 0.0;   // NR·vt
    double vteBRK  = 0.0;   // NBV·vt
    double isat    = 0.0;   // IS(T)
    double isatR   = 0.0;   // ISR(T)
    double tJctPot = 0.0;   // VJ(T)
    double tJctCap = 0.0;   // CJO(T)
    double tDepCap = 0.0;   // FC·VJ(T)
    double F2      = 0.0;   // (1−FC)^(1+M)
    double F3      = 0.0;   // 1−FC·(1+M)
    double tF1     = 0.0;   // depletion charge at the transition point
    double tBrkdwnV = 0.0;  // breakdown voltage fitted from BV/IBV
    double tVcrit  = 0.0;   // critical voltage for PNJLIM
    double gspr    = 0.0;   // 1/RS (0 when RS = 0)
    double ikf     = 0.0;   // IKF at temperature
};

/**
 * @brief Per-instance diode parameters (one netlist line).
 *
 * Declares its own option table (same machinery as card parameters):
 * the parser fills it generically, so new instance options never touch
 * parser code. Range validation stays in the Diode constructor.
 */
struct DiodeInstanceParams {
    double AREA = 1.0;   // area multiplier
    double M    = 1.0;   // multiplicity (parallel devices)
    double TEMP = 0.0;   // instance temperature (Kelvin, when TEMPGiven)
    double DTEMP = 0.0;  // delta temperature (K)
    bool   TEMPGiven = false;

    using Table = std::vector<TypedParamDesc<DiodeInstanceParams>>;
    static const Table& paramTable() {
        static const Table t = {
            makeParam("AREA", &DiodeInstanceParams::AREA,
                      "area multiplier"),
            makeParam("M", &DiodeInstanceParams::M,
                      "multiplicity (parallel devices)"),
            makeCelsiusParam("TEMP", &DiodeInstanceParams::TEMP,
                             "instance temperature",
                             [](DiodeInstanceParams& o) { o.TEMPGiven = true; }),
            makeParam("DTEMP", &DiodeInstanceParams::DTEMP,
                      "delta temperature (K)"),
        };
        return t;
    }
};

class DiodeModel : public DeviceModel {
public:
    explicit DiodeModel(std::string name);

    std::string typeName() const override { return "D"; }

    // ---- card parameters (defaults per SPICE/Xyce level 1) ----
    double IS  = 1e-14;
    double N   = 1.0;
    double RS  = 0.0;
    double CJO = 0.0;
    double VJ  = 1.0;
    double M   = 0.5;
    double FC  = 0.5;
    double BV  = 1e99;
    double IBV = 1e-3;
    double NBV = 1.0;
    double TT  = 0.0;
    double EG  = 1.11;
    double XTI = 3.0;
    double IKF = 0.0;
    double ISR = 0.0;
    double NR  = 2.0;
    double TNOM = 300.15;   // Kelvin
    double TEMP = 300.15;   // Kelvin
    double GMIN = 1e-12;
    bool   BVGiven = false;

    /**
     * @brief Validates + guards + scales to TEMP. See DeviceModel.
     */
    bool processParams(std::string& message) override;

    /**
     * @brief Recomputes thermal_ for T (Kelvin).
     */
    void updateTemperature(double tempKelvin) override;

    /**
     * @brief Computes a thermal snapshot for an arbitrary T without
     *        touching the model (instance TEMP overrides).
     */
    DiodeThermal computeThermal(double tempKelvin) const;

    /**
     * @brief Active temperature-scaled quantities (valid after
     *        processParams; devices read this through a pointer).
     */
    const DiodeThermal& thermal() const { return thermal_; }

    /**
     * @brief Process-wide default model (SPICE defaults). Shared,
     *        never modified after first use.
     */
    static RCP<const DiodeModel> defaultModel();

protected:
    /**
     * @brief Tracks BVGiven under either spelling ("BV" or "VB").
     */
    void onParamGiven(const std::string& key) override {
        if (key == "BV" || key == "VB") BVGiven = true;
    }

private:
    DiodeThermal thermal_;

    void applyGuards();
    DiodeThermal scaleTo(double tempKelvin) const;
    double fitBreakdown(double tempKelvin, double vt, double isat) const;
};

} // namespace e_XSpice

#endif // DIODE_MODEL_H

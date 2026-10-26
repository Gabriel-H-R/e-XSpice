/**
 * @author Gabriel Henrique Silva
 * @file MosfetModel.h
 * @brief Shared MOSFET model card — level 1 (Shichman-Hodges), NMOS + PMOS.
 *
 * ================================================================
 * CURRENT MODEL — Shichman-Hodges (SPICE3 MOS1 / Xyce level 1)
 * ================================================================
 * N-oriented voltages (PMOS: the device remaps with p = −1, the card is
 * polarity-agnostic — same trick as the BJT):
 *   vgs = v_g − v_s,  vds = v_d − v_s,  vbs = v_b − v_s
 * (when vds < 0 source and drain swap roles; the device handles it, the
 * card below always sees vds ≥ 0)
 *
 *   BETA = KP·(W/L_eff),   L_eff = L − 2·LD
 *   Vth  = VTO + GAMMA·(√(max(PHI − vbs, 1e-3)) − √PHI)
 *
 *   cutoff:     vgs ≤ Vth:  Id = 0
 *   linear:     vds < vgs − Vth:
 *     Id = BETA·((vgs − Vth)·vds − vds²/2)·(1 + LAMBDA·vds)
 *   saturation: Id = (BETA/2)·(vgs − Vth)²·(1 + LAMBDA·vds)
 *
 *   gm   = ∂Id/∂vgs,  gds = ∂Id/∂vds,  gmbs = ∂Id/∂vbs (via ∂Vth/∂vbs)
 *   with the exact region derivatives (no chordal approximations), so
 *   Newton keeps quadratic convergence inside each region.
 *
 * Bulk diodes (always present, SPICE semantics):
 *   Ibs = IS·(exp(vbs/vt) − 1) + gmin·vbs   (bulk → source when forward)
 *   Ibd = IS·(exp(vbd/vt) − 1) + gmin·vbd   (bulk → drain when forward)
 * Terminal currents (positive INTO the device, N-oriented):
 *   ig = 0,  id = Id + Ibd,  is = −Id + Ibs,  ib = −(Ibs + Ibd)
 *
 * ================================================================
 * CHARGE MODEL — Meyer (piecewise) + overlap + bulk depletion
 * ================================================================
 * Intrinsic Meyer partition of Cox = ε_ox/TOX·W·L_eff (×M):
 *   cutoff:     Cgb = Cox,            Cgs = Cgd = 0
 *   linear:     Cgs = Cgd = Cox/2,    Cgb = 0
 *   saturation: Cgs = 2·Cox/3,        Cgd = 0, Cgb = 0
 * plus overlap CGSO·W / CGDO·W / CGBO·L (×M). Charges are C·V by
 * construction (Qg = Cgs·vgs + Cgd·vgd + Cgb·vgb, ...), which is EXACT
 * inside each region (the Meyer fractions are piecewise-constant, so
 * dQ/dV = C with no V·dC/dV remainder) and charge-conserving by KCL.
 * Bulk depletion uses the shared junction model (Junction.h, exact Q/C).
 * Three states: Qg, Qd, Qb (Qs = −Qg−Qd−Qb, never stored).
 *
 * ================================================================
 * TEMPERATURE (computed ONCE per card, never per instance)
 * ================================================================
 * IS(T) via the SPICE law (drives the bulk diodes); PHI(T)/PB(T) via the
 * bandgap shift. VTO/KP/LAMBDA/GAMMA are temperature-independent
 * (documented level-1 simplification).
 *
 * NOT in this card (roadmap): RD/RS/RG/RB ohmic resistances, Narrow-width
 * (DELTA), short-channel (VMAX/THETA), subthreshold (NFS), noise, AC-only
 * excess parameters. LEVEL other than 1 is rejected with a clear message.
 */

#ifndef MOSFET_MODEL_H
#define MOSFET_MODEL_H

#include "../core/DeviceModel.h"
#include "../core/ModelCard.h"
#include "Junction.h"
#include <string>

namespace e_XSpice {

/**
 * @brief Temperature-scaled derived quantities of the card (plain data).
 */
struct MosfetThermal {
    double vt     = 0.0;   // kT/q
    double isat   = 0.0;   // IS(T) bulk saturation current
    double phiT   = 0.0;   // PHI(T) surface potential
    double pbT    = 0.0;   // PB(T) bulk-junction potential
    double coxPrime = 0.0;  // ε_ox/TOX (F/m²)
    double vcrit  = 0.0;   // PNJLIM critical voltage for the bulk diodes
    double gmin   = 1e-12;
    JunctionCap bulkCap{};  // bulk-junction capacitance constants (PB/MJ/FC)
};

/**
 * @brief Per-instance geometry/overrides (one netlist line).
 *
 * Own option table (same machinery as card parameters): the parser fills
 * it generically. Range validation stays in the Mosfet constructor.
 */
struct MosfetInstanceParams {
    double W = 100e-6;   // channel width (m)
    double L = 100e-6;   // channel length (m)
    double M = 1.0;      // multiplicity (parallel devices)
    double AD = 0.0;     // drain diffusion area (m², documented/future)
    double AS = 0.0;     // source diffusion area (m², documented/future)
    double NF = 1.0;     // number of fingers (Xyce NF)
    double PD = 0.0;     // drain perimeter (m, Xyce PD)
    double PS = 0.0;     // source perimeter (m, Xyce PS)
    double NRD = 1.0;    // drain squares (Xyce NRD)
    double NRS = 1.0;    // source squares (Xyce NRS)
    double SA = 0.0;     // distance to STI, source side (Xyce SA)
    double SB = 0.0;     // distance to STI, drain side (Xyce SB)
    double SD = 0.0;     // distance to STI, end (Xyce SD; 2*DMCG if not given)
    double SCA = 0.0;    // stress, source (Xyce SCA)
    double SCB = 0.0;    // stress, drain (Xyce SCB)
    double SCC = 0.0;    // stress, channel (Xyce SCC)
    double SC = 0.0;     // stress, combined (Xyce SC)
    double TEMP = 0.0;   // instance temperature (Kelvin, when TEMPGiven)
    double DTEMP = 0.0;  // delta temperature (K)
    bool   TEMPGiven = false;
    double RBDB = 0.0;     // body resistance (drain-bulk)
    double RBSB = 0.0;     // body resistance (source-bulk)
    double RBPB = 0.0;     // body resistance (bulk-bulk)
    double RBPS = 0.0;     // body resistance (bulk-source)
    double RBPD = 0.0;     // body resistance (bulk-drain)
    double DELVTO = 0.0;   // zero-bias threshold voltage variation
    double XGW = 0.0;      // distance from gate contact to device edge
    int    NGCON = 1;      // number of gate contacts
    double IC1 = 0.0;      // initial Vds
    double IC2 = 0.0;      // initial Vgs
    double IC3 = 0.0;      // initial Vbs
    int    TRNQSMOD = 0;   // transient NQS model
    int    ACNQSMOD = 0;   // AC NQS model
    int    RBODYMOD = 0;   // body resistance model (instance-level)
    int    RGATEMOD = 0;   // gate resistance model (instance-level)
    int    GEOMOD = 0;     // geometry-dependent parasitics model
    int    RGEOMOD = 0;    // S/D resistance and contact model
    int    MIN = 0;        // minimize D or S
    int    OFF = 0;        // device initially off
    int    useBerkeleyKernel = 0;  // use Berkeley C kernel via adapter (0/1)

    using Table = std::vector<TypedParamDesc<MosfetInstanceParams>>;
    static const Table& paramTable() {
        static const Table t = {
            makeParam("W", &MosfetInstanceParams::W, "channel width (m)"),
            makeParam("L", &MosfetInstanceParams::L, "channel length (m)"),
            makeParam("M", &MosfetInstanceParams::M,
                      "multiplicity (parallel devices)"),
            makeParam("AD", &MosfetInstanceParams::AD,
                      "drain diffusion area (m², future)"),
            makeParam("AS", &MosfetInstanceParams::AS,
                      "source diffusion area (m², future)"),
            makeParam("NF", &MosfetInstanceParams::NF, "number of fingers"),
            makeParam("PD", &MosfetInstanceParams::PD, "drain perimeter (m)"),
            makeParam("PS", &MosfetInstanceParams::PS, "source perimeter (m)"),
            makeParam("NRD", &MosfetInstanceParams::NRD, "drain squares"),
            makeParam("NRS", &MosfetInstanceParams::NRS, "source squares"),
            makeParam("SA", &MosfetInstanceParams::SA, "distance to STI (source)"),
            makeParam("SB", &MosfetInstanceParams::SB, "distance to STI (drain)"),
            makeParam("SD", &MosfetInstanceParams::SD, "distance to STI (end)"),
            makeParam("SCA", &MosfetInstanceParams::SCA, "stress (source)"),
            makeParam("SCB", &MosfetInstanceParams::SCB, "stress (drain)"),
            makeParam("SCC", &MosfetInstanceParams::SCC, "stress (channel)"),
            makeParam("SC", &MosfetInstanceParams::SC, "stress (combined)"),
            makeCelsiusParam("TEMP", &MosfetInstanceParams::TEMP,
                             "instance temperature",
                             [](MosfetInstanceParams& o) { o.TEMPGiven = true; }),
            makeParam("DTEMP", &MosfetInstanceParams::DTEMP,
                      "delta temperature (K)"),
            makeIntParam("USEBERKELEYKERNEL", &MosfetInstanceParams::useBerkeleyKernel,
                         "use Berkeley C kernel via BSim4Adapter (0/1)"),
            makeParam("RBDB", &MosfetInstanceParams::RBDB, "body resistance drain-bulk"),
            makeParam("RBSB", &MosfetInstanceParams::RBSB, "body resistance source-bulk"),
            makeParam("RBPB", &MosfetInstanceParams::RBPB, "body resistance bulk-bulk"),
            makeParam("RBPS", &MosfetInstanceParams::RBPS, "body resistance bulk-source"),
            makeParam("RBPD", &MosfetInstanceParams::RBPD, "body resistance bulk-drain"),
            makeParam("DELVTO", &MosfetInstanceParams::DELVTO, "threshold voltage variation"),
            makeParam("XGW", &MosfetInstanceParams::XGW, "gate contact distance"),
            makeIntParam("NGCON", &MosfetInstanceParams::NGCON, "number of gate contacts"),
            makeParam("IC1", &MosfetInstanceParams::IC1, "initial Vds"),
            makeParam("IC2", &MosfetInstanceParams::IC2, "initial Vgs"),
            makeParam("IC3", &MosfetInstanceParams::IC3, "initial Vbs"),
            makeIntParam("TRNQSMOD", &MosfetInstanceParams::TRNQSMOD, "transient NQS model"),
            makeIntParam("ACNQSMOD", &MosfetInstanceParams::ACNQSMOD, "AC NQS model"),
            makeIntParam("RBODYMOD", &MosfetInstanceParams::RBODYMOD, "body resistance model"),
            makeIntParam("RGATEMOD", &MosfetInstanceParams::RGATEMOD, "gate resistance model"),
            makeIntParam("GEOMOD", &MosfetInstanceParams::GEOMOD, "geometry parasitics model"),
            makeIntParam("RGEOMOD", &MosfetInstanceParams::RGEOMOD, "S/D resistance model"),
            makeIntParam("MIN", &MosfetInstanceParams::MIN, "minimize D or S"),
            makeIntParam("OFF", &MosfetInstanceParams::OFF, "device initially off"),
        };
        return t;
    }
};

class MosfetModel : public DeviceModel {
public:
    /**
     * @param nmos true → NMOS card, false → PMOS card
     */
    MosfetModel(std::string name, bool nmos);

    std::string typeName() const override { return nmos_ ? "NMOS" : "PMOS"; }

    /**
     * @brief true for an NMOS card (the device flips the terminal signs).
     */
    bool isNmos() const { return nmos_; }

    // ---- card parameters (SPICE3 MOS1 defaults) ----
    int    LEVEL  = 1;        // only 1 supported (Shichman-Hodges)
    double VTO    = 0.0;      // zero-bias threshold voltage (N-oriented)
    double KP     = 2e-5;      // transconductance parameter (A/V²)
    double LAMBDA = 0.0;      // channel-length modulation (1/V)
    double PHI    = 0.6;      // surface potential (V)
    double GAMMA  = 0.0;      // bulk threshold parameter (V^1/2)
    double CBD    = 0.0;      // zero-bias bulk-drain capacitance (total, F)
    double CBS    = 0.0;      // zero-bias bulk-source capacitance (total, F)
    double IS     = 1e-14;    // bulk-junction saturation current (A)
    double PB     = 0.8;      // bulk-junction potential (V)
    double MJ     = 0.5;      // bulk-junction grading coefficient
    double FC     = 0.5;      // bulk-junction forward-bias coefficient
    double CGSO   = 0.0;      // gate-source overlap capacitance per width
    double CGDO   = 0.0;      // gate-drain overlap capacitance per width
    double CGBO   = 0.0;      // gate-bulk overlap capacitance per length
    double TOX    = 1e-7;     // oxide thickness (m)
    double LD     = 0.0;      // lateral diffusion (m, Leff = L − 2·LD)
    double EG     = 1.11;     // bandgap voltage
    double XTI    = 3.0;      // saturation-current temperature exponent
    double TNOM   = 300.15;   // Kelvin
    double TEMP   = 300.15;   // Kelvin
    double GMIN   = 1e-12;    // parallel junction conductance

    bool processParams(std::string& message) override;
    void updateTemperature(double tempKelvin) override;
    MosfetThermal computeThermal(double tempKelvin) const;
    const MosfetThermal& thermal() const { return thermal_; }

    /**
     * @brief Oxide capacitance per unit area ε_ox/TOX (F/m²).
     */
    double coxPrime() const;

    static RCP<const MosfetModel> defaultNmos();
    static RCP<const MosfetModel> defaultPmos();

private:
    bool nmos_;
    MosfetThermal thermal_;

    void applyGuards();
    MosfetThermal scaleTo(double tempKelvin) const;
};

} // namespace e_XSpice

#endif // MOSFET_MODEL_H

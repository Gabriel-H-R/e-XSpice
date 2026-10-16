/**
 * @author Gabriel Henrique Silva
 * @file MosfetBSIM4Model.h
 * @brief Shared MOSFET model card — BSIM4 4.8.3 (level 2), NMOS + PMOS.
 *
 * O card carrega os parâmetros do .MODEL via tabela declarativa (mesma
 * infraestrutura dos demais modelos). As equações de corrente, mobilidade
 * e tensão de limiar são implementadas em MosfetBSIM4Device, seguindo as
 * fórmulas do BSIM4 4.8.3 (Berkeley).
 *
 * Os parâmetros aqui listados são os usados pela implementação de corrente
 * contínua (DC). Parâmetros de ruído, AC e efeitos secundários podem ser
 * adicionados posteriormente na mesma tabela, sem tocar no parser.
 */

#ifndef MOSFET_BSIM4_MODEL_H
#define MOSFET_BSIM4_MODEL_H

#include "../core/DeviceModel.h"
#include "../core/ModelCard.h"
#include "Junction.h"
#include <string>
#include <vector>

namespace e_XSpice {

/**
 * @brief Quantidades escaladas por temperatura do card BSIM4 (plain data).
 * Computadas UMA VEZ por card (setup), nunca dentro de um iteration.
 */
struct MosfetBSIM4Thermal {
    double vt       = 0.0;   // kT/q
    double eg       = 0.0;   // bandgap voltage (T)
    double ni       = 0.0;   // intrinsic carrier concentration (T)
    double u0       = 0.0;   // zero-bias mobility (T)
    double vsat     = 0.0;   // saturation velocity (T)
    double k1ox     = 0.0;   // body-effect coefficient (scaled)
    double kt1      = 0.0;   // Vth temperature coefficient
    double kt1l     = 0.0;
    double kt2      = 0.0;   // Vth temperature coefficient (Vb)
    double ua       = 0.0;   // 1st-order mobility degradation
    double ub       = 0.0;   // 2nd-order mobility degradation
    double uc       = 0.0;   // velocity saturation (mobility-dependent)
    double ud       = 0.0;   // saturation effect
    double pvag       = 0.0;   // Early-like output conductance coeff.
    double theta0vb0  = 0.0;   // mobility temperature coefficient
    double thetaRout  = 0.0;   // output conductance temperature coefficient
    double dwg      = 0.0;   // width offset (gate)
    double dwb      = 0.0;   // width offset (bulk)
    double w0       = 0.0;   // narrow width offset
    double lpe0     = 0.0;   // length offset
    double lpeb     = 0.0;   // length offset (bulk)
    double keta       = 0.0;   // back-bias threshold coeff. (linear)
    double ketac      = 0.0;   // back-bias threshold coeff. (saturation)
    double dsub       = 0.0;   // drain-induced barrier lowering coeff.
    double dvt0       = 0.0;   // short-channel threshold coeff. 1
    double dvt1       = 0.0;   // short-channel threshold coeff. 2
    double dvt2       = 0.0;   // short-channel threshold coeff. 3
    double dvt2w      = 0.0;   // short-channel threshold coeff. width
    double dvtp0      = 0.0;   // pocket implant threshold coeff.
    double dvtp1      = 0.0;
    double a0       = 0.0;   // bulk charge coeff. 1
    double ags      = 0.0;   // gate-bias coeff. (Abulk)
    double a1       = 0.0;   // nonlinearity coefficient
    double a2       = 0.0;
    double k3       = 0.0;   // narrow width coefficient
    double k3b      = 0.0;   // narrow width coefficient (Vb)
    double wl       = 0.0;   // width-length coefficient
    double wlc      = 0.0;   // width-length coefficient (C)
    double wlr      = 0.0;   // width-length coefficient (R)
    double b0       = 0.0;   // bulk charge coefficient
    double b1       = 0.0;
    double cgsl     = 0.0;   // overlap exponential coefficient
    double cgd      = 0.0;   // overlap exponential coefficient
    double cgdl     = 0.0;
    double alpha0   = 0.0;   // substrate current coefficient
    double alpha1   = 0.0;
    double beta0    = 0.0;   // substrate current exponential coefficient
    double prt      = 0.0;   // Rds temperature coefficient
    double ute      = 0.0;   // mobility temperature exponent
    double ua1      = 0.0;   // ua temperature coefficient
    double ub1      = 0.0;
    double uc1      = 0.0;
    double ud1      = 0.0;
    double prwb     = 0.0;   // Rds body-bias coefficient
    double prwg     = 0.0;   // Rds gate-bias coefficient
    double wr       = 0.0;   // width offset (Rds)

    double gmin     = 1e-12; // parallel conductance
    double isat     = 1e-14; // bulk saturation current (per area)
    JunctionCap bulkCap;     // bulk junction C-V (per area)
    double vcrit    = 0.5;   // PNJLIM critical voltage
};

/**
 * @brief Per-instance overrides (one netlist line).
 */
struct MosfetBSIM4InstanceParams {
    double AREA = 1.0;   // multiplicity
    double TEMP = 0.0;   // instance temperature (Kelvin, when TEMPGiven)
    double DTEMP = 0.0;  // delta temperature (K)
    bool   TEMPGiven = false;
};

class MosfetBSIM4Model : public DeviceModel {
public:
    /**
     * @param nmos true → NMOS card, false → PMOS card
     */
    MosfetBSIM4Model(std::string name, bool nmos);

    std::string typeName() const override { return nmos_ ? "NMOS" : "PMOS"; }
    bool isNmos() const { return nmos_; }
    bool isBSIM4() const override { return true; }

    // ---- card parameters (BSIM4 4.8.3 defaults) ----
    int    LEVEL = 2;        // only 2 supported
    double VTH0 = 0.0;       // zero-bias threshold voltage (N-oriented)
    double KP   = 2e-5;      // transconductance parameter
    double LAMBDA = 0.0;     // channel-length modulation (Early)
    double PHI  = 0.6;       // surface potential
    double GAMMA = 0.0;      // bulk threshold parameter
    double CBD  = 0.0;       // zero-bias bulk-drain capacitance (total)
    double CBS  = 0.0;       // zero-bias bulk-source capacitance (total)
    double IS   = 1e-14;     // bulk saturation current
    double PB   = 0.8;       // bulk junction potential
    double MJ   = 0.5;       // bulk junction grading coefficient
    double FC   = 0.5;       // bulk forward-bias coefficient
    double CGSO = 0.0;       // gate-source overlap capacitance per width
    double CGDO = 0.0;       // gate-drain overlap capacitance per width
    double CGBO = 0.0;       // gate-bulk overlap capacitance per length
    double TOX  = 1e-7;      // oxide thickness (m)
    double LD   = 0.0;       // lateral diffusion (m)
    double EG   = 1.11;      // bandgap voltage
    double XTI  = 3.0;       // saturation-current temperature exponent
    double TNOM = 300.15;    // Kelvin
    double TEMP = 300.15;    // Kelvin
    double GMIN = 1e-12;     // parallel conductance

    // BSIM4-specific (level 2)
    double U0   = 0.05;      // zero-bias mobility (m^2/Vs)
    double VSAT = 8e4;       // saturation velocity (m/s)
    double UA   = 1e-9;      // 1st-order mobility degradation (m/V)
    double UB   = 7.8e-19;   // 2nd-order mobility degradation (m^2/V^2)
    double UC   = -4.9e-11;  // saturation velocity coefficient (m/V)
    double UD   = 0.0;       // saturation effect coefficient
    double A0   = 1.0;       // bulk charge coefficient 1
    double B0   = 0.0;       // bulk charge coefficient 2
    double B1   = 0.0;       // bulk charge coefficient 3
    double KETA = -0.047;    // back-bias threshold coefficient
    double ALPHA0 = 0.0;     // substrate current coefficient
    double ALPHA1 = 0.0;     // substrate current exponential coefficient
    double PRT  = 0.0;       // Rds temperature coefficient
    double UTE  = -1.5;      // mobility temperature exponent
    double DVT0 = 2.2;       // short-channel threshold coefficient 1
    double DVT1 = 0.0;       // short-channel threshold coefficient 2
    double DVT2 = 0.0;       // short-channel threshold coefficient 3
    double DVT2W = 0.0;      // short-channel threshold coefficient width
    double DVTP0 = 0.0;      // pocket implant threshold coefficient
    double W0   = 0.0;       // narrow width offset
    double NFS  = 0.0;       // subthreshold slope coeff. (optional)

    bool processParams(std::string& message) override;
    void updateTemperature(double tempKelvin) override;
    MosfetBSIM4Thermal computeThermal(double tempKelvin) const;
    const MosfetBSIM4Thermal& thermal() const { return thermal_; }

    static RCP<DeviceModel> fromCard(const std::string& name,
                                     const ParamMap& params,
                                     bool nmos);
    static RCP<const MosfetBSIM4Model> defaultNmos();
    static RCP<const MosfetBSIM4Model> defaultPmos();

private:
    bool nmos_;
    MosfetBSIM4Thermal thermal_;

    void applyGuards();
    MosfetBSIM4Thermal scaleTo(double tempKelvin) const;
};

} // namespace e_XSpice

#endif // MOSFET_BSIM4_MODEL_H

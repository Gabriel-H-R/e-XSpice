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
#include <set>
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

    // ---- card parameters (BSIM4 4.8.2 addPar defaults where Xyce has the name) ----
    // Legacy names with no Xyce addPar (kept for the current device; step 2
    // replaces them): KP, PHI, GAMMA, CBD, CBS, IS, PB, MJ, FC, TAUN, NFS.
    // Xyce analogs with another name: TOX (TOXE), LD (LINT), EG (BG0SUB),
    // XTI (XTIS/XTID), GMIN (GBMIN). LAMBDA is the Xyce velocity-overshoot
    // name but means channel-length modulation here.
    int    LEVEL = 2;        // only 2 supported
    double VTH0 = 0.0;       // zero-bias threshold voltage (N-oriented); Xyce 0.7 NMOS / -0.7 PMOS if not given
    double VFB  = -1.0;      // flat-band voltage (Xyce B4.C:1105)
    double VGSTEFFVTH = 1.0;  // VgsteffVth for mobMod 3+
    double KP   = 2e-5;      // transconductance parameter
    double LAMBDA = 0.0;     // channel-length modulation (legacy meaning)
    double RDSW   = 0.0;     // source/drain series resistance (ohm*um); kept 0: legacy RDSW/W units (Xyce 200)
    // Body resistance parameters
    double RBDB   = 50.0;    // body-drain resistance
    double RBSB   = 50.0;    // body-source resistance
    double RBPB   = 50.0;    // body-bulk resistance
    double RBPS   = 50.0;    // body-source (prime) resistance
    double RBPD   = 50.0;    // body-drain (prime) resistance
    // Gate resistance parameters
    double XRCRG1 = 12.0;    // gate resistance coefficient 1
    double XRCRG2 = 1.0;     // gate resistance coefficient 2
    // Mod flags
    int    RBODYMOD = 0;     // body resistance model (0=off)
    int    RGATEMOD = 0;     // gate resistance model (0=off)
    double PRWG   = 1.0;     // Rds gate-bias coefficient
    double PRWB   = 0.0;     // Rds body-bias coefficient (Early)
    double PHI  = 0.6;       // surface potential
    double GAMMA = 0.0;      // bulk threshold parameter (k1)
    double K2    = 0.0;      // body-effect coefficient 2
    double CBD  = 0.0;       // zero-bias bulk-drain capacitance (total)
    double CBS  = 0.0;       // zero-bias bulk-source capacitance (total)
    double IS   = 1e-14;     // bulk saturation current
    double PB   = 0.8;       // bulk junction potential
    double MJ   = 0.5;       // bulk junction grading coefficient
    double FC   = 0.5;       // bulk forward-bias coefficient
    double CGSO = 0.0;       // gate-source overlap capacitance per width
    double CGDO = 0.0;       // gate-drain overlap capacitance per width
    double CGBO = 0.0;       // gate-bulk overlap capacitance per length
    double TOX  = 30.0e-10;  // oxide thickness (m)
    double LD   = 0.0;       // lateral diffusion (m)
    double EG   = 1.11;      // bandgap voltage
    double XTI  = 3.0;       // saturation-current temperature exponent
    double TNOM = 300.15;    // Kelvin
    double TEMP = 300.15;    // Kelvin
    double GMIN = 1e-12;     // parallel conductance

    // BSIM4-specific (level 2)
    double U0   = 0.0;       // zero-bias mobility (m^2/Vs); Xyce 0.067 NMOS / 0.025 PMOS if not given
    double VSAT = 8e4;       // saturation velocity (m/s)
    double UA   = 0.0;       // 1st-order mobility degradation (m/V); Xyce 1e-9 if not given
    double UB   = 1.0e-19;   // 2nd-order mobility degradation (m^2/V^2)
    double UC   = 0.0;       // saturation velocity coefficient (m/V); Xyce -4.65e-11 if not given
    double UD   = 0.0;       // saturation effect coefficient
    double A0   = 1.0;       // bulk charge coefficient 1
    double B0   = 0.0;       // bulk charge coefficient 2
    double B1   = 0.0;       // bulk charge coefficient 3
    double KETA = -0.047;    // back-bias threshold coefficient
    double KT1  = -0.11;     // threshold temperature coefficient 1
    double KT2  = 0.022;     // threshold temperature coefficient 2
    double ALPHA0 = 0.0;     // substrate current coefficient
    double AGIDL  = 0.0;     // GIDL coefficient
    double BGIDL  = 2.3e9;   // GIDL exponential coefficient
    double CGIDL  = 0.5;     // GIDL body-effect coefficient
    double EGIDL  = 0.8;     // GIDL activation energy
    // Gate current parameters
    double AIGC   = 1.36e-2; // gate current coeff. (accumulation); Xyce 9.80e-3 PMOS if not given
    double AECHVB = 1.0;      // gate current prefactor
    double BECHVB = 1.0;      // gate current exponent coeff
    double BIGC   = 1.71e-3; // gate current coeff. (accumulation); Xyce 7.59e-4 PMOS if not given
    double CIGC   = 0.075;   // gate current coeff. (accumulation); Xyce 0.03 PMOS if not given
    double AIGS   = 1.36e-2; // gate current coeff. (source)
    double AIGD   = 1.36e-2; // gate current coeff. (drain)
    // Substrate (bulk) current parameters
    double AIGBACC = 1.36e-2; // bulk current coeff. (accumulation)
    double AIGBINV = 1.11e-2; // bulk current coeff. (inversion)
    double NIGC   = 1.0;     // gate current exponent
    int    IGCMod = 0;       // gate current model (0=off)
    int    IGBMod = 0;       // bulk current model (0=off)
    // NQS parameters
    int    TRNQSMOD = 0;     // transient NQS model (0=off)
    int    ACNQSMOD = 0;     // AC NQS model (0=off)
    double TAUN = 0.0;       // NQS time constant (s)
    // Binning parameters
    double LMIN = 0.0;       // minimum channel length for binning
    double LMAX = 1.0;       // maximum channel length for binning
    double WMIN = 0.0;       // minimum channel width for binning
    double WMAX = 1.0;       // maximum channel width for binning
    // Geometry parameters
    double DMCG = 0.0;       // distance from gate to contact
    double DMCI = 0.0;       // distance from gate to contact (inner)
    double DMDG = 0.0;       // distance from gate to contact (drain)
    double DMCGT = 0.0;      // DMCG temperature coefficient
    // Stress parameters
    double SA    = 0.0;      // distance to STI (source side)
    double SB    = 0.0;      // distance to STI (drain side)
    double SD    = 0.0;      // distance to STI (end)
    double SCA   = 0.0;      // stress (source)
    double SCB   = 0.0;      // stress (drain)
    double SCC   = 0.0;      // stress (channel)
    double SC    = 0.0;      // stress (combined)
    // Well proximity parameters
    double KVTH0WE = 0.0;    // well proximity Vth coefficient
    double K2WE    = 0.0;    // well proximity K2 coefficient
    double KU0WE   = 0.0;    // well proximity mobility coefficient
    int    WPEMOD  = 0;      // well proximity model (0=off)
    // Noise parameters
    double KF     = 0.0;     // flicker noise coefficient
    // Subthreshold parameters
    double NFACTOR = 1.0;    // subthreshold ideality factor
    double VOFF   = -0.08;   // threshold voltage offset
    double AF     = 1.0;     // flicker noise exponent
    double EF     = 1.0;     // flicker noise frequency exponent
    double RNOIA  = 0.577;   // thermal noise coefficient (drain)
    double NOIA   = 1e42;    // oxide trap density A (flicker noise)
    double NOIB   = 0.0;     // oxide trap density B
    double NOIC   = 0.0;     // oxide trap density C
    double EM     = 4.1e7;   // flicker noise exponent
    double RNOIB  = 0.5164;  // thermal noise coefficient (source)
    double RNOIC  = 0.395;   // thermal noise coefficient (gate)
    double AGISL  = 0.0;     // GISL coefficient
    double BGISL  = 2.3e-9;  // GISL exponential coefficient; Xyce: =BGIDL if not given
    double CGISL  = 0.5;     // GISL body-effect coefficient
    double EGISL  = 0.8;     // GISL activation energy
    double ALPHA1 = 0.0;     // substrate current exponential coefficient
    double PRT  = 0.0;       // Rds temperature coefficient
    double UTE  = -1.5;      // mobility temperature exponent
    double DVT0 = 2.2;       // short-channel threshold coefficient 1
    double DVT1 = 0.53;      // short-channel threshold coefficient 2
    double DVT2 = -0.032;    // short-channel threshold coefficient 3
    double DVT2W = -0.032;   // short-channel threshold coefficient width
    double DVTP0 = 0.0;      // pocket implant threshold coefficient
    double W0   = 2.5e-6;    // narrow width offset
    double NFS  = 0.0;       // subthreshold slope coeff. (optional)

    // ---- Xyce BSIM4 4.8.2 addPar (N_DEV_MOSFET_B4.C), non-binned ----
    // Initializers are the addPar defaults; card-dependent defaults are
    // applied in applyDefaults() (Model::processParams4p82_).
    double EOT = 15.0e-10;
    double VDDEOT = 1.5;
    double TEMPEOT = 300.15;
    double LEFFEOT = 1e-6;
    double WEFFEOT = 10e-6;
    double ADOS = 1.0;
    double BDOS = 1.0;
    double TOXE = 30.0e-10;
    double TOXP = 30.0e-10;
    double TOXM = 30.0e-10;
    double TOXREF = 30.0e-10;
    double DTOX = 0.0;
    double EPSROX = 3.9;
    double DTYPE = 1.0;
    double EG0 = 1.16;
    double VTFBPHI1 = 0.0;
    double VGSTEFFFVTH = 0.0;
    double CDSC = 2.4e-4;
    double CDSCB = 0.0;
    double CDSCD = 0.0;
    double CIT = 0.0;
    double XJ = 0.15e-6;
    double AT = 3.3e4;
    double AGS = 0.0;
    double A1 = 0.0;
    double A2 = 1.0;
    double PHIG = 4.05;
    double EPSRGATE = 11.7;
    double EASUB = 4.05;
    double EPSRSUB = 11.7;
    double NI0SUB = 1.45e10;
    double BG0SUB = 1.16;
    double TBGASUB = 7.02e-4;
    double TBGBSUB = 1108.0;
    double NSUB = 6.0e16;
    double NDEP = 1.7e17;
    double NSD = 1.0e20;
    double CKAPPAS = 0.6;   // S-side overlap smoothing (Xyce CKAPPAS)
    double CKAPPAD = 0.6;   // D-side overlap smoothing (Xyce CKAPPAD)
    double CF = 0.0;        // fringe capacitance (Xyce default when not given)
    double CLC = 0.1e-6;    // C-V length scale of abulkCV (Xyce CLC)
    double CLE = 0.6;       // C-V length exponent (Xyce CLE)
    double ACDE = 1.0;      // depletion-width factor of CTM (Xyce ACDE)
    double VFBCV = -1.0;    // flat-band voltage of capMod 0 C-V (Xyce VFBCV)
    double VOFFCV = 0.0;    // C-V offset voltage (Xyce VOFFCV)
    double VOFFCVL = 0.0;   // C-V length offset (Xyce VOFFCVL)
    double MINVCV = 0.0;    // C-V inversion-layer parameter (Xyce MINVCV)
    double DLC = 0.0;       // C-V length offset (Xyce DLC)
    double DWC = 0.0;       // C-V width offset (Xyce DWC)
    double CGDL = 0.0;      // drain-gate overlap smoothing charge (Xyce CGDL)
    double CGSL = 0.0;      // source-gate overlap smoothing charge (Xyce CGSL)
    double PHIN = 0.0;
    double NGATE = 0.0;
    double GAMMA1 = 0.0;
    double GAMMA2 = 0.0;
    double VBX = 0.0;
    double VBM = -3.0;
    double XT = 1.55e-7;
    double K1 = 0.0;
    double KT1L = 0.0;
    double K3 = 80.0;
    double K3B = 0.0;
    double DVTP1 = 0.0;
    double DVTP2 = 0.0;
    double DVTP3 = 0.0;
    double DVTP4 = 0.0;
    double DVTP5 = 0.0;
    double LPE0 = 1.74e-7;
    double LPEB = 0.0;
    double DVT0W = 0.0;
    double DVT1W = 5.3e6;
    double DROUT = 0.56;
    double DSUB = 0.0;
    double UA1 = 1.0e-9;
    double UB1 = -1.0e-18;
    double UC1 = 0.0;
    double UD1 = 0.0;
    double UP = 0.0;
    double LP = 1.0e-8;
    double EU = 0.0;
    double UCS = 1.67;
    double UCSTE = -4.775e-3;
    double MINV = 0.0;
    double VOFFL = 0.0;
    double XPART = 0.0;
    double DELTA = 0.01;
    double RSH = 0.0;
    double RDSWMIN = 0.0;
    double RS0 = 0.0;       // source/drain series resistance coeff
    double VFBSD = 0.0;     // flat-band voltage for RDS (source/drain)
    double RSW = 100.0;
    double RDW = 100.0;
    double RDWMIN = 0.0;
    double RSWMIN = 0.0;
    double ETA0 = 0.08;
    double ETAB = -0.07;
    double PCLM = 1.3;
    double PDIBLC1 = 0.39;
    double PDIBLC2 = 0.0086;
    double PDIBLCB = 0.0;
    double FPROUT = 0.0;
    double PDITS = 0.0;
    double PDITSL = 0.0;
    double PDITSD = 0.0;
    double PSCBE1 = 4.24e8;
    double PSCBE2 = 1.0e-5;
    double PVAG = 0.0;
    double JSS = 1e-4;
    double JSWS = 0.0;
    double JSWGS = 0.0;
    double PBS = 1.0;
    double NJS = 1.0;
    double XTIS = 3.0;
    double MJS = 0.5;
    double PBSWS = 1.0;
    double MJSWS = 0.33;
    double PBSWGS = 0.0;
    double MJSWGS = 0.33;
    double CJS = 5e-4;
    double CJSWS = 5e-10;
    double CJSWGS = 0.0;
    double JSD = 1e-4;
    double JSWD = 0.0;
    double JSWGD = 0.0;
    double PBD = 1.0;
    double NJD = 1.0;
    double XTID = 3.0;
    double MJD = 0.5;
    double PBSWD = 1.0;
    double MJSWD = 0.33;
    double PBSWGD = 0.0;
    double MJSWGD = 0.33;
    double CJD = 5e-4;
    double CJSWD = 5e-10;
    double CJSWGD = 0.0;
    double TPB = 0.0;
    double TCJ = 0.0;
    double TPBSW = 0.0;
    double TCJSW = 0.0;
    double TPBSWG = 0.0;
    double TCJSWG = 0.0;
    double MOIN = 15.0;
    double NOFF = 1.0;
    double XGW = 0.0;
    double XGL = 0.0;
    double RSHG = 0.1;
    double NGCON = 1.0;
    double VTL = 2.0e5;
    double LC = 5.0e-9;
    double XN = 3.0;
    double VFBSDOFF = 0.0;
    double TVFBSDOFF = 0.0;
    double TVOFF = 0.0;
    double TNFACTOR = 0.0;
    double TETA0 = 0.0;
    double LINT = 0.0;
    double LL = 0.0;
    double LLN = 1.0;
    double LW = 0.0;
    double LWN = 1.0;
    double LWL = 0.0;
    double WR = 1.0;
    double WINT = 0.0;
    double DWG = 0.0;
    double DWB = 0.0;
    double WL = 0.0;
    double WLN = 1.0;
    double WW = 0.0;
    double WWN = 1.0;
    double WWL = 0.0;
    double XW = 0.0;
    double XL = 0.0;
    double BETA0 = 0.0;
    double RGIDL = 1.0;
    double KGIDL = 0.0;
    double FGIDL = 0.0;
    double RGISL = 1.0;
    double KGISL = 0.0;
    double FGISL = 0.0;
    double AIGSD = 1.36e-2;
    double BIGSD = 1.71e-3;
    double CIGSD = 0.075;
    double BIGS = 1.71e-3;
    double CIGS = 0.075;
    double BIGD = 1.71e-3;
    double CIGD = 0.075;
    double BIGBACC = 1.71e-3;
    double CIGBACC = 0.075;
    double BIGBINV = 9.49e-4;
    double CIGBINV = 0.006;
    double NIGBINV = 3.0;
    double NIGBACC = 1.0;
    double NTOX = 1.0;
    double EIGBINV = 1.1;
    double PIGCD = 1.0;
    double POXEDGE = 1.0;
    double IJTHDFWD = 0.1;
    double IJTHSFWD = 0.1;
    double IJTHDREV = 0.1;
    double IJTHSREV = 0.1;
    double XJBVD = 1.0;
    double XJBVS = 1.0;
    double BVD = 10.0;
    double BVS = 10.0;
    double JTSS = 0.0;
    double JTSD = 0.0;
    double JTSSWS = 0.0;
    double JTSSWD = 0.0;
    double JTSSWGS = 0.0;
    double JTSSWGD = 0.0;
    double JTWEFF = 0.0;
    double NJTS = 20.0;
    double NJTSSW = 20.0;
    double NJTSSWG = 20.0;
    double NJTSD = 20.0;
    double NJTSSWD = 20.0;
    double NJTSSWGD = 20.0;
    double XTSS = 0.02;
    double XTSD = 0.02;
    double XTSSWS = 0.02;
    double XTSSWD = 0.02;
    double XTSSWGS = 0.02;
    double XTSSWGD = 0.02;
    double TNJTS = 0.0;
    double TNJTSSW = 0.0;
    double TNJTSSWG = 0.0;
    double TNJTSD = 0.0;
    double TNJTSSWD = 0.0;
    double TNJTSSWGD = 0.0;
    double VTSS = 10.0;
    double VTSD = 10.0;
    double VTSSWS = 10.0;
    double VTSSWD = 10.0;
    double VTSSWGS = 10.0;
    double VTSSWGD = 10.0;
    double GBMIN = 1.0e-12;
    double RBPS0 = 50.0;
    double RBPSL = 0.0;
    double RBPSW = 0.0;
    double RBPSNF = 0.0;
    double RBPD0 = 50.0;
    double RBPDL = 0.0;
    double RBPDW = 0.0;
    double RBPDNF = 0.0;
    double RBPBX0 = 100.0;
    double RBPBXL = 0.0;
    double RBPBXW = 0.0;
    double RBPBXNF = 0.0;
    double RBPBY0 = 100.0;
    double RBPBYL = 0.0;
    double RBPBYW = 0.0;
    double RBPBYNF = 0.0;
    double RBSBX0 = 100.0;
    double RBSBY0 = 100.0;
    double RBDBX0 = 100.0;
    double RBDBY0 = 100.0;
    double RBSDBXL = 0.0;
    double RBSDBXW = 0.0;
    double RBSDBXNF = 0.0;
    double RBSDBYL = 0.0;
    double RBSDBYW = 0.0;
    double RBSDBYNF = 0.0;
    double SAREF = 1.0e-6;
    double SBREF = 1.0e-6;
    double WLOD = 0.0;
    double KU0 = 0.0;
    double KVSAT = 0.0;
    double KVTH0 = 0.0;
    double LKVTH0 = 0.0;
    double WKVTH0 = 0.0;
    double PKVTH0 = 0.0;
    double TKU0 = 0.0;
    double LLODKU0 = 0.0;
    double WLODKU0 = 0.0;
    double LLODVTH = 0.0;
    double WLODVTH = 0.0;
    double STK2 = 0.0;
    double LODK2 = 1.0;
    double STETA0 = 0.0;
    double LODETA0 = 1.0;
    double WEB = 0.0;
    double WEC = 0.0;
    double SCREF = 1.0e-6;
    double GIDLCLAMP = -1.0e-5;
    double IDOVVDS = 1e-9;
    int CVCHARGEMOD = 0;
    int CAPMOD = 2;
    int DIOMOD = 1;
    int RDSMOD = 0;
    int MOBMOD = 0;
    int PERMOD = 1;
    int GEOMOD = 0;
    int RGEOMOD = 0;
    int FNOIMOD = 1;
    int TNOIMOD = 0;
    int MTRLMOD = 0;
    int MTRLCOMPATMOD = 0;
    int TEMPMOD = 0;
    int GIDLMOD = 0;
    int PARAMCHK = 1;
    int BINUNIT = 1;

    bool processParams(std::string& message) override;
    void updateTemperature(double tempKelvin) override;
    MosfetBSIM4Thermal computeThermal(double tempKelvin) const;
    const MosfetBSIM4Thermal& thermal() const { return thermal_; }

    static RCP<DeviceModel> fromCard(const std::string& name,
                                     const ParamMap& params,
                                     bool nmos);
    static RCP<const MosfetBSIM4Model> defaultNmos();
    static RCP<const MosfetBSIM4Model> defaultPmos();

public:
    /** @brief True when the card explicitly gave `key` (Xyce `given` flag). */
    bool isGiven(const std::string& key) const {
        return given_.count(key) != 0;
    }

protected:
    /**
     * @brief Records every card key given explicitly (drives the
     *        Xyce `given` defaults in applyDefaults).
     */
    void onParamGiven(const std::string& key) override;

private:
    bool nmos_;
    MosfetBSIM4Thermal thermal_;
    std::set<std::string> given_;

    void applyDefaults();
    void checkAndFixVersion();
    bool processInstanceParams();
    void setupVersionPointers();
    void clearTemperatureData();
    void printOutInstances();
    void applyGuards();
    MosfetBSIM4Thermal scaleTo(double tempKelvin) const;
};

} // namespace e_XSpice

#endif // MOSFET_BSIM4_MODEL_H

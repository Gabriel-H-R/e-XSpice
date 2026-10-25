/**
 * @author Gabriel Henrique Silva
 * @file BSim4Adapter.cpp
 * @brief BSIM4 4.8.3 adapter implementation — wraps the Berkeley C load
 *        functions (BSIM4load, BSIM4temp, ...) into e-XSpice's DAE contract.
 *
 * The Berkeley C code is compiled verbatim (as C) from
 * src/models/bsim4/code/. This adapter:
 *
 *   1. Allocates and populates BSIM4model + BSIM4instance from the e-XSpice
 *      model card (MosfetBSIM4Model).
 *   2. Keeps a private bsim4_state_ buffer (NSTATVARS doubles) that the
 *      Berkeley code reads/writes through CKTstate0 — scratch storage for
 *      terminal voltages (vbd, vbs, vgs, vds, ...), NOT the solution vector x.
 *   3. Builds a minimal CKTcircuit (SpiceContext.h) whose CKTstate0 /
 *      CKTrhs / CKTmatrix borrow e-XSpice's vectors/matrix.
 *   4. On each evaluation, copies terminal voltages from x into bsim4_state_
 *      and calls BSIM4load() — unchanged Berkeley code.
 *   5. The SMPaddElt shim translates matrix stamps into pre-resolved slots.
 *
 * Physics is 100% Berkeley: not a single line of b4ld.c / b4temp.c is
 * modified. The adapter only translates interfaces.
 */

#include "BSim4Adapter.h"

#include "../../solver/DaeMatrix.h"
#include "../../core/DaeSolution.h"
#include "../../core/Node.h"
#include "../../models/core/DeviceConstants.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace e_XSpice {

// =========================================================================
//  SPICE stub implementations (C linkage — called by Berkeley C code)
// =========================================================================

extern "C" {

int NIintegrate(void * /*ckt*/, double * /*geq*/, double * /*ceq*/,
                double /*cf*/, int /*state*/) {
    return 0;  // no charge modeling in level-2 subset
}

int DEVpnjlim(double vnew, double vold, double vt,
              double vcrit, double *vlimited) {
    *vlimited = vnew;  // no limiting in level-2 subset
    return 0;
}

int BSIM4PAeffGeo(int /*nf*/, int /*geoMod*/, int /*min*/,
                  double *pseff, double *adeff, double *pdeff) {
    if (pseff) *pseff = 1.0;
    if (adeff) *adeff = 1.0;
    if (pdeff) *pdeff = 1.0;
    return 0;
}

int BSIM4RdseffGeo(int /*nf*/, int /*geoMod*/, int /*min*/,
                   double *wr, double *dwg, double *dwb,
                   double *prwb, double *prwg) {
    if (wr) *wr = 0.0;
    if (dwg) *dwg = 0.0;
    if (dwb) *dwb = 0.0;
    if (prwb) *prwb = 0.0;
    if (prwg) *prwg = 0.0;
    return 0;
}

int BSIM4checkModel(BSIM4model * /*model*/,
                    BSIM4instance * /*inst*/,
                    CKTcircuit * /*ckt*/) {
    return OK;  // always OK
}

/* DEVfetlim / DEVlimvds: SPICE limiting functions (stubs —
 * level-2 subset handles limiting in the Mosfet wrapper). */
int DEVfetlim(double /*vnew*/, double /*vold*/, double /*vt*/,
              double /*vcrit*/, double *vlimited) {
    *vlimited = 0.0;
    return 0;
}

int DEVlimvds(double /*vnew*/, double /*vold*/, double /*vt*/,
              double /*vcrit*/, double *vlimited) {
    *vlimited = 0.0;
    return 0;
}

static struct SPfrontEnd kSpFrontEndStub = {nullptr, nullptr};
struct SPfrontEnd *SPfrontEnd = &kSpFrontEndStub;

} // extern "C"

// =========================================================================
//  Constructor / destructor
// =========================================================================

BSim4Adapter::BSim4Adapter(RCP<const MosfetBSIM4Model> model,
                           const std::vector<Node*>& nodes,
                           double W, double L, double M, double area)
    : SpiceDeviceAdapter("BSIM4", DeviceType::MOSFET, nodes),
      model_(std::move(model)),
      W_(W), L_(L), M_(M), area_(area),
      p_(1),
      bsim4_state_(NSTATVARS, 0.0),  // NSTATVARS state variables per instance
      g_cache_() {

    if (!model_) {
        throw std::invalid_argument("BSIM4 adapter requires a model card");
    }
    if (nodes_.size() < 4) {
        throw std::invalid_argument(
            "BSIM4 MOSFET needs 4 nodes: drain, gate, source, bulk");
    }

    // Determine polarity from the model card.
    p_ = model_->isNmos() ? +1 : -1;

    // Allocate the Berkeley structures (owned by this adapter).
    bsim4_model_ = new BSIM4model();
    bsim4_inst_ = new BSIM4instance();
    std::memset(bsim4_model_, 0, sizeof(BSIM4model));
    std::memset(bsim4_inst_, 0, sizeof(BSIM4instance));
    buildModelCard();
    buildInstance();
    setupStates();
    setup();              // <<< full BSIM4setup(): defaults + pParam + instance fields
    // Link model -> instance so BSIM4load can find this instance
    bsim4_model_->BSIM4instances = bsim4_inst_;
    setup_done_ = true;
}

BSim4Adapter::~BSim4Adapter() {
    delete pParam_;
    delete bsim4_inst_;
    delete bsim4_model_;
}

// =========================================================================
//  Build the Berkeley model card from the e-XSpice card
// =========================================================================

void BSim4Adapter::buildModelCard() {
    auto* m = bsim4_model_;
    const auto& card = *model_;

    // Device type
    m->BSIM4modType = model_->isNmos() ? NMOS : PMOS;
    m->BSIM4type = model_->isNmos() ? NMOS : PMOS;
    m->BSIM4typeGiven = 1;

    // Flags — enable the standard (level-2-compatible) subset.
    m->BSIM4mobModGiven = 1;   m->BSIM4mobMod = 0;
    m->BSIM4cvchargeModGiven = 1; m->BSIM4cvchargeMod = 0;
    m->BSIM4capModGiven = 1;   m->BSIM4capMod = 2;
    m->BSIM4dioModGiven = 1;   m->BSIM4dioMod = 1;
    m->BSIM4rdsModGiven = 1;   m->BSIM4rdsMod = 0;
    m->BSIM4rbodyModGiven = 1; m->BSIM4rbodyMod = 0;
    m->BSIM4rgateModGiven = 1; m->BSIM4rgateMod = 0;
    m->BSIM4perModGiven = 1;   m->BSIM4perMod = 1;
    m->BSIM4geoModGiven = 1;   m->BSIM4geoMod = 0;
    m->BSIM4fnoiModGiven = 1;  m->BSIM4fnoiMod = 0;
    m->BSIM4tnoiModGiven = 1;  m->BSIM4tnoiMod = 0;
    m->BSIM4trnqsModGiven = 1; m->BSIM4trnqsMod = 0;
    m->BSIM4acnqsModGiven = 1; m->BSIM4acnqsMod = 0;
    m->BSIM4mtrlModGiven = 1;  m->BSIM4mtrlMod = 0;
    m->BSIM4igcModGiven = 1;   m->BSIM4igcMod = 0;
    m->BSIM4igbModGiven = 1;   m->BSIM4igbMod = 0;
    m->BSIM4tempModGiven = 1;  m->BSIM4tempMod = 0;

    m->BSIM4paramChkGiven = 1; m->BSIM4paramChk = 0;
    m->BSIM4versionGiven = 1;  m->BSIM4version = 4.83;
    m->BSIM4binUnitGiven = 1;  m->BSIM4binUnit = 1;

    // --- Model parameters that exist in sBSIM4model ---
    m->BSIM4vth0Given = 1;  m->BSIM4vth0 = card.VTH0;
    m->BSIM4u0Given = 1;    m->BSIM4u0 = card.U0;
    m->BSIM4vsatGiven = 1;  m->BSIM4vsat = card.VSAT;
    m->BSIM4uaGiven = 1;    m->BSIM4ua = card.UA;
    m->BSIM4ubGiven = 1;    m->BSIM4ub = card.UB;
    m->BSIM4ucGiven = 1;    m->BSIM4uc = card.UC;
    m->BSIM4udGiven = 1;    m->BSIM4ud = card.UD;
    m->BSIM4eta0Given = 1;  m->BSIM4eta0 = 0.0;
    m->BSIM4kt1Given = 1;   m->BSIM4kt1 = -0.011;
    m->BSIM4kt2Given = 1;   m->BSIM4kt2 = 0.022;
    m->BSIM4k2Given = 1;    m->BSIM4k2 = 0.05;
    m->BSIM4k3Given = 1;    m->BSIM4k3 = 0.0;
    m->BSIM4k3bGiven = 1;   m->BSIM4k3b = 0.0;
    m->BSIM4w0Given = 1;    m->BSIM4w0 = card.W0;
    m->BSIM4dvt0Given = 1;  m->BSIM4dvt0 = card.DVT0;
    m->BSIM4dvt1Given = 1;  m->BSIM4dvt1 = card.DVT1;
    m->BSIM4dvt2Given = 1;  m->BSIM4dvt2 = card.DVT2;
    m->BSIM4dvt0wGiven = 1; m->BSIM4dvt0w = card.DVT2W;
    m->BSIM4dvt1wGiven = 1; m->BSIM4dvt1w = 0.0;
    m->BSIM4dvt2wGiven = 1; m->BSIM4dvt2w = 0.0;
    m->BSIM4dsubGiven = 1;  m->BSIM4dsub = 0.0;
    m->BSIM4a0Given = 1;    m->BSIM4a0 = card.A0;
    m->BSIM4agsGiven = 1;   m->BSIM4ags = 0.0;
    m->BSIM4a1Given = 1;    m->BSIM4a1 = 0.0;
    m->BSIM4a2Given = 1;    m->BSIM4a2 = 1.0;
    m->BSIM4ketaGiven = 1;  m->BSIM4keta = card.KETA;
    m->BSIM4ketacGiven = 1; m->BSIM4ketac = card.KETA;
    m->BSIM4nsubGiven = 1;  m->BSIM4nsub = 6.0e16;
    m->BSIM4phigGiven = 1;  m->BSIM4phig = 0.6;
    m->BSIM4gamma1Given = 1; m->BSIM4gamma1 = card.GAMMA;
    m->BSIM4gamma2Given = 1; m->BSIM4gamma2 = 0.0;
    m->BSIM4vbxGiven = 1;   m->BSIM4vbx = 0.0;
    m->BSIM4vbmGiven = 1;   m->BSIM4vbm = -3.0;
    m->BSIM4xtGiven = 1;    m->BSIM4xt = 1.55e-7;
    m->BSIM4pclmGiven = 1;  m->BSIM4pclm = 1.0;
    m->BSIM4pdibl1Given = 1; m->BSIM4pdibl1 = 0.0;
    m->BSIM4pdibl2Given = 1; m->BSIM4pdibl2 = 0.0;
    m->BSIM4droutGiven = 1; m->BSIM4drout = 0.0;
    m->BSIM4pvagGiven = 1;  m->BSIM4pvag = (card.LAMBDA != 0.0) ? card.LAMBDA : 0.0;
    m->BSIM4lambdaGiven = 1; m->BSIM4lambda = 0.0;
    m->BSIM4wrGiven = 1;    m->BSIM4wr = 0.0;
    m->BSIM4dwgGiven = 1;   m->BSIM4dwg = 0.0;
    m->BSIM4dwbGiven = 1;   m->BSIM4dwb = 0.0;
    m->BSIM4b0Given = 1;    m->BSIM4b0 = card.B0;
    m->BSIM4b1Given = 1;    m->BSIM4b1 = card.B1;
    m->BSIM4alpha0Given = 1; m->BSIM4alpha0 = card.ALPHA0;
    m->BSIM4alpha1Given = 1; m->BSIM4alpha1 = card.ALPHA1;
    m->BSIM4prtGiven = 1;   m->BSIM4prt = card.PRT;
    m->BSIM4uteGiven = 1;   m->BSIM4ute = card.UTE;

    // Physical constants
    m->BSIM4epsroxGiven = 1; m->BSIM4epsrox = 3.9;
    m->BSIM4toxpGiven = 1;  m->BSIM4toxp = card.TOX;
    m->BSIM4coxe = 3.9 * 8.85418e-12 / card.TOX;

    // Overlap capacitances
    m->BSIM4cgsoGiven = 1; m->BSIM4cgso = card.CGSO;
    m->BSIM4cgdoGiven = 1; m->BSIM4cgdo = card.CGDO;
    m->BSIM4cgboGiven = 1; m->BSIM4cgbo = card.CGBO;

    // Temperature
    m->BSIM4tnomGiven = 1; m->BSIM4tnom = card.TNOM;

    // ---- GIDL/GISL ----
    // m->BSIM4Igidl = 0.0;  // BSIM4Igidl is an instance field, not model
    // m->BSIM4Igisl = 0.0;  // BSIM4Igisl is an instance field, not model
    // m->BSIM4csub = 0.0;   // not a model member
}

// =========================================================================
//  Build the Berkeley instance
// =========================================================================

void BSim4Adapter::buildInstance() {
    auto* inst = bsim4_inst_;
    auto* m = bsim4_model_;

    inst->BSIM4modPtr = m;
    inst->BSIM4nextInstance = nullptr;
    inst->BSIM4name = const_cast<char*>("BSIM4_MOSFET");

    // Instance geometry
    inst->BSIM4lGiven = 1;   inst->BSIM4l = L_;
    inst->BSIM4wGiven = 1;   inst->BSIM4w = W_;
    inst->BSIM4mult_iGiven = 1; inst->BSIM4mult_i = M_;
    inst->BSIM4mult_qGiven = 1; inst->BSIM4mult_q = M_;
    inst->BSIM4mult_fnGiven = 1; inst->BSIM4mult_fn = 1.0;
    inst->BSIM4nfGiven = 1;   inst->BSIM4nf = 1.0;
    inst->BSIM4minGiven = 1;  inst->BSIM4min = 1;
    inst->BSIM4off = 0;

    // Areas (scaled by multiplicity)
    double ad = area_ * M_;
    double as = area_ * M_;
    inst->BSIM4drainAreaGiven = 1;   inst->BSIM4drainArea = ad;
    inst->BSIM4sourceAreaGiven = 1;  inst->BSIM4sourceArea = as;
    inst->BSIM4drainSquaresGiven = 1; inst->BSIM4drainSquares = 0.0;
    inst->BSIM4sourceSquaresGiven = 1; inst->BSIM4sourceSquares = 0.0;
    inst->BSIM4drainPerimeterGiven = 1; inst->BSIM4drainPerimeter = 0.0;
    inst->BSIM4sourcePerimeterGiven = 1; inst->BSIM4sourcePerimeter = 0.0;

    inst->BSIM4sourceConductance = 0.0;
    inst->BSIM4drainConductance = 0.0;

    // ---- Instance-level diode parameters ----
    inst->BSIM4grbsb = 0.0;
    inst->BSIM4grbdb = 0.0;
    inst->BSIM4grbpb = 0.0;
    inst->BSIM4rbdb = 0.0;
    inst->BSIM4rbsb = 0.0;
    inst->BSIM4rbpb = 0.0;
    inst->BSIM4rbps = 0.0;
    inst->BSIM4rbpd = 0.0;
    inst->BSIM4cbs = 0.0;
    inst->BSIM4cbd = 0.0;

    // ---- Instance current fields (accessed by BSIM4load) ----
    inst->BSIM4Igidl = 0.0;
    inst->BSIM4Igisl = 0.0;
    inst->BSIM4Igcs = 0.0;
    inst->BSIM4Igcd = 0.0;
    inst->BSIM4Igs = 0.0;
    inst->BSIM4Igd = 0.0;
    inst->BSIM4Igb = 0.0;
    inst->BSIM4csub = 0.0;

    // ---- Instance-level junction diode parameters ----
    inst->BSIM4vjsmFwd = 0.0;
    inst->BSIM4vjsmRev = 0.0;
    inst->BSIM4vjdmFwd = 0.0;
    inst->BSIM4vjdmRev = 0.0;
    inst->BSIM4XExpBVS = 0.0;
    inst->BSIM4XExpBVD = 0.0;
    inst->BSIM4SslpFwd = 0.0;
    inst->BSIM4SslpRev = 0.0;
    inst->BSIM4DslpFwd = 0.0;
    inst->BSIM4DslpRev = 0.0;
    inst->BSIM4IVjsmFwd = 0.0;
    inst->BSIM4IVjsmRev = 0.0;
    inst->BSIM4IVjdmFwd = 0.0;
    inst->BSIM4IVjdmRev = 0.0;
}

// =========================================================================
//  State assignment
// =========================================================================

void BSim4Adapter::setupStates() {
    // Berkeley state indexing: BSIM4vds = BSIM4states + 3, etc.
    // We point CKTstate0 to our private buffer (size NSTATVARS), so
    // BSIM4states must be 0 so that BSIM4states + offset indexes
    // directly into bsim4_state_[0..NSTATVARS-1].
    bsim4_inst_->BSIM4states = 0;
}

// =========================================================================
//  Full BSIM4setup() equivalent: defaults + pParam + instance fields
// =========================================================================

void BSim4Adapter::setup() {
    auto* m = bsim4_model_;
    auto* inst = bsim4_inst_;

    // ------------------------------------------------------------------
    // 1. CKTcircuit fields needed by BSIM4load
    // ------------------------------------------------------------------
    ctx_.CKTgmin = 1.0e-12;          // parallel conductance
    ctx_.CKTtemp = m->BSIM4tnom;     // circuit temperature (Kelvin)
    ctx_.CKTmode = 0x0001;           // MODEDCOP
    ctx_.CKTbypass = 0;              // no bypass on first call
    ctx_.CKTreltol = 1.0e-3;         // relative tolerance
    ctx_.CKTvoltTol = 1.0e-4;        // voltage tolerance
    ctx_.CKTabstol = 1.0e-12;        // absolute tolerance (A)
    // CKTrhs, CKTstate0, CKTmatrix, CKTnumNodes set per-evaluation

    // ------------------------------------------------------------------
    // 2. Model-level defaults (from BSIM4setup in b4set.c)
    //    All parameters not already set by buildModelCard() get defaults.
    // ------------------------------------------------------------------

    // Physical constants
    const double eps0 = 8.85418e-12;   // F/m
    const double epsSi = 1.03594e-10;  // F/m
    const double q = 1.60219e-19;      // C
    const double k = 1.38065e-23;      // J/K
    const double Tnom = m->BSIM4tnom;  // nominal temperature
    const double Vt = k * Tnom / q;    // thermal voltage at Tnom

    // --- Defaults for flags not already set ---
    if (!m->BSIM4mobModGiven)   m->BSIM4mobMod = 0;
    if (!m->BSIM4capModGiven)   m->BSIM4capMod = 2;
    if (!m->BSIM4dioModGiven)   m->BSIM4dioMod = 1;
    if (!m->BSIM4rdsModGiven)   m->BSIM4rdsMod = 0;
    if (!m->BSIM4rbodyModGiven) m->BSIM4rbodyMod = 0;
    if (!m->BSIM4rgateModGiven) m->BSIM4rgateMod = 0;
    if (!m->BSIM4perModGiven)   m->BSIM4perMod = 1;
    if (!m->BSIM4geoModGiven)   m->BSIM4geoMod = 0;
    if (!m->BSIM4trnqsModGiven) m->BSIM4trnqsMod = 0;
    if (!m->BSIM4acnqsModGiven) m->BSIM4acnqsMod = 0;
    if (!m->BSIM4igcModGiven)   m->BSIM4igcMod = 0;
    if (!m->BSIM4igbModGiven)   m->BSIM4igbMod = 0;
    if (!m->BSIM4tempModGiven)  m->BSIM4tempMod = 0;

    // --- Model parameter defaults (subset read by BSIM4load) ---
    // (These Given flags may not exist; set unconditionally since
    // buildModelCard doesn't set them.)
    m->BSIM4vtm = Vt;
    m->BSIM4vtm0 = Vt;
    if (!m->BSIM4tnomGiven)     m->BSIM4tnom = Tnom;

    if (!m->BSIM4xjbvsGiven)    m->BSIM4xjbvs = 1.0;
    if (!m->BSIM4xjbvdGiven)    m->BSIM4xjbvd = m->BSIM4xjbvs;
    if (!m->BSIM4bvsGiven)      m->BSIM4bvs = 10.0;
    if (!m->BSIM4bvdGiven)      m->BSIM4bvd = m->BSIM4bvs;
    m->BSIM4vcrit = 0.5;
    if (!m->BSIM4gbminGiven)    m->BSIM4gbmin = 1.0e-12;

    m->BSIM4factor1 = 1.0;

    // Junction temperature saturation current densities
    if (!m->BSIM4SjctSatCurDensityGiven) m->BSIM4SjctSatCurDensity = 1.0e-4;
    if (!m->BSIM4DjctSatCurDensityGiven) m->BSIM4DjctSatCurDensity = m->BSIM4SjctSatCurDensity;
    if (!m->BSIM4SjctSidewallSatCurDensityGiven) m->BSIM4SjctSidewallSatCurDensity = 0.0;
    if (!m->BSIM4DjctSidewallSatCurDensityGiven) m->BSIM4DjctSidewallSatCurDensity = m->BSIM4SjctSidewallSatCurDensity;
    if (!m->BSIM4SjctGateSidewallSatCurDensityGiven) m->BSIM4SjctGateSidewallSatCurDensity = 0.0;
    if (!m->BSIM4DjctGateSidewallSatCurDensityGiven) m->BSIM4DjctGateSidewallSatCurDensity = m->BSIM4SjctGateSidewallSatCurDensity;
    if (!m->BSIM4SjctEmissionCoeffGiven) m->BSIM4SjctEmissionCoeff = 1.0;
    if (!m->BSIM4DjctEmissionCoeffGiven) m->BSIM4DjctEmissionCoeff = m->BSIM4SjctEmissionCoeff;

    // --- Additional model parameters accessed by BSIM4load ---
    // Junction potentials
    m->BSIM4PhiBS = 1.0;
    m->BSIM4PhiBSWS = 1.0;
    m->BSIM4PhiBSWGS = 1.0;
    m->BSIM4PhiBD = 1.0;
    m->BSIM4PhiBSWD = 1.0;
    m->BSIM4PhiBSWGD = 1.0;

    // Junction capacitances
    if (!m->BSIM4SunitAreaJctCapGiven) m->BSIM4SunitAreaJctCap = 5.0e-4;
    if (!m->BSIM4DunitAreaJctCapGiven) m->BSIM4DunitAreaJctCap = m->BSIM4SunitAreaJctCap;
    m->BSIM4SunitLengthSidewallJctCap = 5.0e-10;
    m->BSIM4DunitLengthSidewallJctCap = m->BSIM4SunitLengthSidewallJctCap;
    m->BSIM4SunitLengthGateSidewallJctCap = m->BSIM4SunitLengthSidewallJctCap;
    m->BSIM4DunitLengthGateSidewallJctCap = m->BSIM4SunitLengthGateSidewallJctCap;

    // Temperature-dependent junction capacitance
    m->BSIM4SunitAreaTempJctCap = m->BSIM4SunitAreaJctCap;
    m->BSIM4DunitAreaTempJctCap = m->BSIM4SunitAreaTempJctCap;
    m->BSIM4SunitLengthSidewallTempJctCap = m->BSIM4SunitLengthSidewallJctCap;
    m->BSIM4DunitLengthSidewallTempJctCap = m->BSIM4SunitLengthSidewallTempJctCap;
    m->BSIM4SunitLengthGateSidewallTempJctCap = m->BSIM4SunitLengthGateSidewallJctCap;
    m->BSIM4DunitLengthGateSidewallTempJctCap = m->BSIM4SunitLengthGateSidewallTempJctCap;

    // GIDL clamp
    if (!m->BSIM4gidlclampGiven) m->BSIM4gidlclamp = -1e-5;
    if (!m->BSIM4idovvdscGiven)  m->BSIM4idovvdsc = 1e-9;

    // DIBL / bodysource
    if (!m->BSIM4bdosGiven) m->BSIM4bdos = 1.0;
    if (!m->BSIM4adosGiven) m->BSIM4ados = 1.0;

    // Junction breakdown voltages
    if (!m->BSIM4vtsdGiven) m->BSIM4vtsd = 10.0;
    if (!m->BSIM4vtssGiven) m->BSIM4vtss = 10.0;
    if (!m->BSIM4vtsswsGiven) m->BSIM4vtssws = 10.0;
    if (!m->BSIM4vtsswdGiven) m->BSIM4vtsswd = 10.0;
    if (!m->BSIM4vtsswgsGiven) m->BSIM4vtsswgs = 10.0;
    if (!m->BSIM4vtsswgdGiven) m->BSIM4vtsswgd = 10.0;

    // Bandgap and eot
    m->BSIM4Eg0 = 1.16;
    if (!m->BSIM4eotGiven) m->BSIM4eot = 15.0e-10;

    // Junction grading coefficients
    m->BSIM4SbulkJctBotGradingCoeff = 0.5;
    m->BSIM4DbulkJctBotGradingCoeff = m->BSIM4SbulkJctBotGradingCoeff;
    m->BSIM4SbulkJctSideGradingCoeff = 0.33;
    m->BSIM4DbulkJctSideGradingCoeff = m->BSIM4SbulkJctSideGradingCoeff;
    m->BSIM4SbulkJctGateSideGradingCoeff = m->BSIM4SbulkJctSideGradingCoeff;
    m->BSIM4DbulkJctGateSideGradingCoeff = m->BSIM4SbulkJctGateSideGradingCoeff;

    // Temperature-dependent emission coefficients
    m->BSIM4njtsdtemp = 20.0;
    m->BSIM4njtsstemp = 20.0;
    m->BSIM4njtsswdtemp = 20.0;
    m->BSIM4njtsswgdtemp = 20.0;
    m->BSIM4njtsswgstemp = 20.0;
    m->BSIM4njtsswstemp = 20.0;

    // Partition coefficient
    if (!m->BSIM4xpartGiven) m->BSIM4xpart = 0.0;

    // --- Physical / oxide parameters ---
    if (!m->BSIM4epsroxGiven) m->BSIM4epsrox = 3.9;
    if (!m->BSIM4toxpGiven)   m->BSIM4toxp = m->BSIM4toxe;
    if (!m->BSIM4toxeGiven)   m->BSIM4toxe = 30.0e-10;

    // Cox = eps_ox / tox
    m->BSIM4coxe = m->BSIM4epsrox * eps0 / m->BSIM4toxp;

    // ------------------------------------------------------------------
    // 3. Create and populate pParam (bsim4SizeDependParam)
    // ------------------------------------------------------------------
    pParam_ = new bsim4SizeDependParam();
    std::memset(pParam_, 0, sizeof(bsim4SizeDependParam));
    inst->pParam = pParam_;

    // Effective channel dimensions (dl = dw = 0 for our subset)
    double Leff = L_;   // L - 2*dl
    double Weff = W_;   // W - 2*dw

    pParam_->Width = Weff;
    pParam_->Length = Leff;
    pParam_->NFinger = 1.0;

    // Copy model parameters into pParam (BSIM4load reads from pParam)
    pParam_->BSIM4cdsc   = m->BSIM4cdsc;
    pParam_->BSIM4nfactor= m->BSIM4nfactor;
    pParam_->BSIM4vsat   = m->BSIM4vsat;
    pParam_->BSIM4u0     = m->BSIM4u0;
    pParam_->BSIM4vth0   = m->BSIM4vth0;
    pParam_->BSIM4k1     = m->BSIM4k1;
    pParam_->BSIM4k2     = m->BSIM4k2;
    pParam_->BSIM4k3     = m->BSIM4k3;
    pParam_->BSIM4k3b    = m->BSIM4k3b;
    pParam_->BSIM4eta0   = m->BSIM4eta0;
    pParam_->BSIM4nsub   = m->BSIM4nsub;
    pParam_->BSIM4gamma1 = m->BSIM4gamma1;
    pParam_->BSIM4gamma2 = m->BSIM4gamma2;
    pParam_->BSIM4vbx    = m->BSIM4vbx;
    pParam_->BSIM4vbm    = m->BSIM4vbm;
    pParam_->BSIM4xt     = m->BSIM4xt;
    pParam_->BSIM4a0     = m->BSIM4a0;
    pParam_->BSIM4ags    = m->BSIM4ags;
    pParam_->BSIM4a1     = m->BSIM4a1;
    pParam_->BSIM4a2     = m->BSIM4a2;
    pParam_->BSIM4keta   = m->BSIM4keta;
    pParam_->BSIM4ketac  = m->BSIM4ketac;
    pParam_->BSIM4dvt0   = m->BSIM4dvt0;
    pParam_->BSIM4dvt1   = m->BSIM4dvt1;
    pParam_->BSIM4dvt2   = m->BSIM4dvt2;
    pParam_->BSIM4dvt0w  = m->BSIM4dvt0w;
    pParam_->BSIM4dvt1w  = m->BSIM4dvt1w;
    pParam_->BSIM4dvt2w  = m->BSIM4dvt2w;
    pParam_->BSIM4dsub   = m->BSIM4dsub;
    pParam_->BSIM4drout  = m->BSIM4drout;
    pParam_->BSIM4pclm   = m->BSIM4pclm;
    pParam_->BSIM4pdibl1 = m->BSIM4pdibl1;
    pParam_->BSIM4pdibl2 = m->BSIM4pdibl2;
    pParam_->BSIM4pdiblb = m->BSIM4pdiblb;
    pParam_->BSIM4ua     = m->BSIM4ua;
    pParam_->BSIM4ua1    = m->BSIM4ua1;
    pParam_->BSIM4ub     = m->BSIM4ub;
    pParam_->BSIM4ub1    = m->BSIM4ub1;
    pParam_->BSIM4uc     = m->BSIM4uc;
    pParam_->BSIM4uc1    = m->BSIM4uc1;
    pParam_->BSIM4ud     = m->BSIM4ud;
    pParam_->BSIM4ud1    = m->BSIM4ud1;
    pParam_->BSIM4up     = m->BSIM4up;
    pParam_->BSIM4lp     = m->BSIM4lp;
    pParam_->BSIM4eu     = m->BSIM4eu;
    pParam_->BSIM4ucs    = m->BSIM4ucs;
    pParam_->BSIM4ute    = m->BSIM4ute;
    pParam_->BSIM4ucste  = m->BSIM4ucste;
    pParam_->BSIM4voff   = m->BSIM4voff;
    pParam_->BSIM4delta  = m->BSIM4delta;
    pParam_->BSIM4minv   = m->BSIM4minv;
    pParam_->BSIM4minvcv = m->BSIM4minvcv;
    pParam_->BSIM4vfb    = m->BSIM4vfb;
    pParam_->BSIM4rdsw   = m->BSIM4rdsw;
    pParam_->BSIM4rsw    = m->BSIM4rsw;
    pParam_->BSIM4rdw    = m->BSIM4rdw;
    pParam_->BSIM4prwg   = m->BSIM4prwg;
    pParam_->BSIM4prwb   = m->BSIM4prwb;
    pParam_->BSIM4prt    = m->BSIM4prt;
    pParam_->BSIM4etab   = m->BSIM4etab;
    pParam_->BSIM4fprout = m->BSIM4fprout;
    pParam_->BSIM4pdits  = m->BSIM4pdits;
    pParam_->BSIM4pditsd = m->BSIM4pditsd;
    pParam_->BSIM4pscbe1 = m->BSIM4pscbe1;
    pParam_->BSIM4pscbe2 = m->BSIM4pscbe2;
    pParam_->BSIM4pvag   = m->BSIM4pvag;
    pParam_->BSIM4wr     = m->BSIM4wr;
    pParam_->BSIM4dwg    = m->BSIM4dwg;
    pParam_->BSIM4dwb    = m->BSIM4dwb;
    pParam_->BSIM4b0     = m->BSIM4b0;
    pParam_->BSIM4b1     = m->BSIM4b1;
    pParam_->BSIM4alpha0 = m->BSIM4alpha0;
    pParam_->BSIM4alpha1 = m->BSIM4alpha1;
    pParam_->BSIM4beta0  = m->BSIM4beta0;
    pParam_->BSIM4agidl  = m->BSIM4agidl;
    pParam_->BSIM4bgidl  = m->BSIM4bgidl;
    pParam_->BSIM4cgidl  = m->BSIM4cgidl;
    pParam_->BSIM4egidl  = m->BSIM4egidl;
    pParam_->BSIM4fgidl  = m->BSIM4fgidl;
    pParam_->BSIM4kgidl  = m->BSIM4kgidl;
    pParam_->BSIM4rgidl  = m->BSIM4rgidl;
    pParam_->BSIM4agisl  = m->BSIM4agisl;
    pParam_->BSIM4bgisl  = m->BSIM4bgisl;
    pParam_->BSIM4cgisl  = m->BSIM4cgisl;
    pParam_->BSIM4egisl  = m->BSIM4egisl;
    pParam_->BSIM4fgisl  = m->BSIM4fgisl;
    pParam_->BSIM4kgisl  = m->BSIM4kgisl;
    pParam_->BSIM4rgisl  = m->BSIM4rgisl;
    pParam_->BSIM4aigc   = m->BSIM4aigc;
    pParam_->BSIM4bigc   = m->BSIM4bigc;
    pParam_->BSIM4cigc   = m->BSIM4cigc;
    pParam_->BSIM4aigsd  = m->BSIM4aigsd;
    pParam_->BSIM4bigsd  = m->BSIM4bigsd;
    pParam_->BSIM4cigsd  = m->BSIM4cigsd;
    pParam_->BSIM4aigs   = m->BSIM4aigs;
    pParam_->BSIM4bigs   = m->BSIM4bigs;
    pParam_->BSIM4cigs   = m->BSIM4cigs;
    pParam_->BSIM4aigd   = m->BSIM4aigd;
    pParam_->BSIM4bigd   = m->BSIM4bigd;
    pParam_->BSIM4cigd   = m->BSIM4cigd;
    pParam_->BSIM4aigbacc= m->BSIM4aigbacc;
    pParam_->BSIM4bigbacc= m->BSIM4bigbacc;
    pParam_->BSIM4cigbacc= m->BSIM4cigbacc;
    pParam_->BSIM4aigbinv= m->BSIM4aigbinv;
    pParam_->BSIM4bigbinv= m->BSIM4bigbinv;
    pParam_->BSIM4cigbinv= m->BSIM4cigbinv;
    pParam_->BSIM4nigc   = m->BSIM4nigc;
    pParam_->BSIM4nigbinv= m->BSIM4nigbinv;
    pParam_->BSIM4nigbacc= m->BSIM4nigbacc;
    pParam_->BSIM4ntox   = m->BSIM4ntox;
    pParam_->BSIM4eigbinv= m->BSIM4eigbinv;
    pParam_->BSIM4pigcd  = m->BSIM4pigcd;
    pParam_->BSIM4poxedge= m->BSIM4poxedge;
    pParam_->BSIM4xrcrg1 = m->BSIM4xrcrg1;
    pParam_->BSIM4xrcrg2 = m->BSIM4xrcrg2;
    pParam_->BSIM4lambda = m->BSIM4lambda;
    pParam_->BSIM4vtl    = m->BSIM4vtl;
    pParam_->BSIM4xn     = m->BSIM4xn;
    pParam_->BSIM4lc     = m->BSIM4lc;
    pParam_->BSIM4vfbsdoff = m->BSIM4vfbsdoff;
    pParam_->BSIM4tvfbsdoff = m->BSIM4tvfbsdoff;
    pParam_->BSIM4ku0    = m->BSIM4ku0;
    pParam_->BSIM4kvth0  = m->BSIM4kvth0;

    pParam_->BSIM4kvth0we= m->BSIM4kvth0we;
    pParam_->BSIM4k2we   = m->BSIM4k2we;
    pParam_->BSIM4ku0we  = m->BSIM4ku0we;
    pParam_->BSIM4cgsl   = m->BSIM4cgsl;
    pParam_->BSIM4cgdl   = m->BSIM4cgdl;
    pParam_->BSIM4ckappas = m->BSIM4ckappas;
    pParam_->BSIM4ckappad = m->BSIM4ckappad;
    pParam_->BSIM4cf     = m->BSIM4cf;
    pParam_->BSIM4clc    = m->BSIM4clc;
    pParam_->BSIM4cle    = m->BSIM4cle;
    pParam_->BSIM4vfbcv  = m->BSIM4vfbcv;
    pParam_->BSIM4noff   = m->BSIM4noff;
    pParam_->BSIM4voffcv = m->BSIM4voffcv;
    pParam_->BSIM4acde   = m->BSIM4acde;
    pParam_->BSIM4moin   = m->BSIM4moin;

    // Pre-calculated constants (computed from model parameters)
    // phi: surface potential ≈ 2 * phi_F where phi_F = Vt * ln(Nsub/ni)
    double ni = 1.45e10;  // intrinsic carrier concentration (1/cm^3)
    double phiF = Vt * std::log(m->BSIM4nsub / ni);
    double phi = 2.0 * phiF;
    if (phi <= 0.0) phi = 0.7;  // fallback

    pParam_->BSIM4phi = phi;
    pParam_->BSIM4sqrtPhi = std::sqrt(phi);
    pParam_->BSIM4phis3 = phi * std::sqrt(phi);

    // Xdep0: depletion width at zero bias
    double Xdep0 = std::sqrt(2.0 * epsSi * phi / (q * m->BSIM4nsub * 1e6));  // nsub in 1/cm^3 -> 1/m^3
    pParam_->BSIM4Xdep0 = Xdep0;
    pParam_->BSIM4sqrtXdep0 = std::sqrt(Xdep0);

    // leff, weff
    pParam_->BSIM4leff = Leff;
    pParam_->BSIM4weff = Weff;
    pParam_->BSIM4leffCV = Leff;
    pParam_->BSIM4weffCV = Weff;
    pParam_->BSIM4weffCJ = Weff;

    // k1: body effect coefficient
    double Cox = m->BSIM4coxe;
    double k1 = std::sqrt(2.0 * q * epsSi * m->BSIM4nsub * 1e6) / Cox;
    pParam_->BSIM4k1 = k1;
    pParam_->BSIM4k1ox = k1 * k1;

    // k2: body effect (from model, default 0.05)
    pParam_->BSIM4k2ox = m->BSIM4k2;

    // vbi: built-in potential
    pParam_->BSIM4vbi = phi;

    // vth0 in pParam
    pParam_->BSIM4vth0 = m->BSIM4vth0;

    // lpe0: narrow width effect parameter
    pParam_->BSIM4lpe0 = 1.74e-7;  // default from BSIM4

    // dvt0, dvt1, dvt2 (from model)
    pParam_->BSIM4dvt0 = m->BSIM4dvt0;
    pParam_->BSIM4dvt1 = m->BSIM4dvt1;
    pParam_->BSIM4dvt2 = m->BSIM4dvt2;
    pParam_->BSIM4dvt0w = m->BSIM4dvt0w;
    pParam_->BSIM4dvt1w = m->BSIM4dvt1w;
    pParam_->BSIM4dvt2w = m->BSIM4dvt2w;

    // Temperature-dependent mobility and saturation velocity
    pParam_->BSIM4u0temp = m->BSIM4u0;
    pParam_->BSIM4vsattemp = m->BSIM4vsat;

    // Cofactors for capacitance
    pParam_->BSIM4cof1 = 0.0;
    pParam_->BSIM4cof2 = 0.0;
    pParam_->BSIM4cof3 = 0.0;
    pParam_->BSIM4cof4 = 0.0;

    // cdep0: zero-bias depletion capacitance
    pParam_->BSIM4cdep0 = std::sqrt(q * epsSi * m->BSIM4nsub * 1e6 / (2.0 * phi));

    // ToxRatio
    pParam_->BSIM4ToxRatio = m->BSIM4toxe / m->BSIM4toxp;

    // GIDL/GISL defaults
    if (pParam_->BSIM4bgidl <= 0.0) pParam_->BSIM4bgidl = 2.3e9;
    if (pParam_->BSIM4cgidl <= 0.0) pParam_->BSIM4cgidl = 0.5;
    if (pParam_->BSIM4egidl <= 0.0) pParam_->BSIM4egidl = 0.8;
    if (pParam_->BSIM4fgidl <= 0.0) pParam_->BSIM4fgidl = 1.0;
    if (pParam_->BSIM4bgisl <= 0.0) pParam_->BSIM4bgisl = pParam_->BSIM4bgidl;
    if (pParam_->BSIM4cgisl <= 0.0) pParam_->BSIM4cgisl = pParam_->BSIM4cgidl;
    if (pParam_->BSIM4egisl <= 0.0) pParam_->BSIM4egisl = pParam_->BSIM4egidl;
    if (pParam_->BSIM4fgisl <= 0.0) pParam_->BSIM4fgisl = pParam_->BSIM4fgidl;

    // Gate current defaults
    if (pParam_->BSIM4aigc <= 0.0) pParam_->BSIM4aigc = (m->BSIM4type == NMOS) ? 1.36e-2 : 9.80e-3;
    if (pParam_->BSIM4bigc <= 0.0) pParam_->BSIM4bigc = (m->BSIM4type == NMOS) ? 1.71e-3 : 7.59e-4;
    if (pParam_->BSIM4cigc <= 0.0) pParam_->BSIM4cigc = (m->BSIM4type == NMOS) ? 0.075 : 0.03;
    if (pParam_->BSIM4aigs <= 0.0) pParam_->BSIM4aigs = pParam_->BSIM4aigc;
    if (pParam_->BSIM4bigs <= 0.0) pParam_->BSIM4bigs = pParam_->BSIM4bigc;
    if (pParam_->BSIM4cigs <= 0.0) pParam_->BSIM4cigs = pParam_->BSIM4cigc;
    if (pParam_->BSIM4aigd <= 0.0) pParam_->BSIM4aigd = pParam_->BSIM4aigc;
    if (pParam_->BSIM4bigd <= 0.0) pParam_->BSIM4bigd = pParam_->BSIM4bigc;
    if (pParam_->BSIM4cigd <= 0.0) pParam_->BSIM4cigd = pParam_->BSIM4cigc;

    pParam_->pNext = nullptr;

    // ------------------------------------------------------------------
    // 4. Instance-level fields read by BSIM4load
    // ------------------------------------------------------------------

    // Junction area/perimeter (used for diode current)
    inst->BSIM4Aseff = inst->BSIM4sourceArea;   // source area
    inst->BSIM4Pseff = inst->BSIM4sourcePerimeter;  // source perimeter
    inst->BSIM4Adeff = inst->BSIM4drainArea;    // drain area
    inst->BSIM4Pdeff = inst->BSIM4drainPerimeter;  // drain perimeter

    // Threshold voltage (instance-level, used for init conditions)
    inst->BSIM4vth0 = m->BSIM4vth0;

    // Body effect coefficient (instance)
    inst->BSIM4vbsc = 0.0;  // body effect threshold (0 = normal operation)

    // von: turn-on voltage (used for limiting)
    inst->BSIM4von = m->BSIM4vth0 + 0.1;

    // Diode junction parameters (for bulk diode model)
    // Forward bias
    inst->BSIM4vjsmFwd = 0.8;   // junction forward voltage (smoothed)
    inst->BSIM4IVjsmFwd = 1.0e-3;  // current at vjsmFwd
    inst->BSIM4XExpBVS = 0.0;   // exponential term
    // Reverse bias
    inst->BSIM4vjsmRev = -2.0;  // junction reverse voltage
    inst->BSIM4IVjsmRev = 1.0e-12;
    inst->BSIM4SslpRev = 0.0;
    // Drain diode forward bias
    inst->BSIM4vjdmFwd = 0.8;
    inst->BSIM4IVjdmFwd = 1.0e-3;
    inst->BSIM4XExpBVD = 0.0;
    // Drain diode reverse bias
    inst->BSIM4vjdmRev = -2.0;
    inst->BSIM4IVjdmRev = 1.0e-12;
    inst->BSIM4DslpRev = 0.0;
    // Forward slope
    inst->BSIM4SslpFwd = 0.0;
    inst->BSIM4DslpFwd = 0.0;

    // Temperature-dependent reverse saturation currents
    inst->BSIM4SjctTempRevSatCur = 0.0;
    inst->BSIM4DjctTempRevSatCur = 0.0;
    inst->BSIM4SswTempRevSatCur = 0.0;
    inst->BSIM4DswTempRevSatCur = 0.0;
    inst->BSIM4SswgTempRevSatCur = 0.0;
    inst->BSIM4DswgTempRevSatCur = 0.0;

    // Body resistance conductances (instance-level)
    inst->BSIM4grbdb = 0.0;
    inst->BSIM4grbsb = 0.0;
    inst->BSIM4grbpb = 0.0;
    inst->BSIM4grbpd = 0.0;
    inst->BSIM4grbps = 0.0;

    // Operating mode (set by BSIM4load based on Vds sign)
    inst->BSIM4mode = 0;

    // ------------------------------------------------------------------
    // 5. State variable initialization
    // ------------------------------------------------------------------
    // All state variables start at 0 (will be set per-evaluation)
    // bsim4_state_ is already zeroed in constructor

    // ------------------------------------------------------------------
    // 6. Link model -> instance
    // ------------------------------------------------------------------
    m->BSIM4instances = inst;
    inst->BSIM4modPtr = m;
}

// =========================================================================
//  Set node matrix rows (called after topology is assigned)
// =========================================================================

void BSim4Adapter::setNodeRows(int rd, int rg, int rs, int rb) {
    dNode_ = rd; gNode_ = rg; sNode_ = rs; bNode_ = rb;

    // BSIM4load uses these directly as indices into CKTrhs/CKTmatrix.
    // For the level-2 subset (no rbody/rgate), prime nodes = main nodes.
    bsim4_inst_->BSIM4dNode = rd;
    bsim4_inst_->BSIM4dNodePrime = rd;
    bsim4_inst_->BSIM4gNodeExt = rg;
    bsim4_inst_->BSIM4gNodePrime = rg;
    bsim4_inst_->BSIM4gNodeMid = rg;
    bsim4_inst_->BSIM4sNode = rs;
    bsim4_inst_->BSIM4sNodePrime = rs;
    bsim4_inst_->BSIM4bNode = rb;
    bsim4_inst_->BSIM4bNodePrime = rb;
}

// =========================================================================
//  SpiceDeviceAdapter pure virtual implementations
// =========================================================================

void BSim4Adapter::discoverPattern(std::vector<int>& rows,
                                    std::vector<int>& cols) const {
    // Declare the 4x4 Jacobian pattern on terminal nodes (d, g, s, b).
    // Ground references use -1 (skipped by emit()).
    int term_rows[4] = {
        rowOf(nodes_[0]),  // drain
        rowOf(nodes_[1]),  // gate
        rowOf(nodes_[2]),  // source
        rowOf(nodes_[3]),  // bulk
    };

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            rows.push_back(term_rows[i]);
            cols.push_back(term_rows[j]);
        }
    }
}

void BSim4Adapter::evaluate(const CKTcircuit& /*ctx*/, void* /*inst*/) const {
    // Not used: BSim4Adapter uses its own evaluate(const DaeSolution&, double)
    // which calls BSIM4load() directly. This SpiceDeviceAdapter interface
    // method is kept for ABI compatibility but is never invoked.
}

// =========================================================================
//  Connect BSIM4*Ptr matrix-element pointers to g_cache_ slots
// =========================================================================

void BSim4Adapter::setJacobianSlots(const std::vector<int>& slots) {
    SpiceDeviceAdapter::setJacobianSlots(slots);

    // Initialize the SMPmatrix slot_of table so that the SMPaddElt shim
    // can translate (row, col) stamps into DaeMatrix slot indices.
    int size = static_cast<int>(nodes_.size());
    smp_.slot_of.assign(static_cast<size_t>(size) * static_cast<size_t>(size), -1);

    // The Jacobian pattern is a full 4x4 stamp on terminal nodes, declared
    // in row-major order by registerJacobianPattern(). The k-th entry maps
    // to (term_rows[k/4], term_rows[k%4]) with slot = slots[k].
    int term_rows[4] = {
        rowOf(nodes_[0]),  // drain
        rowOf(nodes_[1]),  // gate
        rowOf(nodes_[2]),  // source
        rowOf(nodes_[3]),  // bulk
    };
    for (size_t k = 0; k < slots.size(); ++k) {
        int r = static_cast<int>(k / 4);
        int c = static_cast<int>(k % 4);
        if (term_rows[r] >= 0 && term_rows[c] >= 0) {
            smp_.slot_of[static_cast<size_t>(term_rows[r]) * static_cast<size_t>(size)
                         + static_cast<size_t>(term_rows[c])] = slots[k];
        }
    }

    // The Berkeley BSIM4 code uses pre-computed matrix-element pointers
    // (BSIM4*Ptr fields) that BSIM4setup normally fills via SMPmakeElt().
    // Since we don't call BSIM4setup, we set them manually to point into
    // our g_cache_ buffer at the slot corresponding to each (row,col).
    //
    // Slot layout (from registerJacobianPattern, row-major 4x4):
    //   slot = row*4 + col, where row,col in {d=0,g=1,s=2,b=3}
    //
    // Terminal index mapping for the level-2 subset:
    //   D,DP,DB → 0 (drain), G,GP,GM,GE → 1 (gate),
    //   S,SP,SB → 2 (source), B,BP → 3 (bulk), Q → -1 (skip)
    //
    // A pointer BSIM4<ROW><COL>Ptr maps to slot = t(ROW)*4 + t(COL).

    auto slotOf = [this](int row, int col) -> double* {
        if (row < 0 || col < 0) return nullptr;
        int slot = row * 4 + col;
        if (slot >= static_cast<int>(g_cache_.size())) return nullptr;
        return &g_cache_[static_cast<size_t>(slot)];
    };

    auto T = [](const char* prefix) -> int {
        if (prefix[0] == 'D' || prefix[0] == 'd') return 0;   // D, DP, DB
        if (prefix[0] == 'G' || prefix[0] == 'g') return 1;   // G, GP, GM, GE
        if (prefix[0] == 'S' || prefix[0] == 's') return 2;   // S, SP, SB
        if (prefix[0] == 'B' || prefix[0] == 'b') return 3;   // B, BP
        return -1;                                            // Q or unknown
    };

    // Helper: set a BSIM4*Ptr field from its row/col prefix strings.
    auto setPtr = [&](double* BSIM4instance::* field,
                      const char* rowPref, const char* colPref) {
        int r = T(rowPref), c = T(colPref);
        bsim4_inst_->*field = slotOf(r, c);
    };

    // --- D-row pointers ---
    setPtr(&BSIM4instance::BSIM4DPdPtr,   "DP", "d");
    setPtr(&BSIM4instance::BSIM4DPdpPtr,  "DP", "dp");
    setPtr(&BSIM4instance::BSIM4DPgpPtr,  "DP", "gp");
    setPtr(&BSIM4instance::BSIM4DPgmPtr,  "DP", "gm");
    setPtr(&BSIM4instance::BSIM4DPspPtr,  "DP", "sp");
    setPtr(&BSIM4instance::BSIM4DPbpPtr,  "DP", "bp");
    setPtr(&BSIM4instance::BSIM4DPdbPtr,  "DP", "db");
    setPtr(&BSIM4instance::BSIM4DdPtr,    "D",  "d");
    setPtr(&BSIM4instance::BSIM4DdpPtr,   "D",  "dp");

    // --- G-row pointers ---
    setPtr(&BSIM4instance::BSIM4GPdpPtr,  "GP", "dp");
    setPtr(&BSIM4instance::BSIM4GPgpPtr,  "GP", "gp");
    setPtr(&BSIM4instance::BSIM4GPgmPtr,  "GP", "gm");
    setPtr(&BSIM4instance::BSIM4GPgePtr,  "GP", "ge");
    setPtr(&BSIM4instance::BSIM4GPspPtr,  "GP", "sp");
    setPtr(&BSIM4instance::BSIM4GPbpPtr,  "GP", "bp");
    setPtr(&BSIM4instance::BSIM4GMdpPtr,  "GM", "dp");
    setPtr(&BSIM4instance::BSIM4GMgpPtr,  "GM", "gp");
    setPtr(&BSIM4instance::BSIM4GMgmPtr,  "GM", "gm");
    setPtr(&BSIM4instance::BSIM4GMgePtr,  "GM", "ge");
    setPtr(&BSIM4instance::BSIM4GMspPtr,  "GM", "sp");
    setPtr(&BSIM4instance::BSIM4GMbpPtr,  "GM", "bp");
    setPtr(&BSIM4instance::BSIM4GEdpPtr,  "GE", "dp");
    setPtr(&BSIM4instance::BSIM4GEgpPtr,  "GE", "gp");
    setPtr(&BSIM4instance::BSIM4GEgmPtr,  "GE", "gm");
    setPtr(&BSIM4instance::BSIM4GEgePtr,  "GE", "ge");
    setPtr(&BSIM4instance::BSIM4GEspPtr,  "GE", "sp");
    setPtr(&BSIM4instance::BSIM4GEbpPtr,  "GE", "bp");

    // --- S-row pointers ---
    setPtr(&BSIM4instance::BSIM4SPdpPtr,  "SP", "dp");
    setPtr(&BSIM4instance::BSIM4SPgpPtr,  "SP", "gp");
    setPtr(&BSIM4instance::BSIM4SPgmPtr,  "SP", "gm");
    setPtr(&BSIM4instance::BSIM4SPsPtr,   "SP", "s");
    setPtr(&BSIM4instance::BSIM4SPspPtr,  "SP", "sp");
    setPtr(&BSIM4instance::BSIM4SPbpPtr,  "SP", "bp");
    setPtr(&BSIM4instance::BSIM4SPsbPtr,  "SP", "sb");
    setPtr(&BSIM4instance::BSIM4SspPtr,   "S",  "sp");
    setPtr(&BSIM4instance::BSIM4SsPtr,    "S",  "s");

    // --- B-row pointers ---
    setPtr(&BSIM4instance::BSIM4BPdpPtr,  "BP", "dp");
    setPtr(&BSIM4instance::BSIM4BPgpPtr,  "BP", "gp");
    setPtr(&BSIM4instance::BSIM4BPgmPtr,  "BP", "gm");
    setPtr(&BSIM4instance::BSIM4BPspPtr,  "BP", "sp");
    setPtr(&BSIM4instance::BSIM4BPdbPtr,  "BP", "db");
    setPtr(&BSIM4instance::BSIM4BPbPtr,   "BP", "b");
    setPtr(&BSIM4instance::BSIM4BPsbPtr,  "BP", "sb");
    setPtr(&BSIM4instance::BSIM4BPbpPtr,  "BP", "bp");
    setPtr(&BSIM4instance::BSIM4DBdpPtr,  "DB", "dp");
    setPtr(&BSIM4instance::BSIM4DBdbPtr,  "DB", "db");
    setPtr(&BSIM4instance::BSIM4DBbpPtr,  "DB", "bp");
    setPtr(&BSIM4instance::BSIM4DBbPtr,   "DB", "b");
    setPtr(&BSIM4instance::BSIM4SBspPtr,  "SB", "sp");
    setPtr(&BSIM4instance::BSIM4SBbpPtr,  "SB", "bp");
    setPtr(&BSIM4instance::BSIM4SBbPtr,   "SB", "b");
    setPtr(&BSIM4instance::BSIM4SBsbPtr,  "SB", "sb");
    setPtr(&BSIM4instance::BSIM4BdbPtr,   "B",  "db");
    setPtr(&BSIM4instance::BSIM4BbpPtr,   "B",  "bp");
    setPtr(&BSIM4instance::BSIM4BsbPtr,   "B",  "sb");
    setPtr(&BSIM4instance::BSIM4BbPtr,    "B",  "b");

    // --- Cross D/S/B pointers ---
    setPtr(&BSIM4instance::BSIM4DgpPtr,   "D",  "gp");
    setPtr(&BSIM4instance::BSIM4DspPtr,   "D",  "sp");
    setPtr(&BSIM4instance::BSIM4DbpPtr,   "D",  "bp");
    setPtr(&BSIM4instance::BSIM4SdpPtr,   "S",  "dp");
    setPtr(&BSIM4instance::BSIM4SgpPtr,   "S",  "gp");
    setPtr(&BSIM4instance::BSIM4SbpPtr,   "S",  "bp");

    // --- Q-row pointers (charge node, not used in level-2 subset) ---
    bsim4_inst_->BSIM4QdpPtr = nullptr;
    bsim4_inst_->BSIM4QgpPtr = nullptr;
    bsim4_inst_->BSIM4QspPtr = nullptr;
    bsim4_inst_->BSIM4QbpPtr = nullptr;
    bsim4_inst_->BSIM4QqPtr  = nullptr;
    bsim4_inst_->BSIM4DPqPtr = nullptr;
    bsim4_inst_->BSIM4GPqPtr = nullptr;
    bsim4_inst_->BSIM4SPqPtr = nullptr;
}

// =========================================================================
//  Jacobian pattern
// =========================================================================

void BSim4Adapter::registerJacobianPattern(std::vector<int>& rows,
                                           std::vector<int>& cols) const {
    // Full 4x4 stamp on terminal nodes (d, g, s, b).
    // Ground references use -1 (skipped by emit()).
    int term_rows[4] = {
        rowOf(nodes_[0]),  // drain
        rowOf(nodes_[1]),  // gate
        rowOf(nodes_[2]),  // source
        rowOf(nodes_[3]),  // bulk
    };

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            rows.push_back(term_rows[i]);
            cols.push_back(term_rows[j]);
        }
    }
}

// =========================================================================
//  Evaluation: call BSIM4load() (Berkeley, unchanged)
// =========================================================================

void BSim4Adapter::evaluate(const DaeSolution& x, double time) const {
    (void)time;

    // --- Read terminal node voltages from solution ------------------------
    const std::vector<double>& xvec = x.data();
    const double* xdata = xvec.data();

    int rd = rowOf(nodes_[0]);  // drain
    int rg = rowOf(nodes_[1]);  // gate
    int rs = rowOf(nodes_[2]);  // source
    int rb = rowOf(nodes_[3]);  // bulk

    double vd = (rd >= 0) ? xdata[rd] : 0.0;
    double vg = (rg >= 0) ? xdata[rg] : 0.0;
    double vs = (rs >= 0) ? xdata[rs] : 0.0;
    double vb = (rb >= 0) ? xdata[rb] : 0.0;

    // Berkeley BSIM4 expects PHYSICAL terminal voltages (not N-oriented).
    // The model's BSIM4type field (NMOS/PMOS) handles polarity internally.
    // vgs = vg - vs, vds = vd - vs, vbs = vb - vs (physical frame)
    double vgs_phys = vg - vs;
    double vds_phys = vd - vs;
    double vbs_phys = vb - vs;
    double vbd_phys = vbs_phys - vds_phys;

    // Write terminal voltages into the private state buffer (CKTstate0).
    // BSIM4states = 0, so BSIM4vbd=0, BSIM4vbs=1, BSIM4vgs=2, BSIM4vds=3
    // index directly into bsim4_state_[0..3].
    double* s0 = bsim4_state_.data();
    s0[0] = vbd_phys;  // BSIM4vbd
    s0[1] = vbs_phys;  // BSIM4vbs
    s0[2] = vgs_phys;  // BSIM4vgs
    s0[3] = vds_phys;  // BSIM4vds

    // --- Build the shimmed CKTcircuit -------------------------------------
    ctx_.CKTstate0 = bsim4_state_.data();  // BSIM4 reads terminal voltages here
    ctx_.CKTmatrix = &smp_;
    // Initialize SMPmatrix size for CKTnumNodes
    smp_.size = static_cast<int>(nodes_.size());

    ctx_.CKTnumNodes = smp_.size;
    ctx_.CKTmode = 0x0001;  // MODEDCOP (DC operating point)
    ctx_.CKTstate1 = bsim4_state_.data();
    ctx_.CKTstate2 = bsim4_state_.data();
    ctx_.CKTstate3 = bsim4_state_.data();
    ctx_.CKTdelta = 1.0e-12;
    for (int i = 0; i < 8; ++i) ctx_.CKTdeltaOld[i] = 1.0e-12;

    // Point the SMPaddElt shim at the Jacobian cache buffer.
    // g_cache_ is sized to the number of slots; slot indices from
    // setJacobianSlots index into it directly.
    if (static_cast<int>(g_cache_.size()) < static_cast<int>(slots_.size())) {
        g_cache_.assign(slots_.size(), 0.0);
    }
    smp_.values = g_cache_.data();

    // Zero the cache before stamping (BSIM4load uses +=).
    std::fill(g_cache_.begin(), g_cache_.end(), 0.0);

    // Call BSIM4load (Berkeley, unchanged)
    // BSIM4load stamps:
    //   CKTrhs[row]   += ceq   (current contributions → RHS)
    //   SMPaddElt(...) += g    (conductances → g_cache_ via smp_.values)
    //
    // CKTrhs was set by setRhsPtr() (points to F.data() during loadF).
    // SMPaddElt writes to g_cache_[slot].
    //
    // BSIM4load takes GENmodel* — BSIM4model has a different layout, so
    // we reinterpret_cast the pointer (BSIM4load casts it back internally).
    int rc = -1;
    rc = BSIM4load(reinterpret_cast<::GENmodel*>(bsim4_model_), &ctx_);
    if (rc != OK) {
        // Non-fatal: continue with cached values.
        (void)rc;
    }

    // --- Read results from the instance ----------------------------------
    readResults();
    cache_valid_ = true;
    cached_x_ = &x;
}

// =========================================================================
//  Read results from the Berkeley instance (const)
// =========================================================================

void BSim4Adapter::readResults() const {
    auto* inst = bsim4_inst_;

    id_ = inst->BSIM4cd;
    is_ = -id_;
    ib_ = 0.0;

    gm_   = inst->BSIM4gm;
    gds_  = inst->BSIM4gds;
    gmbs_ = inst->BSIM4gmbs;

    vth_ = (inst->BSIM4vth0 != 0.0) ? inst->BSIM4vth0 : model_->VTH0;
}

// =========================================================================
//  DAE contract: loadF — residual
// =========================================================================

void BSim4Adapter::loadF(const DaeSolution& x, double time,
                         std::vector<double>& F) const {
    // Ensure CKTrhs and CKTrhsOld point to F so BSIM4load can read/write them.
    ctx_.CKTrhs = F.data();
    ctx_.CKTrhsOld = F.data();

    if (!cache_valid_ || cached_x_ != &x) {
        evaluate(x, time);
    }
    // After BSIM4load, F already contains the device's ceq contribution
    // (stamped directly via CKTrhs = F.data()). Nothing more to do.
}

// =========================================================================
//  DAE contract: loadQ — charges (no charge model in this subset)
// =========================================================================

void BSim4Adapter::loadQ(const DaeSolution& x, double time,
                         std::vector<double>& Q) const {
    (void)x; (void)time; (void)Q;
    // Level-2 subset: no charge storage (no Meyer/NQS/capMod).
}

// =========================================================================
//  DAE contract: loadB — independent sources (none for a MOSFET)
// =========================================================================

void BSim4Adapter::loadB(double time, std::vector<double>& B) const {
    (void)time; (void)B;
}

// =========================================================================
//  DAE contract: loadJacobianF — dF/dx
// =========================================================================

void BSim4Adapter::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                                 double time) const {
    (void)time;

    // Ensure evaluation is current (BSIM4load stamps conductances into
    // g_cache_, which we now emit to J via the pre-resolved slots).
    if (!cache_valid_ || cached_x_ != &x) {
        evaluate(x, time);
    }

    // Emit cached conductances to J using pre-resolved slots.
    // slots_[i] is the DaeMatrix NNZ position; g_cache_[i] is the value.
    for (size_t i = 0; i < slots_.size(); ++i) {
        int slot = slots_[i];
        if (slot < 0) continue;  // ground reference, skip
        double g = g_cache_[i];
        J.accumulateSlot(slot, g);
    }
}

// =========================================================================
//  DAE contract: loadJacobianQ — dQ/dx (no charges)
// =========================================================================

void BSim4Adapter::loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                                 double time) const {
    (void)J; (void)x; (void)time;
}

// =========================================================================
//  Temperature notification
// =========================================================================

void BSim4Adapter::setCircuitTemperature(double kelvin) {
    (void)kelvin;
    if (bsim4_model_) {
        BSIM4temp(reinterpret_cast<::GENmodel*>(bsim4_model_), &ctx_);
    }
}

} // namespace e_XSpice

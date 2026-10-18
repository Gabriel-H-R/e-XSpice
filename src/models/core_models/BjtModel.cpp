/**
 * @author Gabriel Henrique Silva
 * @file BjtModel.cpp
 * @brief Shared BJT model card — level 1 (Ebers-Moll) and level 2
 *        (Gummel-Poon) for NPN and PNP.
 */

#include "BjtModel.h"
#include "../core/DeviceConstants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace e_XSpice {

BjtModel::BjtModel(std::string name, bool pnp)
    : DeviceModel(std::move(name)), pnp_(pnp) {
    // Card parameter table (canonical first, aliases after — canonical
    // wins when a card spells both). "FC" here is the BASE-COLLECTOR
    // coefficient (member FC); the BE coefficient is "FE". (An older
    // revision mapped card FC→FE and left member FC unsettable; no
    // in-repo card used FC, so nothing changes.)
    declareParams({
        makeIntParam("LEVEL", &BjtModel::LEVEL, "1=Ebers-Moll 2=Gummel-Poon"),
        makeParam("IS", &BjtModel::IS, "transport saturation current"),
        makeParam("BF", &BjtModel::BF, "ideal forward beta"),
        makeParam("BR", &BjtModel::BR, "ideal reverse beta"),
        makeParam("NF", &BjtModel::NF, "forward emission coefficient"),
        makeParam("NR", &BjtModel::NR, "reverse emission coefficient"),
        makeParam("VAF", &BjtModel::VAF, "forward Early voltage"),
        makeAlias("VA", &BjtModel::VAF),
        makeParam("VAR", &BjtModel::VAR, "reverse Early voltage"),
        makeAlias("VB", &BjtModel::VAR),
        makeParam("TF", &BjtModel::TF, "forward transit time"),
        makeParam("TR", &BjtModel::TR, "reverse transit time"),
        makeParam("CJE", &BjtModel::CJE, "BE zero-bias capacitance"),
        makeAlias("CJ", &BjtModel::CJE),
        makeParam("VJE", &BjtModel::VJE, "BE built-in potential"),
        makeAlias("PE", &BjtModel::VJE),
        makeParam("MJE", &BjtModel::MJE, "BE grading coefficient"),
        makeAlias("ME", &BjtModel::MJE),
        makeParam("FE", &BjtModel::FE, "BE forward-bias coefficient"),
        makeParam("CJC", &BjtModel::CJC, "BC zero-bias capacitance"),
        makeParam("VJC", &BjtModel::VJC, "BC built-in potential"),
        makeAlias("PC", &BjtModel::VJC),
        makeParam("MJC", &BjtModel::MJC, "BC grading coefficient"),
        makeAlias("MC", &BjtModel::MJC),
        makeParam("FC", &BjtModel::FC, "BC forward-bias coefficient"),
        makeParam("IKF", &BjtModel::IKF, "forward high-injection knee"),
        makeAlias("IK", &BjtModel::IKF),
        makeParam("IKR", &BjtModel::IKR, "reverse high-injection knee"),
        makeParam("ISE", &BjtModel::ISE, "BE non-ideal saturation current"),
        makeAlias("C2", &BjtModel::ISE),
        makeParam("NE", &BjtModel::NE, "BE non-ideal emission coefficient"),
        makeParam("ISC", &BjtModel::ISC, "BC non-ideal saturation current"),
        makeAlias("C4", &BjtModel::ISC),
        makeParam("NC", &BjtModel::NC, "BC non-ideal emission coefficient"),
        makeParam("RB", &BjtModel::RB, "base ohmic resistance"),
        makeAlias("RBM", &BjtModel::RB),
        makeParam("RC", &BjtModel::RC, "collector ohmic resistance"),
        makeAlias("RCC", &BjtModel::RC),
        makeParam("RE", &BjtModel::RE, "emitter ohmic resistance"),
        makeAlias("REE", &BjtModel::RE),
        makeParam("ISS", &BjtModel::ISS, "substrate saturation current"),
        makeParam("CJS", &BjtModel::CJS, "C-S zero-bias capacitance"),
        makeParam("VJS", &BjtModel::VJS, "C-S built-in potential"),
        makeAlias("PS", &BjtModel::VJS),
        makeParam("MJS", &BjtModel::MJS, "C-S grading coefficient"),
        makeAlias("MS", &BjtModel::MJS),
        makeParam("PTF", &BjtModel::PTF, "excess phase (AC-only)"),
        makeParam("EG", &BjtModel::EG, "bandgap voltage"),
        makeParam("XTI", &BjtModel::XTI, "saturation-current temp exponent"),
        makeParam("GMIN", &BjtModel::GMIN, "parallel conductance"),
        makeCelsiusParam("TNOM", &BjtModel::TNOM, "reference temperature"),
        makeCelsiusParam("TEMP", &BjtModel::TEMP, "card temperature"),
    });
}

void BjtModel::applyGuards() {
    // History-independent guards (Xyce Model::processParams).
    MJE = std::min(MJE, 0.9);
    MJC = std::min(MJC, 0.9);
    FE  = std::min(FE, 0.95);
    FC  = std::min(FC, 0.95);
    EG  = std::max(EG, 0.1);
    if (MJS <= 0.0) MJS = 0.33;   // SPICE default when unspecified
    MJS = std::min(MJS, 0.9);
}

bool BjtModel::processParams(std::string& message) {
    if (LEVEL != 1 && LEVEL != 2) {
        message = "LEVEL must be 1 (Ebers-Moll) or 2 (Gummel-Poon)";
        return false;
    }
    if (IS  <= 0.0) { message = "IS must be positive";  return false; }
    if (BF  <= 0.0) { message = "BF must be positive";  return false; }
    if (BR  <= 0.0) { message = "BR must be positive";  return false; }
    if (NF  <= 0.0) { message = "NF must be positive";  return false; }
    if (NR  <= 0.0) { message = "NR must be positive";  return false; }
    if (NE  <= 0.0) { message = "NE must be positive";  return false; }
    if (NC  <= 0.0) { message = "NC must be positive";  return false; }
    if (CJE <  0.0) { message = "CJE must be non-negative"; return false; }
    if (CJC <  0.0) { message = "CJC must be non-negative"; return false; }
    if (CJS <  0.0) { message = "CJS must be non-negative"; return false; }
    if (VJE <= 0.0) { message = "VJE must be positive"; return false; }
    if (VJC <= 0.0) { message = "VJC must be positive"; return false; }
    if (VJS <= 0.0) { message = "VJS must be positive"; return false; }
    if (TF  <  0.0) { message = "TF must be non-negative"; return false; }
    if (TR  <  0.0) { message = "TR must be non-negative"; return false; }
    if (IKF <  0.0) { message = "IKF must be non-negative"; return false; }
    if (IKR <  0.0) { message = "IKR must be non-negative"; return false; }
    if (ISE <  0.0) { message = "ISE must be non-negative"; return false; }
    if (ISC <  0.0) { message = "ISC must be non-negative"; return false; }
    if (ISS <  0.0) { message = "ISS must be non-negative"; return false; }
    if (RB  <  0.0) { message = "RB must be non-negative"; return false; }
    if (RC  <  0.0) { message = "RC must be non-negative"; return false; }
    if (RE  <  0.0) { message = "RE must be non-negative"; return false; }
    if (TEMP <= 0.0 || TNOM <= 0.0) {
        message = "TEMP/TNOM must be positive (Kelvin)";
        return false;
    }
    if (LEVEL == 1 && (IKF > 0.0 || IKR > 0.0 || ISE > 0.0 || ISC > 0.0 ||
                       RB > 0.0 || RC > 0.0 || RE > 0.0 || ISS > 0.0 ||
                       CJS > 0.0)) {
        message = "LEVEL=1 card carries level-2 parameters "
                  "(IKF/IKR/ISE/ISC/RB/RC/RE/ISS/CJS need LEVEL=2)";
        return false;
    }

    applyGuards();
    updateTemperature(TEMP);
    return true;
}

void BjtModel::updateTemperature(double tempKelvin) {
    thermal_ = scaleTo(tempKelvin);
}

BjtThermal BjtModel::computeThermal(double tempKelvin) const {
    return scaleTo(tempKelvin);
}

BjtThermal BjtModel::scaleTo(double Temp) const {
    BjtThermal t;

    t.vt  = Const::KoverQ * Temp;
    t.vte = NF * t.vt;
    t.vtc = NR * t.vt;
    t.vteE = NE * t.vt;
    t.vtcC = NC * t.vt;

    // Saturation current per junction: the SPICE formula depends on the
    // emission coefficient of that junction (NF for BE, NR for BC).
    t.beIsat = scaleSatCurrent(IS, NF, EG, XTI, Temp, TNOM);
    t.bcIsat = scaleSatCurrent(IS, NR, EG, XTI, Temp, TNOM);
    t.betaf  = t.beIsat / BF;
    t.betar  = t.bcIsat / BR;

    // Level-2 non-ideal saturation currents scale with their own emission
    // coefficients (SPICE3: same IS(T) law, NE/NC in place of N).
    t.iseT = (ISE > 0.0) ? scaleSatCurrent(ISE, NE, EG, XTI, Temp, TNOM) : 0.0;
    t.iscT = (ISC > 0.0) ? scaleSatCurrent(ISC, NC, EG, XTI, Temp, TNOM) : 0.0;
    t.issT = (ISS > 0.0) ? scaleSatCurrent(ISS, 1.0, EG, XTI, Temp, TNOM) : 0.0;

    // High-injection knees are temperature-independent (see header).
    t.ikf = IKF;
    t.ikr = IKR;

    // Junction capacitance constants (shared temperature scaling).
    t.beCap = scaleJunctionCap(VJE, MJE, FE, CJE, Temp, TNOM);
    t.bcCap = scaleJunctionCap(VJC, MJC, FC, CJC, Temp, TNOM);
    t.csCap = scaleJunctionCap(VJS, MJS, FC, CJS, Temp, TNOM);

    // PNJLIM critical voltage per junction (SPICE3f5 value).
    t.beVcrit = t.vte * std::log(t.vte / (Const::Root2 * t.beIsat));
    t.bcVcrit = t.vtc * std::log(t.vtc / (Const::Root2 * t.bcIsat));

    t.tf   = TF;
    t.tr   = TR;
    t.gmin = GMIN;
    return t;
}

RCP<const BjtModel> BjtModel::defaultModel() {
    static RCP<BjtModel> instance;
    if (!instance) {
        instance = make_rcp<BjtModel>("$DEFAULT_NPN", false);
        std::string message;
        if (!instance->processParams(message)) {
            throw std::logic_error("default BJT model invalid: " + message);
        }
    }
    return instance;
}

} // namespace e_XSpice

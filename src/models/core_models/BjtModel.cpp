/**
 * @author Gabriel Henrique Silva
 * @file BjtModel.cpp
 * @brief Shared BJT model card — transport (Ebers-Moll) model, level 1.
 */

#include "BjtModel.h"
#include "../core/DeviceConstants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace e_XSpice {

BjtModel::BjtModel(std::string name, bool pnp)
    : DeviceModel(std::move(name)), pnp_(pnp) {}

void BjtModel::applyGuards() {
    // History-independent guards (Xyce Model::processParams).
    MJE = std::min(MJE, 0.9);
    MJC = std::min(MJC, 0.9);
    FE  = std::min(FE, 0.95);
    FC  = std::min(FC, 0.95);
    EG  = std::max(EG, 0.1);
}

bool BjtModel::processParams(std::string& message) {
    if (IS  <= 0.0) { message = "IS must be positive";  return false; }
    if (BF  <= 0.0) { message = "BF must be positive";  return false; }
    if (BR  <= 0.0) { message = "BR must be positive";  return false; }
    if (NF  <= 0.0) { message = "NF must be positive";  return false; }
    if (NR  <= 0.0) { message = "NR must be positive";  return false; }
    if (CJE <  0.0) { message = "CJE must be non-negative"; return false; }
    if (CJC <  0.0) { message = "CJC must be non-negative"; return false; }
    if (VJE <= 0.0) { message = "VJE must be positive"; return false; }
    if (VJC <= 0.0) { message = "VJC must be positive"; return false; }
    if (TF  <  0.0) { message = "TF must be non-negative"; return false; }
    if (TR  <  0.0) { message = "TR must be non-negative"; return false; }
    if (TEMP <= 0.0 || TNOM <= 0.0) {
        message = "TEMP/TNOM must be positive (Kelvin)";
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

    // Saturation current per junction: the SPICE formula depends on the
    // emission coefficient of that junction (NF for BE, NR for BC).
    t.beIsat = scaleSatCurrent(IS, NF, EG, XTI, Temp, TNOM);
    t.bcIsat = scaleSatCurrent(IS, NR, EG, XTI, Temp, TNOM);
    t.betaf  = t.beIsat / BF;
    t.betar  = t.bcIsat / BR;

    // Junction capacitance constants (shared temperature scaling).
    t.beCap = scaleJunctionCap(VJE, MJE, FE, CJE, Temp, TNOM);
    t.bcCap = scaleJunctionCap(VJC, MJC, FC, CJC, Temp, TNOM);

    // PNJLIM critical voltage per junction (SPICE3f5 value).
    t.beVcrit = t.vte * std::log(t.vte / (Const::Root2 * t.beIsat));
    t.bcVcrit = t.vtc * std::log(t.vtc / (Const::Root2 * t.bcIsat));

    t.tf   = TF;
    t.tr   = TR;
    t.gmin = GMIN;
    return t;
}

RCP<DeviceModel> BjtModel::fromCard(const std::string& name,
                                   const ParamMap& params, bool pnp) {
    ModelCard card;
    card.name = name;
    card.type = pnp ? "PNP" : "NPN";
    card.params = params;

    RCP<BjtModel> m = make_rcp<BjtModel>(name, pnp);
    m->IS  = card.get("IS", m->IS);
    m->BF  = card.get("BF", m->BF);
    m->BR  = card.get("BR", m->BR);
    m->NF  = card.get("NF", m->NF);
    m->NR  = card.get("NR", m->NR);
    m->VAF = card.get("VAF", card.get("VA", m->VAF));
    m->VAR = card.get("VAR", card.get("VB", m->VAR));
    m->TF  = card.get("TF", m->TF);
    m->TR  = card.get("TR", m->TR);

    m->CJE = card.get("CJE", card.get("CJ", m->CJE));
    m->VJE = card.get("VJE", card.get("PE", m->VJE));
    m->MJE = card.get("MJE", card.get("ME", m->MJE));
    m->FE  = card.get("FC", card.get("FE", m->FE));

    m->CJC = card.get("CJC", card.get("CJS", m->CJC));
    m->VJC = card.get("VJC", card.get("PC", m->VJC));
    m->MJC = card.get("MJC", card.get("MC", m->MJC));

    m->EG   = card.get("EG", m->EG);
    m->XTI  = card.get("XTI", m->XTI);
    m->GMIN = card.get("GMIN", m->GMIN);
    // SPICE writes TNOM/TEMP in Celsius; the card stores Kelvin.
    if (card.given("TNOM")) m->TNOM = card.get("TNOM", 27.0) + Const::CtoK;
    if (card.given("TEMP")) m->TEMP = card.get("TEMP", 27.0) + Const::CtoK;
    return m;
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

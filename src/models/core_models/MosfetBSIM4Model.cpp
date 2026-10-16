/**
 * @author Gabriel Henrique Silva
 * @file MosfetBSIM4Model.cpp
 * @brief BSIM4 4.8.3 model card — level 2 MOSFET (NMOS + PMOS).
 */

#include "MosfetBSIM4Model.h"
#include "../core/DeviceConstants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace e_XSpice {

MosfetBSIM4Model::MosfetBSIM4Model(std::string name, bool nmos)
    : DeviceModel(std::move(name)), nmos_(nmos) {
    // Card parameter table (canonical first, aliases after — canonical wins
    // when a card spells both). Adding a parameter = one line here.
    declareParams({
        makeIntParam("LEVEL", &MosfetBSIM4Model::LEVEL, "2=BSIM4"),
        makeParam("VTH0", &MosfetBSIM4Model::VTH0, "zero-bias threshold voltage"),
        makeParam("KP", &MosfetBSIM4Model::KP, "transconductance parameter"),
        makeParam("LAMBDA", &MosfetBSIM4Model::LAMBDA,
                  "channel-length modulation"),
        makeParam("PHI", &MosfetBSIM4Model::PHI, "surface potential"),
        makeParam("GAMMA", &MosfetBSIM4Model::GAMMA,
                  "bulk threshold parameter"),
        makeParam("CBD", &MosfetBSIM4Model::CBD,
                  "bulk-drain capacitance (total)"),
        makeParam("CBS", &MosfetBSIM4Model::CBS,
                  "bulk-source capacitance (total)"),
        makeParam("IS", &MosfetBSIM4Model::IS, "bulk saturation current"),
        makeParam("PB", &MosfetBSIM4Model::PB, "bulk junction potential"),
        makeParam("MJ", &MosfetBSIM4Model::MJ, "bulk junction grading"),
        makeParam("FC", &MosfetBSIM4Model::FC,
                  "bulk forward-bias coefficient"),
        makeParam("CGSO", &MosfetBSIM4Model::CGSO,
                  "gate-source overlap / width"),
        makeParam("CGDO", &MosfetBSIM4Model::CGDO,
                  "gate-drain overlap / width"),
        makeParam("CGBO", &MosfetBSIM4Model::CGBO,
                  "gate-bulk overlap / length"),
        makeParam("TOX", &MosfetBSIM4Model::TOX, "oxide thickness (m)"),
        makeParam("LD", &MosfetBSIM4Model::LD, "lateral diffusion (m)"),
        makeParam("EG", &MosfetBSIM4Model::EG, "bandgap voltage"),
        makeParam("XTI", &MosfetBSIM4Model::XTI,
                  "saturation-current temperature exponent"),
        makeParam("GMIN", &MosfetBSIM4Model::GMIN, "parallel conductance"),
        // BSIM4-specific
        makeParam("U0", &MosfetBSIM4Model::U0, "zero-bias mobility"),
        makeParam("VSAT", &MosfetBSIM4Model::VSAT, "saturation velocity"),
        makeParam("UA", &MosfetBSIM4Model::UA,
                  "1st-order mobility degradation"),
        makeParam("UB", &MosfetBSIM4Model::UB,
                  "2nd-order mobility degradation"),
        makeParam("UC", &MosfetBSIM4Model::UC,
                  "saturation velocity coefficient"),
        makeParam("UD", &MosfetBSIM4Model::UD, "saturation effect coeff."),
        makeParam("A0", &MosfetBSIM4Model::A0, "bulk charge coeff. 1"),
        makeParam("B0", &MosfetBSIM4Model::B0, "bulk charge coeff. 2"),
        makeParam("B1", &MosfetBSIM4Model::B1, "bulk charge coeff. 3"),
        makeParam("KETA", &MosfetBSIM4Model::KETA,
                  "back-bias threshold coeff."),
        makeParam("ALPHA0", &MosfetBSIM4Model::ALPHA0,
                  "substrate current coefficient"),
        makeParam("ALPHA1", &MosfetBSIM4Model::ALPHA1,
                  "substrate current exponential coefficient"),
        makeParam("PRT", &MosfetBSIM4Model::PRT,
                  "Rds temperature coefficient"),
        makeParam("UTE", &MosfetBSIM4Model::UTE,
                  "mobility temperature exponent"),
        makeParam("DVT0", &MosfetBSIM4Model::DVT0,
                  "short-channel threshold coefficient 1"),
        makeParam("DVT1", &MosfetBSIM4Model::DVT1,
                  "short-channel threshold coefficient 2"),
        makeParam("DVT2", &MosfetBSIM4Model::DVT2,
                  "short-channel threshold coefficient 3"),
        makeParam("DVT2W", &MosfetBSIM4Model::DVT2W,
                  "short-channel threshold coeff. width"),
        makeParam("DVTP0", &MosfetBSIM4Model::DVTP0,
                  "pocket implant threshold coeff."),
        makeParam("W0", &MosfetBSIM4Model::W0, "narrow width offset"),
        makeParam("NFS", &MosfetBSIM4Model::NFS,
                  "subthreshold slope coefficient (optional)"),
        makeCelsiusParam("TNOM", &MosfetBSIM4Model::TNOM,
                         "reference temperature"),
        makeCelsiusParam("TEMP", &MosfetBSIM4Model::TEMP, "card temperature"),
    });
}

void MosfetBSIM4Model::applyGuards() {
    EG   = std::max(EG, 0.1);
    PHI  = std::max(PHI, 0.05);
    PB   = std::max(PB, 0.05);
    TOX  = std::max(TOX, 1e-10);
    U0   = std::max(U0, 1e-6);
    VSAT = std::max(VSAT, 1e3);
    UA   = std::max(UA, 0.0);
    UB   = std::max(UB, 0.0);
    KETA = std::max(KETA, 0.0);
    GMIN = std::max(GMIN, 1e-15);
}

bool MosfetBSIM4Model::processParams(std::string& message) {
    if (LEVEL != 2) {
        message = "Only LEVEL=2 (BSIM4) is supported on this card";
        return false;
    }
    if (VTH0 == 0.0) { message = "VTH0 must be nonzero"; return false; }
    if (KP <= 0.0) { message = "KP must be positive"; return false; }
    if (TOX <= 0.0) { message = "TOX must be positive"; return false; }
    if (PHI <= 0.0) { message = "PHI must be positive"; return false; }
    if (TEMP <= 0.0 || TNOM <= 0.0) {
        message = "TEMP/TNOM must be positive (Kelvin)";
        return false;
    }
    applyGuards();
    updateTemperature(TEMP);
    return true;
}

void MosfetBSIM4Model::updateTemperature(double tempKelvin) {
    thermal_ = scaleTo(tempKelvin);
}

MosfetBSIM4Thermal MosfetBSIM4Model::computeThermal(double tempKelvin) const {
    return scaleTo(tempKelvin);
}

MosfetBSIM4Thermal MosfetBSIM4Model::scaleTo(double Temp) const {
    MosfetBSIM4Thermal t;

    t.vt  = Const::KoverQ * Temp;
    double TRatio = Temp / TNOM;

    // Bandgap voltage (SPICE3 bandgap shift).
    double egfet = EG - (7.02e-4 * Temp * Temp) / (Temp + 1108.0);
    t.eg = egfet;

    // Intrinsic carrier concentration (ni).
    t.ni = 1.45e10 * (Temp / 300.15) * std::sqrt(Temp / 300.15)
         * std::exp(EG / (2.0 * Const::KoverQ) *
                    (1.0 / 300.15 - 1.0 / Temp));

    // Zero-bias mobility (temperature dependence via UTE).
    t.u0 = U0 * std::pow(TRatio, -UTE);

    // Saturation velocity (temperature dependence).
    t.vsat = VSAT * (1.0 - 2.0e-3 * (Temp - 300.15));

    // Body-effect coefficient scaled by TOX.
    t.k1ox = GAMMA * TOX / (Const::KoverQ * Temp);

    // Temperature coefficients for threshold (kept for compatibility).
    t.kt1 = 0.0;
    t.kt1l = 0.0;
    t.kt2 = 0.0;
    t.theta0vb0 = 0.0;

    // Mobility degradation coefficients.
    t.ua = UA;
    t.ub = UB;
    t.uc = UC;
    t.ud = UD;

    // Short-channel threshold coefficients.
    t.dvt0 = DVT0;
    t.dvt1 = DVT1;
    t.dvt2 = DVT2;
    t.dvt2w = DVT2W;

    t.gmin = GMIN;
    return t;
}

RCP<DeviceModel> MosfetBSIM4Model::fromCard(const std::string& name,
                                           const ParamMap& params,
                                           bool nmos) {
    ModelCard card;
    card.name = name;
    card.type = nmos ? "NMOS" : "PMOS";
    card.params = params;

    RCP<MosfetBSIM4Model> m = make_rcp<MosfetBSIM4Model>(name, nmos);
    m->LEVEL = static_cast<int>(card.get("LEVEL", 2.0));
    m->VTH0 = card.get("VTH0", m->VTH0);
    m->KP = card.get("KP", m->KP);
    m->LAMBDA = card.get("LAMBDA", m->LAMBDA);
    m->PHI = card.get("PHI", m->PHI);
    m->GAMMA = card.get("GAMMA", m->GAMMA);
    m->CBD = card.get("CBD", m->CBD);
    m->CBS = card.get("CBS", m->CBS);
    m->IS = card.get("IS", m->IS);
    m->PB = card.get("PB", m->PB);
    m->MJ = card.get("MJ", m->MJ);
    m->FC = card.get("FC", m->FC);
    m->CGSO = card.get("CGSO", m->CGSO);
    m->CGDO = card.get("CGDO", m->CGDO);
    m->CGBO = card.get("CGBO", m->CGBO);
    m->TOX = card.get("TOX", m->TOX);
    m->LD = card.get("LD", m->LD);
    m->EG = card.get("EG", m->EG);
    m->XTI = card.get("XTI", m->XTI);
    m->TNOM = card.get("TNOM", m->TNOM);
    m->TEMP = card.get("TEMP", m->TEMP);
    m->GMIN = card.get("GMIN", m->GMIN);
    // BSIM4-specific
    m->U0 = card.get("U0", m->U0);
    m->VSAT = card.get("VSAT", m->VSAT);
    m->UA = card.get("UA", m->UA);
    m->UB = card.get("UB", m->UB);
    m->UC = card.get("UC", m->UC);
    m->UD = card.get("UD", m->UD);
    m->A0 = card.get("A0", m->A0);
    m->B0 = card.get("B0", m->B0);
    m->B1 = card.get("B1", m->B1);
    m->KETA = card.get("KETA", m->KETA);
    m->ALPHA0 = card.get("ALPHA0", m->ALPHA0);
    m->ALPHA1 = card.get("ALPHA1", m->ALPHA1);
    m->PRT = card.get("PRT", m->PRT);
    m->UTE = card.get("UTE", m->UTE);
    m->DVT0 = card.get("DVT0", m->DVT0);
    m->DVT1 = card.get("DVT1", m->DVT1);
    m->DVT2 = card.get("DVT2", m->DVT2);
    m->DVT2W = card.get("DVT2W", m->DVT2W);
    m->DVTP0 = card.get("DVTP0", m->DVTP0);
    m->W0 = card.get("W0", m->W0);
    m->NFS = card.get("NFS", m->NFS);
    if (card.given("TNOM")) m->TNOM = card.get("TNOM", 27.0) + Const::CtoK;
    if (card.given("TEMP")) m->TEMP = card.get("TEMP", 27.0) + Const::CtoK;
    return m;
}

RCP<const MosfetBSIM4Model> MosfetBSIM4Model::defaultNmos() {
    static RCP<MosfetBSIM4Model> instance;
    if (!instance) {
        instance = make_rcp<MosfetBSIM4Model>("$DEFAULT_NMOS_BSIM4", true);
        std::string message;
        if (!instance->processParams(message)) {
            throw std::logic_error("default NMOS BSIM4 model invalid: " +
                                   message);
        }
    }
    return instance;
}

RCP<const MosfetBSIM4Model> MosfetBSIM4Model::defaultPmos() {
    static RCP<MosfetBSIM4Model> instance;
    if (!instance) {
        instance = make_rcp<MosfetBSIM4Model>("$DEFAULT_PMOS_BSIM4", false);
        std::string message;
        if (!instance->processParams(message)) {
            throw std::logic_error("default PMOS BSIM4 model invalid: " +
                                   message);
        }
    }
    return instance;
}

} // namespace e_XSpice

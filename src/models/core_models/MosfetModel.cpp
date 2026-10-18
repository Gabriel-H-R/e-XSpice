/**
 * @author Gabriel Henrique Silva
 * @file MosfetModel.cpp
 * @brief Shared MOSFET model card — level 1 (Shichman-Hodges), NMOS + PMOS.
 */

#include "MosfetModel.h"
#include "../core/DeviceConstants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace e_XSpice {

// ε_SiO2 = 3.9·ε_0
static constexpr double EpsOx = 3.9 * 8.854214871e-12;

MosfetModel::MosfetModel(std::string name, bool nmos)
    : DeviceModel(std::move(name)), nmos_(nmos) {
    // Card parameter table (canonical first, aliases after — canonical
    // wins when a card spells both). Adding a parameter = one line here.
    declareParams({
        makeIntParam("LEVEL", &MosfetModel::LEVEL, "1=Shichman-Hodges"),
        makeParam("VTO", &MosfetModel::VTO, "zero-bias threshold (SPICE sign)"),
        makeAlias("VT0", &MosfetModel::VTO),
        makeParam("KP", &MosfetModel::KP, "transconductance parameter"),
        makeAlias("K", &MosfetModel::KP),
        makeParam("LAMBDA", &MosfetModel::LAMBDA, "channel-length modulation"),
        makeParam("PHI", &MosfetModel::PHI, "surface potential"),
        makeAlias("PHIB", &MosfetModel::PHI),
        makeParam("GAMMA", &MosfetModel::GAMMA, "bulk threshold parameter"),
        makeAlias("K1", &MosfetModel::GAMMA),
        makeParam("CBD", &MosfetModel::CBD, "bulk-drain capacitance"),
        makeAlias("CBDB", &MosfetModel::CBD),
        makeParam("CBS", &MosfetModel::CBS, "bulk-source capacitance"),
        makeAlias("CBSB", &MosfetModel::CBS),
        makeParam("IS", &MosfetModel::IS, "bulk saturation current"),
        makeAlias("JS", &MosfetModel::IS),
        makeParam("PB", &MosfetModel::PB, "bulk-junction potential"),
        makeParam("MJ", &MosfetModel::MJ, "bulk grading coefficient"),
        makeParam("FC", &MosfetModel::FC, "bulk forward-bias coefficient"),
        makeParam("CGSO", &MosfetModel::CGSO, "gate-source overlap / width"),
        makeParam("CGDO", &MosfetModel::CGDO, "gate-drain overlap / width"),
        makeParam("CGBO", &MosfetModel::CGBO, "gate-bulk overlap / length"),
        makeParam("TOX", &MosfetModel::TOX, "oxide thickness"),
        makeAlias("TOXE", &MosfetModel::TOX),
        makeParam("LD", &MosfetModel::LD, "lateral diffusion"),
        makeAlias("LDIFF", &MosfetModel::LD),
        makeParam("EG", &MosfetModel::EG, "bandgap voltage"),
        makeParam("XTI", &MosfetModel::XTI, "saturation-current temp exponent"),
        makeParam("GMIN", &MosfetModel::GMIN, "parallel conductance"),
        makeCelsiusParam("TNOM", &MosfetModel::TNOM, "reference temperature"),
        makeCelsiusParam("TEMP", &MosfetModel::TEMP, "card temperature"),
    });
}

double MosfetModel::coxPrime() const {
    return EpsOx / TOX;
}

void MosfetModel::applyGuards() {
    MJ = std::min(MJ, 0.9);
    FC = std::min(FC, 0.95);
    EG = std::max(EG, 0.1);
    PHI = std::max(PHI, 0.05);
    PB = std::max(PB, 0.05);
}

bool MosfetModel::processParams(std::string& message) {
    if (LEVEL != 1) {
        message = "Only LEVEL=1 (Shichman-Hodges) is supported";
        return false;
    }
    if (KP <= 0.0)  { message = "KP must be positive";  return false; }
    if (TOX <= 0.0) { message = "TOX must be positive"; return false; }
    if (LAMBDA < 0.0) { message = "LAMBDA must be non-negative"; return false; }
    if (GAMMA  < 0.0) { message = "GAMMA must be non-negative";  return false; }
    if (CBD < 0.0) { message = "CBD must be non-negative"; return false; }
    if (CBS < 0.0) { message = "CBS must be non-negative"; return false; }
    if (IS  < 0.0) { message = "IS must be non-negative"; return false; }
    if (CGSO < 0.0 || CGDO < 0.0 || CGBO < 0.0) {
        message = "overlap capacitances must be non-negative";
        return false;
    }
    if (LD < 0.0) { message = "LD must be non-negative"; return false; }
    if (TEMP <= 0.0 || TNOM <= 0.0) {
        message = "TEMP/TNOM must be positive (Kelvin)";
        return false;
    }

    applyGuards();
    updateTemperature(TEMP);
    return true;
}

void MosfetModel::updateTemperature(double tempKelvin) {
    thermal_ = scaleTo(tempKelvin);
}

MosfetThermal MosfetModel::computeThermal(double tempKelvin) const {
    return scaleTo(tempKelvin);
}

MosfetThermal MosfetModel::scaleTo(double Temp) const {
    MosfetThermal t;

    t.vt = Const::KoverQ * Temp;
    t.isat = (IS > 0.0) ? scaleSatCurrent(IS, 1.0, EG, XTI, Temp, TNOM) : 0.0;

    // Surface potential and bulk-junction potential follow the bandgap
    // shift (same construction as the diode VJ(T); PHI/PB play the role
    // of VJ with the card MJ).
    t.phiT = scaleJunctionCap(PHI, 0.0, 0.5, 0.0, Temp, TNOM).vj;
    t.pbT = scaleJunctionCap(PB, MJ, FC, 1.0, Temp, TNOM).vj;
    t.bulkCap = scaleJunctionCap(PB, MJ, FC, 1.0, Temp, TNOM);
    t.coxPrime = EpsOx / TOX;
    t.gmin = GMIN;

    // PNJLIM critical voltage for the bulk diodes (SPICE3f5 value, N = 1).
    const double isat = (t.isat > 0.0) ? t.isat : IS;
    t.vcrit = t.vt * std::log(t.vt / (Const::Root2 * std::max(isat, 1e-99)));

    return t;
}

RCP<const MosfetModel> MosfetModel::defaultNmos() {
    static RCP<MosfetModel> instance;
    if (!instance) {
        instance = make_rcp<MosfetModel>("$DEFAULT_NMOS", true);
        std::string message;
        if (!instance->processParams(message)) {
            throw std::logic_error("default NMOS model invalid: " + message);
        }
    }
    return instance;
}

RCP<const MosfetModel> MosfetModel::defaultPmos() {
    static RCP<MosfetModel> instance;
    if (!instance) {
        instance = make_rcp<MosfetModel>("$DEFAULT_PMOS", false);
        std::string message;
        if (!instance->processParams(message)) {
            throw std::logic_error("default PMOS model invalid: " + message);
        }
    }
    return instance;
}

} // namespace e_XSpice

/**
 * @author Gabriel Henrique Silva
 * @file DiodeModel.cpp
 * @brief Shared PN-junction model card (Xyce level-1 physics).
 */

#include "DiodeModel.h"
#include "../core/DeviceConstants.h"
#include "Junction.h"

#include <cmath>
#include <stdexcept>

namespace e_XSpice {

using std::exp;
using std::log;
using std::min;
using std::max;

DiodeModel::DiodeModel(std::string name) : DeviceModel(std::move(name)) {
    // Card parameter table (canonical first, aliases after — canonical
    // wins when a card spells both). Adding a parameter = one line here.
    declareParams({
        makeParam("IS", &DiodeModel::IS, "saturation current"),
        makeAlias("JS", &DiodeModel::IS),
        makeParam("N", &DiodeModel::N, "emission coefficient"),
        makeParam("RS", &DiodeModel::RS, "ohmic series resistance"),
        makeParam("CJO", &DiodeModel::CJO, "zero-bias capacitance"),
        makeAlias("CJ", &DiodeModel::CJO),
        makeAlias("CJ0", &DiodeModel::CJO),
        makeParam("VJ", &DiodeModel::VJ, "junction potential"),
        makeAlias("PB", &DiodeModel::VJ),
        makeParam("M", &DiodeModel::M, "grading coefficient"),
        makeParam("FC", &DiodeModel::FC, "forward-bias coefficient"),
        makeParam("BV", &DiodeModel::BV, "reverse breakdown voltage"),
        makeAlias("VB", &DiodeModel::BV),
        makeParam("IBV", &DiodeModel::IBV, "breakdown knee current"),
        makeParam("NBV", &DiodeModel::NBV, "breakdown emission coefficient"),
        makeParam("TT", &DiodeModel::TT, "transit time"),
        makeParam("EG", &DiodeModel::EG, "bandgap voltage"),
        makeParam("XTI", &DiodeModel::XTI, "saturation-current temp exponent"),
        makeParam("IKF", &DiodeModel::IKF, "high-injection knee"),
        makeParam("ISR", &DiodeModel::ISR, "recombination saturation current"),
        makeParam("NR", &DiodeModel::NR, "recombination emission coefficient"),
        makeParam("GMIN", &DiodeModel::GMIN, "parallel conductance"),
        makeCelsiusParam("TNOM", &DiodeModel::TNOM, "reference temperature"),
        makeCelsiusParam("TEMP", &DiodeModel::TEMP, "card temperature"),
    });
}

void DiodeModel::applyGuards() {
    // History-dependent parameter guards (Xyce Model::processParams).
    M  = min(M, 0.9);
    EG = max(EG, 0.1);
    FC = min(FC, 0.95);
}

bool DiodeModel::processParams(std::string& message) {
    if (IS <= 0.0)  { message = "IS must be positive"; return false; }
    if (N  <= 0.0)  { message = "N must be positive"; return false; }
    if (CJO < 0.0)  { message = "CJO must be non-negative"; return false; }
    if (VJ <= 0.0)  { message = "VJ must be positive"; return false; }
    if (TEMP <= 0.0 || TNOM <= 0.0) {
        message = "TEMP/TNOM must be positive (Kelvin)";
        return false;
    }

    applyGuards();
    updateTemperature(TEMP);
    return true;
}

void DiodeModel::updateTemperature(double tempKelvin) {
    thermal_ = scaleTo(tempKelvin);
}

DiodeThermal DiodeModel::computeThermal(double tempKelvin) const {
    return scaleTo(tempKelvin);
}

DiodeThermal DiodeModel::scaleTo(double Temp) const {
    DiodeThermal t;

    const double vt = Const::KoverQ * Temp;
    t.vt     = vt;
    t.vte    = N  * vt;
    t.vteR   = NR * vt;
    t.vteBRK = (NBV > 0.0 ? NBV : 1.0) * vt;

    // Junction potential/capacitance at temperature + continuity constants
    // (shared with every other junction device, see Junction.h).
    const JunctionCap cap = scaleJunctionCap(VJ, M, FC, CJO, Temp, TNOM);
    t.tJctPot = cap.vj;
    t.tJctCap = cap.cjo;
    t.tDepCap = cap.tdep;
    t.F2      = cap.f2;
    t.F3      = cap.f3;
    t.tF1     = cap.f1;

    // saturation currents at temperature (SPICE3 formula)
    t.isat  = scaleSatCurrent(IS,  N,  EG, XTI, Temp, TNOM);
    t.isatR = scaleSatCurrent(ISR, NR, EG, XTI, Temp, TNOM);

    t.gspr = (RS > 0.0) ? 1.0 / RS : 0.0;
    t.ikf  = IKF;

    // PNJLIM critical voltage: vcrit = vte·ln(vte/(√2·IS))
    t.tVcrit = t.vte * std::log(t.vte / (Const::Root2 * t.isat));

    // breakdown fit (only when BV was given)
    t.tBrkdwnV = 0.0;
    if (BVGiven) {
        t.tBrkdwnV = fitBreakdown(Temp, vt, t.isat);
    }
    return t;
}

double DiodeModel::fitBreakdown(double Temp, double vt, double isat) const {
    // Fixed-point fit of tBrkdwnV so that I_D = −IBV at Vd = −BV (Xyce).
    const double tempBV = BV * (1.0 + (Temp - TNOM) * (0.0 + 0.0 * (Temp - TNOM)));
    const double reltol = 1.0e-3;

    const double cbv = IBV;
    double xbv = 0.0;

    if (cbv < isat * tempBV / vt) {
        xbv = tempBV;
    } else {
        const double tol = reltol * cbv;
        xbv = tempBV - NBV * vt * log(1.0 + cbv / isat);
        for (int i = 0; i < 25; ++i) {
            xbv = tempBV - NBV * vt * log(cbv / isat + 1.0 - xbv / vt);
            const double xcbv =
                isat * (exp((tempBV - xbv) / (NBV * vt)) - 1.0 + xbv / vt);
            if (std::abs(xcbv - cbv) <= tol) break;
        }
    }
    return xbv;
}

RCP<const DiodeModel> DiodeModel::defaultModel() {
    static RCP<DiodeModel> instance;
    if (!instance) {
        instance = make_rcp<DiodeModel>("$DEFAULT_D");
        std::string message;
        if (!instance->processParams(message)) {
            throw std::logic_error("default diode model invalid: " + message);
        }
    }
    return instance;
}

} // namespace e_XSpice

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

RCP<DeviceModel> DiodeModel::fromCard(const std::string& name,
                                     const ParamMap& params) {
    ModelCard card;
    card.name = name;
    card.type = "D";
    card.params = params;

    RCP<DiodeModel> m = make_rcp<DiodeModel>(name);
    m->IS  = card.get("IS", m->IS);
    m->IS  = card.get("JS", m->IS);
    m->N   = card.get("N", m->N);
    m->RS  = card.get("RS", m->RS);
    m->CJO = card.get("CJO", m->CJO);
    m->CJO = card.get("CJ", m->CJO);
    m->CJO = card.get("CJ0", m->CJO);
    m->VJ  = card.get("VJ", m->VJ);
    m->VJ  = card.get("PB", m->VJ);
    m->M   = card.get("M", m->M);
    m->FC  = card.get("FC", m->FC);
    if (card.given("BV") || card.given("VB")) {
        m->BV = card.get("BV", card.get("VB", m->BV));
        m->BVGiven = true;
    }
    m->IBV = card.get("IBV", m->IBV);
    m->NBV = card.get("NBV", m->NBV);
    m->TT  = card.get("TT", m->TT);
    m->EG  = card.get("EG", m->EG);
    m->XTI = card.get("XTI", m->XTI);
    m->IKF = card.get("IKF", m->IKF);
    m->ISR = card.get("ISR", m->ISR);
    m->NR  = card.get("NR", m->NR);
    m->GMIN = card.get("GMIN", m->GMIN);
    // SPICE writes TNOM/TEMP in Celsius; the card stores Kelvin.
    if (card.given("TNOM")) m->TNOM = card.get("TNOM", 27.0) + Const::CtoK;
    if (card.given("TEMP")) m->TEMP = card.get("TEMP", 27.0) + Const::CtoK;
    return m;
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

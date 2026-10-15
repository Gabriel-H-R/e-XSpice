/**
 * @author Gabriel Henrique Silva
 * @file PassiveModels.h
 * @brief Model cards for passive devices (R, C, L).
 *
 * Even linear passives benefit from the model/instance split: the card
 * carries manufacturing/temperature behavior shared by every instance
 * (multiplier, TC1/TC2, TNOM), while the instance line keeps its own
 * nominal value. A null model handle means "ideal, no temperature
 * behavior" — every existing netlist keeps working unchanged.
 *
 * Effective value with a model (Xyce convention):
 *   R(T) = R_inst · Rmult · (1 + TC1·ΔT + TC2·ΔT²),  ΔT = T − TNOM
 * (same form for C and L). T is in Kelvin; TNOM defaults to 300.15 K.
 */

#ifndef PASSIVE_MODELS_H
#define PASSIVE_MODELS_H

#include "../core/DeviceModel.h"
#include "../core/ModelCard.h"
#include <string>

namespace e_XSpice {

/**
 * @brief Shared resistor model card.
 */
class ResistorModel : public DeviceModel {
public:
    explicit ResistorModel(std::string name) : DeviceModel(std::move(name)) {}

    std::string typeName() const override { return "R"; }

    double Rmult = 1.0;    // resistance multiplier (model card R=)
    double TC1   = 0.0;    // linear temperature coefficient (1/K)
    double TC2   = 0.0;    // quadratic temperature coefficient (1/K²)
    double TNOM  = 300.15; // reference temperature (K)

    bool processParams(std::string& message) override {
        if (Rmult <= 0.0) { message = "R multiplier must be positive"; return false; }
        return true;
    }

    /**
     * @brief Temperature/load factor applied to the instance value.
     */
    double factor(double tempKelvin) const {
        const double dT = tempKelvin - TNOM;
        return Rmult * (1.0 + TC1 * dT + TC2 * dT * dT);
    }

    static RCP<DeviceModel> fromCard(const std::string& name,
                                     const ParamMap& params) {
        ModelCard card;
        card.name = name;
        card.type = "R";
        card.params = params;
        RCP<ResistorModel> m = make_rcp<ResistorModel>(name);
        m->Rmult = card.get("R", m->Rmult);
        m->TC1   = card.get("TC1", m->TC1);
        m->TC2   = card.get("TC2", m->TC2);
        if (card.given("TNOM")) m->TNOM = card.get("TNOM", 27.0) + 273.15;
        return m;
    }
};

/**
 * @brief Shared capacitor model card (wiring follows the resistor path).
 */
class CapacitorModel : public DeviceModel {
public:
    explicit CapacitorModel(std::string name) : DeviceModel(std::move(name)) {}

    std::string typeName() const override { return "C"; }

    double Cmult = 1.0;
    double TC1   = 0.0;
    double TC2   = 0.0;
    double TNOM  = 300.15;

    bool processParams(std::string& message) override {
        if (Cmult <= 0.0) { message = "C multiplier must be positive"; return false; }
        return true;
    }

    double factor(double tempKelvin) const {
        const double dT = tempKelvin - TNOM;
        return Cmult * (1.0 + TC1 * dT + TC2 * dT * dT);
    }

    static RCP<DeviceModel> fromCard(const std::string& name,
                                     const ParamMap& params) {
        ModelCard card;
        card.name = name;
        card.type = "C";
        card.params = params;
        RCP<CapacitorModel> m = make_rcp<CapacitorModel>(name);
        m->Cmult = card.get("C", m->Cmult);
        m->TC1   = card.get("TC1", m->TC1);
        m->TC2   = card.get("TC2", m->TC2);
        if (card.given("TNOM")) m->TNOM = card.get("TNOM", 27.0) + 273.15;
        return m;
    }
};

/**
 * @brief Shared inductor model card (wiring follows the resistor path).
 */
class InductorModel : public DeviceModel {
public:
    explicit InductorModel(std::string name) : DeviceModel(std::move(name)) {}

    std::string typeName() const override { return "L"; }

    double Lmult = 1.0;
    double TC1   = 0.0;
    double TC2   = 0.0;
    double TNOM  = 300.15;

    bool processParams(std::string& message) override {
        if (Lmult <= 0.0) { message = "L multiplier must be positive"; return false; }
        return true;
    }

    double factor(double tempKelvin) const {
        const double dT = tempKelvin - TNOM;
        return Lmult * (1.0 + TC1 * dT + TC2 * dT * dT);
    }

    static RCP<DeviceModel> fromCard(const std::string& name,
                                     const ParamMap& params) {
        ModelCard card;
        card.name = name;
        card.type = "L";
        card.params = params;
        RCP<InductorModel> m = make_rcp<InductorModel>(name);
        m->Lmult = card.get("L", m->Lmult);
        m->TC1   = card.get("TC1", m->TC1);
        m->TC2   = card.get("TC2", m->TC2);
        if (card.given("TNOM")) m->TNOM = card.get("TNOM", 27.0) + 273.15;
        return m;
    }
};

} // namespace e_XSpice

#endif // PASSIVE_MODELS_H

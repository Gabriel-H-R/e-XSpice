/**
 * @author Gabriel Henrique Silva
 * @file PassiveModels.cpp
 * @brief Model cards for passive devices (R, C, L) — constructors with the
 *        declarative parameter tables (see ModelParams.h).
 */

#include "PassiveModels.h"

namespace e_XSpice {

ResistorModel::ResistorModel(std::string name) : DeviceModel(std::move(name)) {
    declareParams({
        makeParam("R", &ResistorModel::Rmult, "resistance multiplier"),
        makeParam("TC1", &ResistorModel::TC1, "linear temp coefficient"),
        makeParam("TC2", &ResistorModel::TC2, "quadratic temp coefficient"),
        makeCelsiusParam("TNOM", &ResistorModel::TNOM, "reference temperature"),
    });
}

CapacitorModel::CapacitorModel(std::string name) : DeviceModel(std::move(name)) {
    declareParams({
        makeParam("C", &CapacitorModel::Cmult, "capacitance multiplier"),
        makeParam("TC1", &CapacitorModel::TC1, "linear temp coefficient"),
        makeParam("TC2", &CapacitorModel::TC2, "quadratic temp coefficient"),
        makeCelsiusParam("TNOM", &CapacitorModel::TNOM, "reference temperature"),
    });
}

InductorModel::InductorModel(std::string name) : DeviceModel(std::move(name)) {
    declareParams({
        makeParam("L", &InductorModel::Lmult, "inductance multiplier"),
        makeParam("TC1", &InductorModel::TC1, "linear temp coefficient"),
        makeParam("TC2", &InductorModel::TC2, "quadratic temp coefficient"),
        makeCelsiusParam("TNOM", &InductorModel::TNOM, "reference temperature"),
    });
}

} // namespace e_XSpice

/**
 * @author Gabriel Henrique Silva
 * @file ModelRegistry.cpp
 * @brief Central registry of model cards + typed model factories.
 */

#include "ModelRegistry.h"
#include "../core_models/DiodeModel.h"
#include "../core_models/PassiveModels.h"
#include "../core_models/BjtModel.h"
#include "../core_models/MosfetModel.h"
#include "../core_models/MosfetBSIM4Model.h"

#include <cctype>
#include <stdexcept>

namespace e_XSpice {

ModelRegistry::ModelRegistry() {
    // Built-in model factories (Xyce-style startup registration).
    // A new model = one subclass (with its parameter table) + one line
    // here. Factories only construct defaults; the registry fills card
    // values generically, so factories never parse and never change when
    // parameters are added.
    registerModelType("D", [](const std::string& n) {
        return make_rcp<DiodeModel>(n);
    });
    registerModelType("R", [](const std::string& n) {
        return make_rcp<ResistorModel>(n);
    });
    registerModelType("C", [](const std::string& n) {
        return make_rcp<CapacitorModel>(n);
    });
    registerModelType("L", [](const std::string& n) {
        return make_rcp<InductorModel>(n);
    });
    registerModelType("NPN", [](const std::string& n) {
        return make_rcp<BjtModel>(n, false);
    });
    registerModelType("PNP", [](const std::string& n) {
        return make_rcp<BjtModel>(n, true);
    });
    registerModelType("NMOS", [](const std::string& n) {
        return make_rcp<MosfetModel>(n, true);
    });
    registerModelType("PMOS", [](const std::string& n) {
        return make_rcp<MosfetModel>(n, false);
    });
    registerModelType("NMOS2", [](const std::string& n) {
        return make_rcp<MosfetBSIM4Model>(n, true);
    });
    registerModelType("PMOS2", [](const std::string& n) {
        return make_rcp<MosfetBSIM4Model>(n, false);
    });
    registerModelType("NPN2", [](const std::string& n) {
        return make_rcp<BjtModel>(n, false);
    });
    registerModelType("PNP2", [](const std::string& n) {
        return make_rcp<BjtModel>(n, true);
    });
}

std::string ModelRegistry::upper(const std::string& s) {
    std::string out = s;
    for (auto& c : out) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return out;
}

void ModelRegistry::registerModelType(const std::string& type,
                                      ModelFactory factory) {
    factories_[upper(type)] = std::move(factory);
}

bool ModelRegistry::hasModelType(const std::string& type) const {
    return factories_.find(upper(type)) != factories_.end();
}

std::vector<std::string> ModelRegistry::registeredTypes() const {
    std::vector<std::string> types;
    types.reserve(factories_.size());
    for (const auto& kv : factories_) types.push_back(kv.first);
    return types;   // std::map order = sorted
}

static std::string joinNames(const std::vector<std::string>& names) {
    std::string out;
    for (const auto& n : names) {
        if (!out.empty()) out += " ";
        out += n;
    }
    return out;
}

RCP<DeviceModel> ModelRegistry::createModel(const ModelCard& card) {
    const std::string type = upper(card.type);
    const std::string name = upper(card.name);

    if (name.empty()) {
        throw std::invalid_argument(".MODEL requires a name");
    }

    const auto fit = factories_.find(type);
    if (fit == factories_.end()) {
        throw std::invalid_argument("Unknown .MODEL type '" + card.type +
                                    "' (registered types: " +
                                    joinNames(registeredTypes()) + ")");
    }
    if (models_.find(name) != models_.end()) {
        throw std::invalid_argument("Duplicate .MODEL '" + name + "'");
    }

    RCP<DeviceModel> model = fit->second(name);

    // Generic parameter fill from the declared table (aliases, °C keys,
    // canonical-wins precedence). Unknown keys fail with the valid list —
    // a typo'd parameter must never silently become a default.
    const std::vector<std::string> unknown =
        model->applyCardParams(card.params);
    if (!unknown.empty()) {
        throw std::invalid_argument(".MODEL '" + name + "' (" + type +
                                    "): unknown parameter(s) '" +
                                    joinNames(unknown) +
                                    "' (valid parameters: " +
                                    model->validParamList() + ")");
    }

    std::string message;
    if (!model->processParams(message)) {
        throw std::invalid_argument(".MODEL '" + name + "': " + message);
    }

    models_[name] = model;   // shared ownership (one atomic inc, setup)
    return model;
}

void ModelRegistry::addModel(RCP<DeviceModel> model) {
    if (!model) {
        throw std::invalid_argument("Cannot register a null model");
    }
    const std::string name = upper(model->name());
    if (name.empty()) {
        throw std::invalid_argument("Model requires a name");
    }
    if (models_.find(name) != models_.end()) {
        throw std::invalid_argument("Duplicate model '" + name + "'");
    }
    models_[name] = std::move(model);
}

RCP<DeviceModel> ModelRegistry::find(const std::string& name) const {
    const auto it = models_.find(upper(name));
    if (it == models_.end()) return RCP<DeviceModel>();
    return it->second;   // copy handle (one atomic inc, setup)
}

std::vector<std::string> ModelRegistry::modelNames() const {
    std::vector<std::string> names;
    names.reserve(models_.size());
    for (const auto& kv : models_) names.push_back(kv.first);
    return names;
}

} // namespace e_XSpice

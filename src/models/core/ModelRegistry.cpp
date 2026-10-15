/**
 * @author Gabriel Henrique Silva
 * @file ModelRegistry.cpp
 * @brief Central registry of model cards + typed model factories.
 */

#include "ModelRegistry.h"
#include "../core_models/DiodeModel.h"
#include "../core_models/PassiveModels.h"
#include "../core_models/BjtModel.h"

#include <cctype>
#include <stdexcept>

namespace e_XSpice {

ModelRegistry::ModelRegistry() {
    // Built-in model factories (Xyce-style startup registration).
    // Custom physics (MOSFET, JFET, ...) registers here at runtime.
    registerModelType("D", &DiodeModel::fromCard);
    registerModelType("R", &ResistorModel::fromCard);
    registerModelType("C", &CapacitorModel::fromCard);
    registerModelType("L", &InductorModel::fromCard);
    registerModelType("NPN", [](const std::string& n, const ParamMap& p) {
        return BjtModel::fromCard(n, p, false);
    });
    registerModelType("PNP", [](const std::string& n, const ParamMap& p) {
        return BjtModel::fromCard(n, p, true);
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

RCP<DeviceModel> ModelRegistry::createModel(const ModelCard& card) {
    const std::string type = upper(card.type);
    const std::string name = upper(card.name);

    if (name.empty()) {
        throw std::invalid_argument(".MODEL requires a name");
    }

    const auto fit = factories_.find(type);
    if (fit == factories_.end()) {
        throw std::invalid_argument("Unknown .MODEL type '" + card.type +
                                    "' (no factory registered for '" + type + "')");
    }
    if (models_.find(name) != models_.end()) {
        throw std::invalid_argument("Duplicate .MODEL '" + name + "'");
    }

    RCP<DeviceModel> model = fit->second(name, card.params);

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

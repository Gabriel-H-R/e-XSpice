/**
 * @author Gabriel Henrique Silva
 * @file ModelRegistry.h
 * @brief Central registry of model cards + typed model factories.
 *
 * Two jobs, both O(log n):
 *  1. Own every `.MODEL` card of a netlist (`RCP<DeviceModel>` by UPPER name)
 *     and hand out shared handles to referencing devices.
 *  2. Convert generic ModelCards into typed models through factories
 *     registered per device letter. Unknown letters fail with a clear
 *     message instead of a parser switch-case — this is what makes the
 *     package open for extension (MOSFET/BJT = new factory, zero parser
 *     changes).
 *
 * Lifetime: the registry typically belongs to the NetlistParser. Devices
 * keep their own RCP handle, so discarding the registry after parsing can
 * never dangle a running simulation.
 */

#ifndef MODEL_REGISTRY_H
#define MODEL_REGISTRY_H

#include "DeviceModel.h"
#include "ModelCard.h"
#include <functional>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

namespace e_XSpice {

/**
 * @brief Builds a typed model with default parameters.
 *        The registry fills the card values generically afterwards
 *        (DeviceModel::applyCardParams), so factories never parse.
 *        Throws std::invalid_argument on invalid parameter combinations.
 */
using ModelFactory =
    std::function<RCP<DeviceModel>(const std::string& name)>;

class ModelRegistry {
public:
    /**
     * @brief Constructs the registry with the built-in model factories
     *        (D, R, C, L) already registered — Xyce-style startup
     *        registration, no parser involvement needed.
     */
    ModelRegistry();

    /**
     * @brief Registers a factory for a device letter ("D", "R", ...).
     *        Re-registering replaces the previous factory.
     */
    void registerModelType(const std::string& type, ModelFactory factory);

    /**
     * @brief True when a factory exists for this device letter.
     */
    bool hasModelType(const std::string& type) const;

    /**
     * @brief Registered model types ("D", "NPN", "NMOS", ...), sorted.
     *        Used in unknown-type errors so a typo shows what's legal.
     */
    std::vector<std::string> registeredTypes() const;

    /**
     * @brief Creates a typed model from a card (factory + generic parameter
     *        fill + processParams) and stores it. Throws on unknown type
     *        (lists registered types), duplicate name, unknown parameter
     *        (lists the valid parameters of that type), or invalid
     *        parameter combinations.
     * @return shared handle to the stored model
     */
    RCP<DeviceModel> createModel(const ModelCard& card);

    /**
     * @brief Stores an already-built model (programmatic API).
     */
    void addModel(RCP<DeviceModel> model);

    /**
     * @brief Looks a model up by name (case-insensitive). Null RCP when
     *        absent — the caller reports the "unknown model" error.
     */
    RCP<DeviceModel> find(const std::string& name) const;

    /**
     * @brief Typed lookup with a checked downcast. Null when absent OR when
     *        the card is of a different type ("model X is not a D model").
     */
    template <typename M>
    RCP<const M> findAs(const std::string& name) const {
        static_assert(std::is_base_of<DeviceModel, M>::value,
                      "findAs requires a DeviceModel subclass");
        RCP<DeviceModel> base = find(name);
        if (!base) return RCP<const M>();
        M* typed = dynamic_cast<M*>(base.get());
        if (!typed) return RCP<const M>();
        // Share ownership with the stored handle (one atomic inc, setup).
        return RCP<const M>(static_cast<const M*>(base.get()));
    }

    /**
     * @brief Number of stored models.
     */
    size_t size() const { return models_.size(); }

    /**
     * @brief Names of all stored models (debugging / .PRINT support).
     */
    std::vector<std::string> modelNames() const;

private:
    static std::string upper(const std::string& s);

    std::map<std::string, ModelFactory> factories_;   // "D" -> factory
    std::map<std::string, RCP<DeviceModel>> models_;  // NAME -> model
};

} // namespace e_XSpice

#endif // MODEL_REGISTRY_H

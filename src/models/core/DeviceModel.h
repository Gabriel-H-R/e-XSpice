/**
 * @author Gabriel Henrique Silva
 * @file DeviceModel.h
 * @brief Abstract base class for all device models (model cards).
 *
 * ================================================================
 * MODEL vs INSTANCE (Xyce architecture, modernized)
 * ================================================================
 *
 * A *model* (`.MODEL name TYPE (...)`) holds the parameters SHARED by every
 * device instance that references it: physics constants, temperature
 * behavior, geometry scaling rules. An *instance* (one netlist line) holds
 * only its own terminals plus per-instance overrides (AREA, M, IC, ...).
 *
 * Why this split matters:
 *  1. Memory: one model object serves N instances (no per-instance copies
 *     of dozens of doubles).
 *  2. Performance: temperature-dependent derived quantities are computed
 *     ONCE per model (updateTemperature) instead of once per instance.
 *  3. Extensibility: new physics = new DeviceModel subclass + one factory
 *     registration. Devices stay agnostic: they only read through a
 *     `RCP<const SpecificModel>` handle.
 *
 * Ownership: models live in RCP<T> handles shared between the
 * ModelRegistry and every referencing device. Dropping the registry never
 * invalidates a running simulation — textbook shared ownership, setup-time
 * atomic cost only, zero traffic on solver hot paths (devices read through
 * raw observers / cached plain doubles).
 */

#ifndef DEVICE_MODEL_H
#define DEVICE_MODEL_H

#include "e_XSpice_rcp.h"   // RCP<T>: intrusive refcounted ownership
#include "ModelParams.h"     // declarative parameter tables (addPar-style)
#include "ModelCard.h"       // ParamMap (card key=value pairs)
#include "DeviceConstants.h" // Const::CtoK (°C → K for TNOM/TEMP keys)
#include <string>
#include <vector>

namespace e_XSpice {

/**
 * @brief Base class for every model card (D, R, C, L, Q, M, ...).
 */
class DeviceModel : public RefCounted {
public:
    explicit DeviceModel(std::string name) : name_(std::move(name)) {}
    virtual ~DeviceModel() = default;

    // Non-copyable (shared by handle, never by value).
    DeviceModel(const DeviceModel&) = delete;
    DeviceModel& operator=(const DeviceModel&) = delete;

    /**
     * @brief Model card name as written in `.MODEL` (stored UPPER-CASED).
     */
    const std::string& name() const { return name_; }

    /**
     * @brief SPICE device letter this model applies to ("D", "R", ...).
     */
    virtual std::string typeName() const = 0;

    /**
     * @brief Validates parameters, applies guards and precomputes every
     *        temperature-dependent derived quantity. Called once after the
     *        card is filled (and again if temperature stepping changes T).
     * @return true on success (false + message on invalid combination)
     */
    virtual bool processParams(std::string& message) {
        (void)message;
        return true;
    }

    /**
     * @brief Recomputes temperature-dependent internals for T (Kelvin).
     *        Called by processParams(); re-called only when T changes, so
     *        per-instance evaluation never pays for temperature math.
     */
    virtual void updateTemperature(double tempKelvin) { (void)tempKelvin; }

    /**
     * @brief True when the card carried its own `TEMP=`, which by SPICE
     *        convention wins over the circuit temperature (`.TEMP`). Cards
     *        without it follow the circuit temperature.
     */
    bool cardTempGiven() const { return card_temp_given_; }

    /**
     * @brief Applies the circuit temperature to this card.
     *
     * Called once per `.TEMP`/temperature-sweep value. A no-op when the
     * card pinned its own TEMP (SPICE precedence: instance TEMP > card TEMP
     * > circuit TEMP; instances are handled by the device, which owns its
     * own thermal snapshot and never asks the card).
     */
    virtual void setOperatingTemperature(double tempKelvin) {
        if (card_temp_given_) return;
        updateTemperature(tempKelvin);
    }

    /**
     * @brief Returns true if this model card is a BSIM4 level-2 card.
     *        Base class returns false; MosfetBSIM4Model overrides to true.
     */
    virtual bool isBSIM4() const { return false; }

    // ---------------------------------------------------------------
    // Declarative parameters (see ModelParams.h)
    // ---------------------------------------------------------------

    /**
     * @brief Looks a card key up in the declared table (canonical+aliases).
     */
    bool hasParam(const std::string& key) const {
        for (const auto& d : params_) {
            if (d.name == key) return true;
        }
        return false;
    }

    /**
     * @brief Read access to the declared table (introspection / docs).
     */
    const std::vector<ModelParamDesc>& paramTable() const { return params_; }

    /**
     * @brief Fills every declared parameter present in `params`.
     *        Aliases apply first so the canonical spelling always wins.
     *        Celsius keys get the °C→K offset; absent keys keep defaults.
     * @return card keys with no table entry (unknown/typo), in card order
     */
    std::vector<std::string> applyCardParams(const ParamMap& params) {
        std::vector<std::string> unknown =
            applyParamTable(*this, params_, params);
        for (const auto& kv : params) {
            if (hasParam(kv.first)) {
                // TEMP on a card pins its temperature against `.TEMP`
                // (SPICE precedence). One generic rule in the base class,
                // so every card type inherits it without extra code.
                if (kv.first == "TEMP") card_temp_given_ = true;
                onParamGiven(kv.first);
            }
        }
        return unknown;
    }

    /**
     * @brief Human-readable parameter list (canonical first, then aliases).
     *        Used in unknown-parameter errors so a typo shows what's legal.
     */
    std::string validParamList() const {
        return e_XSpice::validParamList(params_);
    }

protected:
    /**
     * @brief Declares the card parameter table (called once, subclass ctor).
     */
    void declareParams(std::vector<ModelParamDesc> table) {
        params_ = std::move(table);
    }

    /**
     * @brief Set when the card carried TEMP= (see cardTempGiven()).
     */
    bool card_temp_given_ = false;

    /**
     * @brief Fires once per APPLIED card key (canonical or alias spelling).
     *        Backs `given` side effects (e.g. DiodeModel sets BVGiven when
     *        "BV" or "VB" arrives) without exposing the raw card.
     */
    virtual void onParamGiven(const std::string& key) { (void)key; }

    std::string name_;

private:
    std::vector<ModelParamDesc> params_;
};

} // namespace e_XSpice

#endif // DEVICE_MODEL_H

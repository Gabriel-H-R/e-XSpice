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
#include <string>

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

protected:
    std::string name_;
};

} // namespace e_XSpice

#endif // DEVICE_MODEL_H

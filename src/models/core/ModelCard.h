/**
 * @author Gabriel Henrique Silva
 * @file ModelCard.h
 * @brief Generic parsed `.MODEL` card (type-erased parameters).
 *
 * The netlist parser produces ModelCards without knowing any device
 * physics: `.MODEL DMOD D (IS=1e-14 N=1.0)` becomes
 * `{name="DMOD", type="D", params={IS:1e-14, N:1.0}}`. The ModelRegistry
 * then converts the card into a typed DeviceModel through the factory
 * registered for that letter. Adding a new device type (Q, M, ...) never
 * touches the parser — only a factory registration.
 *
 * All keys are stored UPPER-CASED; values are already unit-scaled doubles
 * (the parser applies k/MEG/u/n/... suffixes).
 */

#ifndef MODEL_CARD_H
#define MODEL_CARD_H

#include <map>
#include <string>

namespace e_XSpice {

/**
 * @brief Parameter dictionary of a model card (UPPER-CASE keys).
 */
using ParamMap = std::map<std::string, double>;

/**
 * @brief Type-erased model card straight from the netlist.
 */
struct ModelCard {
    std::string name;    // card name (upper-cased by the parser)
    std::string type;    // device letter (upper-cased: "D", "R", "C", ...)
    ParamMap    params;  // key=value pairs

    /**
     * @brief Reads a key with a default (missing key, not an error).
     */
    double get(const std::string& key, double fallback) const {
        const auto it = params.find(key);
        return (it != params.end()) ? it->second : fallback;
    }

    /**
     * @brief True when the key was explicitly given on the card.
     *        Backs the SPICE `given` semantics (BVGiven, ICGiven, ...).
     */
    bool given(const std::string& key) const {
        return params.find(key) != params.end();
    }
};

} // namespace e_XSpice

#endif // MODEL_CARD_H

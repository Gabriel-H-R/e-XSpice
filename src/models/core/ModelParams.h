/**
 * @author Gabriel Henrique Silva
 * @file ModelParams.h
 * @brief Declarative parameter tables (Xyce `addPar` idea, modernized).
 *
 * Every parameter owner declares its keys ONCE, as data instead of code —
 * model cards in the model ctor, per-instance line options in the instance
 * struct:
 *
 *     declareParams({
 *         makeParam("IS", &DiodeModel::IS, "saturation current"),
 *         makeAlias("JS", &DiodeModel::IS),          // SPICE alias
 *         makeCelsiusParam("TEMP", &DiodeModel::TEMP, "card temperature"),
 *     });
 *
 * One generic fill (applyParamTable) serves BOTH levels — Xyce's
 * ParametricData serves loadModelParameters and loadInstanceParameters
 * alike. Consequences:
 *
 *  1. New models/devices: declare the tables + write processParams().
 *     Factories stay one-line lambdas, the parser never changes, and
 *     unknown keys fail with the VALID list instead of being silently
 *     ignored (catches `BTEA`-vs-`BETA` typos at parse time, SPICE-style).
 *  2. Precedence is deterministic and documented: when both a canonical
 *     key and one of its aliases are present, the CANONICAL value wins.
 *  3. `celsius` keys (TNOM/TEMP) arrive in °C (SPICE convention) and are
 *     stored Kelvin; absent keys keep the (already Kelvin) defaults.
 *  4. `given` semantics survive two ways: the `noted` hook per descriptor
 *     (instance TEMPGiven) and the DeviceModel::onParamGiven virtual
 *     (card BVGiven) for the card level.
 *
 * Cost: one pointer + a few strings per entry, setup-time only. Hot paths
 * never touch tables (devices read plain member doubles).
 */

#ifndef MODEL_PARAMS_H
#define MODEL_PARAMS_H

#include <functional>
#include <string>
#include <type_traits>
#include <vector>

#include "ModelCard.h"   // ParamMap (standalone, no cycle)

namespace e_XSpice {

class DeviceModel;   // full type in DeviceModel.h (setters cast there)

template <class O> struct TypedParamDesc;   // defined below

/**
 * @brief One declarative parameter, bound to an owner object of type O.
 *
 * O is DeviceModel for card parameters (alias ModelParamDesc below) and an
 * instance-params struct (Diode/Bjt/MosfetInstanceParams) for per-device
 * line options.
 */
template <class O>
struct TypedParamDesc {
    TypedParamDesc() = default;
    std::string name;    // UPPER-CASE key ("IS", "AREA", "TEMP", ...)
    std::string help;    // one-line doc (shown in valid-parameter listings)
    bool isAlias = false;
    bool celsius = false;   // value in °C on the card/line, stored Kelvin
    std::function<void(O&, double)> set;
    std::function<double(const O&)> get;
    // Optional side effect on apply (e.g. instance TEMPGiven = true).
    // Card-level models use the DeviceModel::onParamGiven virtual instead.
    std::function<void(O&)> noted = nullptr;

    /**
     * @brief Narrows a subclass descriptor to the card-level one, so model
     *        tables (written with the concrete member-holder type) convert
     *        elementwise into vector<ModelParamDesc>. Enabled only when M
     *        derives from DeviceModel; instance tables need no conversion.
     */
    template <class M,
              typename std::enable_if<std::is_base_of<DeviceModel, M>::value,
                                      int>::type = 0>
    TypedParamDesc(const TypedParamDesc<M>& d)
        : name(d.name)
        , help(d.help)
        , isAlias(d.isAlias)
        , celsius(d.celsius)
        , set([s = d.set](DeviceModel& m, double v) {
              s(static_cast<M&>(m), v);
          })
        , get([g = d.get](const DeviceModel& m) -> double {
              return g(static_cast<const M&>(m));
          }) {
        if (d.noted) {
            auto n = d.noted;
            noted = [n](DeviceModel& m) { n(static_cast<M&>(m)); };
        }
    }
};

/**
 * @brief Card-level descriptor (owner is always DeviceModel here).
 */
using ModelParamDesc = TypedParamDesc<DeviceModel>;

/**
 * @brief Canonical parameter bound to a `double` member of any owner.
 *        Optional 4th argument: a `noted` side-effect hook.
 */
template <typename O>
inline TypedParamDesc<O> makeParam(const char* name, double O::* member,
                                  const char* help = "") {
    TypedParamDesc<O> d;
    d.name = name;
    d.help = help ? help : "";
    d.set = [member](O& o, double v) { (o.*member) = v; };
    d.get = [member](const O& o) -> double { return (o.*member); };
    return d;
}

template <typename O, typename N>
inline TypedParamDesc<O> makeParam(const char* name, double O::* member,
                                  const char* help, N noted) {
    TypedParamDesc<O> d = makeParam(name, member, help);
    d.noted = std::move(noted);
    return d;
}

/**
 * @brief SPICE alias for a `double` member (applied first; the canonical
 *        spelling always wins when both are present).
 */
template <typename O>
inline TypedParamDesc<O> makeAlias(const char* name, double O::* member,
                                  const char* help = "") {
    TypedParamDesc<O> d = makeParam(name, member, help);
    d.isAlias = true;
    return d;
}

template <typename O, typename N>
inline TypedParamDesc<O> makeAlias(const char* name, double O::* member,
                                  const char* help, N noted) {
    TypedParamDesc<O> d = makeParam(name, member, help, std::move(noted));
    d.isAlias = true;
    return d;
}

/**
 * @brief Temperature key: value in °C (SPICE), stored Kelvin.
 */
template <typename O>
inline TypedParamDesc<O> makeCelsiusParam(const char* name, double O::* member,
                                         const char* help = "") {
    TypedParamDesc<O> d = makeParam(name, member, help);
    d.celsius = true;
    return d;
}

template <typename O, typename N>
inline TypedParamDesc<O> makeCelsiusParam(const char* name, double O::* member,
                                         const char* help, N noted) {
    TypedParamDesc<O> d = makeParam(name, member, help, std::move(noted));
    d.celsius = true;
    return d;
}

/**
 * @brief Canonical parameter bound to an `int` member (LEVEL).
 */
template <typename O>
inline TypedParamDesc<O> makeIntParam(const char* name, int O::* member,
                                     const char* help = "") {
    TypedParamDesc<O> d;
    d.name = name;
    d.help = help ? help : "";
    d.set = [member](O& o, double v) { (o.*member) = static_cast<int>(v); };
    d.get = [member](const O& o) -> double {
        return static_cast<double>((o.*member));
    };
    return d;
}

template <typename O, typename N>
inline TypedParamDesc<O> makeIntParam(const char* name, int O::* member,
                                     const char* help, N noted) {
    TypedParamDesc<O> d = makeIntParam(name, member, help);
    d.noted = std::move(noted);
    return d;
}

/**
 * @brief Generic table fill (cards AND instance lines): aliases first so
 *        the canonical spelling always wins; °C keys get the Kelvin offset;
 *        `noted` hooks fire per applied key.
 * @return keys with no table entry, in input order (typos).
 */
template <class O>
inline std::vector<std::string> applyParamTable(
    O& owner, const std::vector<TypedParamDesc<O>>& table,
    const ParamMap& params) {
    for (int pass = 0; pass < 2; ++pass) {
        for (const auto& d : table) {
            if ((d.isAlias ? 0 : 1) != pass) continue;
            const auto it = params.find(d.name);
            if (it == params.end()) continue;
            double v = it->second;
            if (d.celsius) v += 273.15;
            d.set(owner, v);
            if (d.noted) d.noted(owner);
        }
    }
    std::vector<std::string> unknown;
    for (const auto& kv : params) {
        bool known = false;
        for (const auto& d : table) {
            if (d.name == kv.first) { known = true; break; }
        }
        if (!known) unknown.push_back(kv.first);
    }
    return unknown;
}

/**
 * @brief Human-readable list (canonical first, then aliases).
 */
template <class O>
inline std::string validParamList(const std::vector<TypedParamDesc<O>>& table) {
    std::string canon, alias;
    for (const auto& d : table) {
        std::string& dst = d.isAlias ? alias : canon;
        if (!dst.empty()) dst += " ";
        dst += d.name;
    }
    if (!alias.empty()) canon += "  (aliases: " + alias + ")";
    return canon;
}

} // namespace e_XSpice

#endif // MODEL_PARAMS_H

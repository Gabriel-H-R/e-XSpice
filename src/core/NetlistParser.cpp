/**
 * @author Gabriel Henrique Silva
 * @file NetlistParser.cpp
 * @brief SPICE netlist parser implementation.
 */

#include "NetlistParser.h"
#include "../devices/Resistor.h"
#include "../devices/VoltageSource.h"
#include "../devices/CurrentSource.h"
#include "../devices/Capacitor.h"
#include "../devices/Inductor.h"
#include "../devices/Diode.h"
#include "../devices/SourceWaveform.h"
#include "../devices/Bjt.h"
#include "../devices/Mosfet.h"
#include "../models/core_models/PassiveModels.h"
#include "../models/core_models/MosfetBSIM4Model.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <stdexcept>

namespace e_XSpice {

NetlistParser::NetlistParser()
    : verbose_(false)
    , line_number_(0)
    , next_node_id_(1)  // 0 is ground
{
    // Common ground names
    node_name_to_id_["0"] = 0;
    node_name_to_id_["GND"] = 0;
    node_name_to_id_["gnd"] = 0;

    // Model factories live in the registry itself (built-ins D/R/C/L);
    // only device-letter dispatch is wired here (extensible via
    // registerDeviceCreator for Q/M/...).
    // ---- device-letter dispatch (extensible via registerDeviceCreator) ----
    creators_['R'] = [this](const std::vector<std::string>& t, Circuit& c) {
        return parseResistor(t, c);
    };
    creators_['C'] = [this](const std::vector<std::string>& t, Circuit& c) {
        return parseCapacitor(t, c);
    };
    creators_['L'] = [this](const std::vector<std::string>& t, Circuit& c) {
        return parseInductor(t, c);
    };
    creators_['V'] = [this](const std::vector<std::string>& t, Circuit& c) {
        return parseVoltageSource(t, c);
    };
    creators_['I'] = [this](const std::vector<std::string>& t, Circuit& c) {
        return parseCurrentSource(t, c);
    };
    creators_['D'] = [this](const std::vector<std::string>& t, Circuit& c) {
        return parseDiode(t, c);
    };
    creators_['Q'] = [this](const std::vector<std::string>& t, Circuit& c) {
        return parseBjt(t, c);
    };
    creators_['M'] = [this](const std::vector<std::string>& t, Circuit& c) {
        return parseMosfet(t, c);
    };
}

void NetlistParser::registerDeviceCreator(char letter, DeviceCreator creator) {
    creators_[static_cast<char>(std::toupper(static_cast<unsigned char>(letter)))] =
        std::move(creator);
}

bool NetlistParser::parse(const std::string& filename, Circuit& circuit) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        setError("Could not open file: " + filename);
        return false;
    }
    
    std::string netlist_content;
    std::string line;
    while (std::getline(file, line)) {
        netlist_content += line + "\n";
    }
    file.close();
    
    return parseString(netlist_content, circuit);
}

bool NetlistParser::parseString(const std::string& netlist_str, Circuit& circuit) {
    std::istringstream stream(netlist_str);
    std::string line;
    line_number_ = 0;
    
    if (verbose_) {
        std::cout << "Parsing netlist..." << std::endl;
    }
    
    while (std::getline(stream, line)) {
        line_number_++;
        
        line = cleanLine(line);
        if (line.empty()) {
            continue;  // Empty line
        }
        
        // Comment
        if (line[0] == '*') {
            continue;
        }
        
        // .END command
        if (toUpper(line).find(".END") == 0) {
            break;
        }
        
        // Continuation lines ("+ ...") are appended to the previous card
        if (line[0] == '+') {
            if (continuation_buffer_.empty()) {
                setError("Continuation line '+' without a preceding card");
                return false;
            }
            continuation_buffer_ += " " + line.substr(1);
            continue;
        }

        if (!continuation_buffer_.empty()) {
            if (!parseLine(continuation_buffer_, circuit)) return false;
            continuation_buffer_.clear();
        }
        // The line may itself be the beginning of a next card; keep it in the
        // buffer so a '+' line that follows can extend it.
        continuation_buffer_ = line;
    }
    if (!continuation_buffer_.empty()) {
        if (!parseLine(continuation_buffer_, circuit)) return false;
    }
    
    if (verbose_) {
        std::cout << "Netlist parsed successfully" << std::endl;
        std::cout << "  Nodes: " << circuit.getNumNodes() << std::endl;
        std::cout << "  Devices: " << circuit.getNumDevices() << std::endl;
    }
    
    return true;
}

bool NetlistParser::parseLine(const std::string& line, Circuit& circuit) {
    auto tokens = tokenize(line);
    if (tokens.empty()) {
        return true;
    }

    const char first_char =
        static_cast<char>(std::toupper(static_cast<unsigned char>(tokens[0][0])));

    if (first_char == '.') {
        return parseCommand(tokens);
    }

    // Device-letter dispatch: the parser core never switch-cases physics.
    // Unknown letters (Q, M, ...) fail here until someone registers them.
    const auto it = creators_.find(first_char);
    if (it == creators_.end()) {
        setError("Unknown component type: " + tokens[0]);
        return false;
    }
    return it->second(tokens, circuit);
}

bool NetlistParser::parseResistor(const std::vector<std::string>& tokens, Circuit& circuit) {
    // R<name> <n+> <n-> [model] <value>
    if (tokens.size() < 4) {
        setError("Resistor requires: R<name> <n+> <n-> [model] <value>");
        return false;
    }

    std::string name = tokens[0];
    int node_plus = getNodeId(tokens[1]);
    int node_minus = getNodeId(tokens[2]);

    // Optional shared model card: a token naming a registered R model.
    // (Checked before numeric parsing so `RMOD` never hits parseValue.)
    RCP<const ResistorModel> model;
    size_t value_at = 3;
    if (tokens.size() > 4) {
        RCP<const ResistorModel> candidate =
            models_.findAs<ResistorModel>(tokens[3]);
        if (candidate) {
            model = std::move(candidate);
            value_at = 4;
        }
        // else: plain positional value (checked below); a mistyped model
        // name surfaces as a "cannot parse value" error with its name.
    }
    if (tokens.size() <= value_at) {
        setError("Resistor requires a value");
        return false;
    }

    double value = 0.0;
    try {
        value = parseValue(tokens[value_at]);
    } catch (const std::exception& e) {
        setError(std::string("Resistor value: ") + e.what());
        return false;
    }

    // Trailing card name (`R1 n+ n- <value> <model>`), the Xyce/ngspice
    // spelling; the card-before-value form above is the other accepted one.
    if (!model && tokens.size() > value_at + 1) {
        model = models_.findAs<ResistorModel>(tokens[value_at + 1]);
    }

    if (value <= 0) {
        setError("Resistor value must be positive");
        return false;
    }

    Node* n_plus = getOrCreateNode(circuit, tokens[1], node_plus);
    Node* n_minus = getOrCreateNode(circuit, tokens[2], node_minus);

    try {
        circuit.addDevice(make_rcp<Resistor>(name, n_plus, n_minus, value,
                                             std::move(model)));
    } catch (const std::exception& e) {
        setError(e.what());
        return false;
    }

    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2]
                  << " = " << value << " ohm" << std::endl;
    }

    return true;
}

bool NetlistParser::parseCapacitor(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 4) {
        setError("Capacitor requires: C<name> <n+> <n-> <value> [IC=<v>]");
        return false;
    }
    
    std::string name = tokens[0];
    int node_plus = getNodeId(tokens[1]);
    int node_minus = getNodeId(tokens[2]);
    double value = parseValue(tokens[3]);
    
    if (value <= 0) {
        setError("Capacitor value must be positive");
        return false;
    }
    
    // Initial condition (IC=) — used by the transient initial-condition
    // operating point. (Previously silently ignored.)
    double ic = 0.0;
    bool ic_given = false;
    for (size_t i = 4; i < tokens.size(); ++i) {
        std::string token_upper = toUpper(tokens[i]);
        if (token_upper.find("IC=") == 0) {
            ic = parseValue(tokens[i].substr(3));
            ic_given = true;
        }
    }
    
    Node* n_plus = getOrCreateNode(circuit, tokens[1], node_plus);
    Node* n_minus = getOrCreateNode(circuit, tokens[2], node_minus);
    
    circuit.addDevice(make_rcp<Capacitor>(name, n_plus, n_minus, value,
                                    ic, ic_given));
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2] 
                  << " = " << value << " F" << (ic_given ? " (IC)" : "") << std::endl;
    }
    
    return true;
}

bool NetlistParser::parseInductor(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 4) {
        setError("Inductor requires: L<name> <n+> <n-> <value> [IC=<i>]");
        return false;
    }
    
    std::string name = tokens[0];
    int node_plus = getNodeId(tokens[1]);
    int node_minus = getNodeId(tokens[2]);
    double value = parseValue(tokens[3]);
    
    if (value <= 0) {
        setError("Inductor value must be positive");
        return false;
    }
    
    double ic = 0.0;
    bool ic_given = false;
    for (size_t i = 4; i < tokens.size(); ++i) {
        std::string token_upper = toUpper(tokens[i]);
        if (token_upper.find("IC=") == 0) {
            ic = parseValue(tokens[i].substr(3));
            ic_given = true;
        }
    }
    
    Node* n_plus = getOrCreateNode(circuit, tokens[1], node_plus);
    Node* n_minus = getOrCreateNode(circuit, tokens[2], node_minus);
    
    circuit.addDevice(make_rcp<Inductor>(name, n_plus, n_minus, value,
                                   ic, ic_given));
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2] 
                  << " = " << value << " H" << (ic_given ? " (IC)" : "") << std::endl;
    }
    
    return true;
}

bool NetlistParser::parseVoltageSource(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 5) {
        setError("Voltage source requires: V<name> <n+> <n-> <waveform>");
        return false;
    }
    
    std::string name = tokens[0];
    int node_plus = getNodeId(tokens[1]);
    int node_minus = getNodeId(tokens[2]);
    
    Node* n_plus = getOrCreateNode(circuit, tokens[1], node_plus);
    Node* n_minus = getOrCreateNode(circuit, tokens[2], node_minus);
    
    SourceWaveform waveform;
    if (!parseSourceSpec(tokens, 3, waveform, circuit)) {
        return false;
    }
    
    circuit.addDevice(make_rcp<VoltageSource>(name, n_plus, n_minus, waveform));
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2]
                  << " (" << tokens[3] << ")" << std::endl;
    }
    
    return true;
}

bool NetlistParser::parseCurrentSource(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 5) {
        setError("Current source requires: I<name> <n+> <n-> <waveform>");
        return false;
    }
    
    std::string name = tokens[0];
    int node_plus = getNodeId(tokens[1]);
    int node_minus = getNodeId(tokens[2]);
    
    Node* n_plus = getOrCreateNode(circuit, tokens[1], node_plus);
    Node* n_minus = getOrCreateNode(circuit, tokens[2], node_minus);
    
    SourceWaveform waveform;
    if (!parseSourceSpec(tokens, 3, waveform, circuit)) {
        return false;
    }
    
    circuit.addDevice(make_rcp<CurrentSource>(name, n_plus, n_minus, waveform));
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2]
                  << " (" << tokens[3] << ")" << std::endl;
    }
    
    return true;
}

bool NetlistParser::parseDiode(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 4) {
        setError("Diode requires: D<name> <anode> <cathode> <model> [AREA=<a>] [M=<m>]");
        return false;
    }

    std::string name = tokens[0];
    int node_anode = getNodeId(tokens[1]);
    int node_cathode = getNodeId(tokens[2]);
    const std::string& model_name = tokens[3];

    // The device is agnostic to the card contents: untyped lookup first
    // (for the "missing .MODEL" message), then a checked downcast.
    RCP<const DiodeModel> model = models_.findAs<DiodeModel>(model_name);
    if (!model) {
        if (models_.find(model_name)) {
            setError("Model '" + model_name + "' is not a diode (D) model");
        } else {
            setError("Unknown diode model '" + model_name + "' (missing .MODEL card)");
        }
        return false;
    }

    DiodeInstanceParams inst;
    if (!parseInstanceTail(tokens, 4, inst, DiodeInstanceParams::paramTable(),
                           "Diode")) {
        return false;
    }

    Node* n_anode = getOrCreateNode(circuit, tokens[1], node_anode);
    Node* n_cathode = getOrCreateNode(circuit, tokens[2], node_cathode);

    try {
        circuit.addDevice(make_rcp<Diode>(name, n_anode, n_cathode,
                                          std::move(model), inst));
    } catch (const std::exception& e) {
        setError(e.what());
        return false;
    }

    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " -> " << tokens[2]
                  << " model " << model_name << std::endl;
    }

    return true;
}

bool NetlistParser::parseBjt(const std::vector<std::string>& tokens, Circuit& circuit) {
    // Q<name> <collector> <base> <emitter> [substrate] <model>
    //         [AREA=<a>] [TEMP=<degC>]
    // The 4-node form carries the substrate terminal (needs LEVEL=2).
    if (tokens.size() < 5) {
        setError("Transistor requires: Q<name> <C> <B> <E> [S] <model> "
                 "[AREA=<a>] [TEMP=<c>]");
        return false;
    }

    std::string name = tokens[0];

    // Disambiguate the 3-node vs 4-node (substrate) form: scan for the model
    // as the first token that names a registered card; the node tokens
    // before it are the terminals (3 or 4). When the card is missing (error
    // path) fall back to the classic slots so the message names the model.
    size_t model_at = 0;
    {
        for (size_t i = 1; i < tokens.size(); ++i) {
            if (tokens[i].find('=') != std::string::npos) break;
            if (models_.find(tokens[i])) { model_at = i; break; }
            if (i > 5) break;
        }
        if (model_at == 0) {
            model_at = (tokens.size() >= 6 &&
                        tokens[5].find('=') == std::string::npos) ? 5 : 4;
        }
        const size_t n_nodes = model_at - 1;
        if (n_nodes != 3 && n_nodes != 4) {
            setError("Transistor requires: Q<name> <C> <B> <E> [S] <model> "
                     "[AREA=<a>] [TEMP=<c>]");
            return false;
        }
    }
    const bool has_sub = (model_at - 1 == 4);
    const std::string& model_name = tokens[model_at];

    // The device is agnostic to the card contents: untyped lookup first
    // (for the "missing .MODEL" message), then the checked downcast.
    RCP<const BjtModel> model = models_.findAs<BjtModel>(model_name);
    if (!model) {
        if (models_.find(model_name)) {
            setError("Model '" + model_name + "' is not a BJT model "
                     "(expected .MODEL " + model_name + " NPN|PNP)");
        } else {
            setError("Unknown transistor model '" + model_name +
                     "' (missing .MODEL card)");
        }
        return false;
    }

    BjtInstanceParams inst;
    if (!parseInstanceTail(tokens, model_at + 1, inst,
                           BjtInstanceParams::paramTable(), "Transistor")) {
        return false;
    }

    Node* n_collector = getOrCreateNode(circuit, tokens[1], getNodeId(tokens[1]));
    Node* n_base      = getOrCreateNode(circuit, tokens[2], getNodeId(tokens[2]));
    Node* n_emitter   = getOrCreateNode(circuit, tokens[3], getNodeId(tokens[3]));
    Node* n_sub = nullptr;
    if (has_sub) {
        n_sub = getOrCreateNode(circuit, tokens[4], getNodeId(tokens[4]));
    }

    try {
        circuit.addDevice(make_rcp<Bjt>(name, n_collector, n_base, n_emitter,
                                        std::move(model), inst, n_sub));
    } catch (const std::exception& e) {
        setError(e.what());
        return false;
    }

    if (verbose_) {
        std::cout << "  " << name << ": c=" << tokens[1] << " b=" << tokens[2]
                  << " e=" << tokens[3] << " model " << model_name
                  << (inst.AREA != 1.0 ? " (AREA)" : "") << std::endl;
    }

    return true;
}

bool NetlistParser::parseMosfet(const std::vector<std::string>& tokens,
                                Circuit& circuit) {
    // M<name> <drain> <gate> <source> <bulk> <model>
    //         [W=<w>] [L=<l>] [M=<m>] [AD=<a>] [AS=<a>] [TEMP=<degC>]
    //         [LEVEL=2]  -> BSIM4 model (NMOS2/PMOS2 card)
    if (tokens.size() < 6) {
        setError("MOSFET requires: M<name> <D> <G> <S> <B> <model> "
                 "[W=<w>] [L=<l>] [M=<m>]");
        return false;
    }

    std::string name = tokens[0];
    std::string model_name = tokens[5];

    // Detect LEVEL=2 -> BSIM4 model. The card must be .MODEL X NMOS2/PMOS2.
    bool is_bsim4 = false;
    for (size_t i = 6; i < tokens.size(); ++i) {
        const std::string& t = tokens[i];
        if (t.find("LEVEL=") == 0) {
            try {
                int level = std::stoi(t.substr(6));
                if (level == 2) {
                    is_bsim4 = true;
                    if (model_name == "NMOS") model_name = "NMOS2";
                    else if (model_name == "PMOS") model_name = "PMOS2";
                }
            } catch (...) {
                // fall through: invalid LEVEL value -> error below
            }
            break;
        }
    }

    RCP<const MosfetModel> model = models_.findAs<MosfetModel>(model_name);
    RCP<const MosfetBSIM4Model> bsim_model = models_.findAs<MosfetBSIM4Model>(model_name);
    if (!model && !bsim_model) {
        if (models_.find(model_name)) {
            setError("Model '" + model_name + "' is not a MOSFET model "
                     "(expected .MODEL " + model_name + " NMOS|PMOS|NMOS2|PMOS2)");
        } else {
            setError("Unknown MOSFET model '" + model_name +
                     "' (missing .MODEL card)");
        }
        return false;
    }

    MosfetInstanceParams inst;
    if (!parseInstanceTail(tokens, 6, inst, MosfetInstanceParams::paramTable(),
                           "MOSFET")) {
        return false;
    }

    Node* n_drain  = getOrCreateNode(circuit, tokens[1], getNodeId(tokens[1]));
    Node* n_gate   = getOrCreateNode(circuit, tokens[2], getNodeId(tokens[2]));
    Node* n_source = getOrCreateNode(circuit, tokens[3], getNodeId(tokens[3]));
    Node* n_bulk   = getOrCreateNode(circuit, tokens[4], getNodeId(tokens[4]));

    try {
        if (bsim_model) {
            circuit.addDevice(make_rcp<Mosfet>(name, n_drain, n_gate, n_source,
                                               n_bulk, bsim_model, inst));
        } else {
            circuit.addDevice(make_rcp<Mosfet>(name, n_drain, n_gate, n_source,
                                               n_bulk, std::move(model), inst));
        }
    } catch (const std::exception& e) {
        setError(e.what());
        return false;
    }

    if (verbose_) {
        std::cout << "  " << name << ": d=" << tokens[1] << " g=" << tokens[2]
                  << " s=" << tokens[3] << " b=" << tokens[4]
                  << " model " << model_name
                  << (is_bsim4 ? " [BSIM4]" : "") << std::endl;
    }

    return true;
}

// ===============================================================
// .MODEL + waveform helpers
// ===============================================================
bool NetlistParser::parseModel(const std::vector<std::string>& tokens) {
    // .MODEL <name> <TYPE> ( P1=V1 P2=V2 ... ) — TYPE dispatches through
    // the registry, so the parser never hardcodes device physics.
    if (tokens.size() < 3) {
        setError(".MODEL requires: .MODEL <name> <type> (...)");
        return false;
    }
    const std::string name = tokens[1];
    const std::string type = tokens[2];

    if (!models_.hasModelType(type)) {
        std::string valid;
        for (const auto& t : models_.registeredTypes()) {
            if (!valid.empty()) valid += " ";
            valid += t;
        }
        setError("Unknown .MODEL type '" + type + "' (registered types: " +
                 valid + ")");
        return false;
    }

    // Re-join the remainder and split at '(' ...
    std::string joined;
    for (size_t i = 3; i < tokens.size(); ++i) {
        joined += " " + tokens[i];
    }
    const size_t lp = joined.find('(');
    const size_t rp = joined.rfind(')');
    if (lp == std::string::npos || rp == std::string::npos || rp < lp) {
        setError(".MODEL parameters must be enclosed in parentheses");
        return false;
    }
    std::string body = joined.substr(lp + 1, rp - lp - 1);

    ParamMap params;
    for (const auto& p : tokenize(body)) {
        const size_t eq = p.find('=');
        if (eq == std::string::npos) {
            setError("Malformed .MODEL parameter '" + p + "' (want KEY=value)");
            return false;
        }
        const std::string key = toUpper(p.substr(0, eq));
        try {
            params[key] = parseValue(p.substr(eq + 1));
        } catch (const std::exception& e) {
            setError(std::string(".MODEL parameter '") + key + "': " + e.what());
            return false;
        }
    }

    try {
        ModelCard card;
        card.name = name;
        card.type = type;
        card.params = std::move(params);
        models_.createModel(card);
    } catch (const std::exception& e) {
        setError(e.what());
        return false;
    }

    if (verbose_) {
        std::cout << "  .MODEL " << name << " (" << type << ")" << std::endl;
    }
    return true;
}

bool NetlistParser::parseSourceSpec(const std::vector<std::string>& tokens,
                                    size_t first_param,
                                    SourceWaveform& waveform,
                                    Circuit& circuit) {
    (void)circuit;

    // SPICE gives a source TWO independent specs on the same line:
    //     V1 n+ n- DC 5              -> DC operating point value
    //     V1 n+ n- AC 1 90           -> small-signal phasor for .AC
    //     V1 n+ n- PULSE(0 5 1m) AC 1
    // and "AC mag" alone means DC 0 (the source only excites small signal).
    // The time-domain part is parsed first; an optional trailing "AC mag
    // [phase]" is then attached to the same waveform.
    size_t i = first_param;
    std::string kw = toUpper(tokens[i]);

    if (kw == "AC") {
        // AC-only source: no DC contribution to the operating point.
        waveform = SourceWaveform::dc(0.0);
        ++i;
        return parseACSpec(tokens, i, waveform);
    }

    if (kw == "DC") {
        std::vector<double> p = collectNumbers(tokens, i + 1);
        if (p.empty()) {
            setError("DC source requires a value");
            return false;
        }
        waveform = SourceWaveform::dc(p[0]);
        i = first_param + 2;   // consumed "DC <value>"
        // Anything after may still be the AC spec.
        for (size_t k = i; k < tokens.size(); ++k) {
            if (toUpper(tokens[k]) == "AC") return parseACSpec(tokens, k + 1, waveform);
        }
        return true;
    }

    // Time-domain waveforms: PULSE / SIN / PWL (numbers follow the keyword).
    std::vector<double> p = collectNumbers(tokens, i + 1);
    if (kw == "PULSE") {
        if (p.size() < 2) {
            setError("PULSE requires at least V1 V2");
            return false;
        }
        const double v1 = p[0], v2 = p[1];
        const double td = p.size() > 2 ? p[2] : 0.0;
        const double tr = p.size() > 3 ? p[3] : 0.0;
        const double tf = p.size() > 4 ? p[4] : 0.0;
        const double pw = p.size() > 5 ? p[5] : 0.0;
        const double per = p.size() > 6 ? p[6] : 0.0;
        waveform = SourceWaveform::pulse(v1, v2, td, tr, tf, pw, per);
    } else if (kw == "SIN") {
        if (p.size() < 3) {
            setError("SIN requires VO VA FREQ");
            return false;
        }
        const double v0 = p[0], va = p[1], freq = p[2];
        const double td = p.size() > 3 ? p[3] : 0.0;
        const double theta = p.size() > 4 ? p[4] : 0.0;
        const double phase = p.size() > 5 ? p[5] : 0.0;
        waveform = SourceWaveform::sin(v0, va, freq, td, theta, phase);
    } else if (kw == "PWL") {
        if (p.size() < 2 || (p.size() % 2) != 0) {
            setError("PWL requires (t1 v1 t2 v2 ...)");
            return false;
        }
        std::vector<std::pair<double, double>> pts;
        for (size_t k = 0; k + 1 < p.size(); k += 2) {
            pts.emplace_back(p[k], p[k + 1]);
        }
        waveform = SourceWaveform::pwl(std::move(pts));
    } else {
        setError("Unsupported source specification: " + kw);
        return false;
    }

    // Trailing AC spec on a time-domain source (allowed by SPICE; the AC
    // analysis ignores the waveform and uses only the phasor).
    for (size_t k = i + 1; k < tokens.size(); ++k) {
        if (toUpper(tokens[k]) == "AC") return parseACSpec(tokens, k + 1, waveform);
    }
    return true;
}

bool NetlistParser::parseACSpec(const std::vector<std::string>& tokens,
                                size_t from,
                                SourceWaveform& waveform) {
    const std::vector<double> p = collectNumbers(tokens, from);
    if (p.empty()) {
        setError("AC source specification requires a magnitude");
        return false;
    }
    const double mag   = p[0];
    const double phase = (p.size() > 1) ? p[1] : 0.0;   // degrees, SPICE
    waveform.setAC(mag, phase);
    return true;
}

std::vector<double> NetlistParser::collectNumbers(const std::vector<std::string>& tokens,
                                                  size_t from) {
    std::vector<double> out;
    for (size_t i = from; i < tokens.size(); ++i) {
        std::string t = tokens[i];
        if (t.empty()) continue;
        // Strip a trailing AC/DB/other SPICE qualifiers; only value tokens
        // are numeric.
        const char c = t.front();
        if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+' || c == '.')) {
            continue;
        }
        try {
            out.push_back(parseValue(t));
        } catch (...) {
            // ignore non-numeric token
        }
    }
    return out;
}

bool NetlistParser::parseCommand(const std::vector<std::string>& tokens) {
    std::string cmd = toUpper(tokens[0]);
    
    if (cmd == ".OP") {
        analysis_params_.type = AnalysisType::DC_OP;
        analysis_params_.want_op = true;
        if (verbose_) std::cout << "  Analysis: DC Operating Point" << std::endl;
    } else if (cmd == ".TRAN") {
        if (tokens.size() < 3) {
            setError(".TRAN requires: .TRAN <tstep> <tstop> [tstart] [tmax]");
            return false;
        }
        analysis_params_.type = AnalysisType::TRANSIENT;
        analysis_params_.want_tran = true;
        analysis_params_.tran_tstep = parseValue(tokens[1]);
        analysis_params_.tran_tstop = parseValue(tokens[2]);
        if (tokens.size() > 3) {
            analysis_params_.tran_tstart = parseValue(tokens[3]);
        }
        if (verbose_) {
            std::cout << "  Analysis: Transient (" << analysis_params_.tran_tstep 
                      << "s to " << analysis_params_.tran_tstop << "s)" << std::endl;
        }
    } else if (cmd == ".DC") {
        // Xyce 7.8 §2.1.3:
        //   .DC [LIN|DEC|OCT] <var> <start> <stop> <step|points>
        //          [[LIN|DEC|OCT] <var> <start> <stop> <step|points>]
        // LIN is the default when the keyword is omitted (.DC VIN -10 15 1);
        // for DEC/OCT the last field is points-per-decade/octave. The second
        // (optional) specification is the INNER loop (Xyce allows two levels).
        std::string variation[2];
        std::string source[2];
        double vstart[2] = {0.0, 0.0};
        double vstop[2]  = {0.0, 0.0};
        double vstep[2]  = {0.0, 0.0};
        int level = 0;
        size_t i = 1;
        while (level < 2 && i < tokens.size()) {
            std::string var = "LIN";
            const std::string kw = toUpper(tokens[i]);
            if (kw == "LIN" || kw == "DEC" || kw == "OCT") {
                var = kw;
                ++i;
            }
            if (tokens.size() < i + 4) {
                setError(".DC requires: .DC [LIN|DEC|OCT] <source> "
                         "<start> <stop> <step|points>");
                return false;
            }
            variation[level] = var;
            source[level]    = tokens[i];
            vstart[level]    = parseValue(tokens[i + 1]);
            vstop[level]     = parseValue(tokens[i + 2]);
            vstep[level]     = parseValue(tokens[i + 3]);
            i += 4;
            ++level;
        }
        if (i < tokens.size()) {
            setError(".DC accepts at most two sweep specifications "
                     "(nested sweeps go one level deep)");
            return false;
        }

        analysis_params_.type = AnalysisType::DC_SWEEP;
        analysis_params_.want_dc_sweep = true;
        analysis_params_.dc_variation = variation[0];
        analysis_params_.dc_source    = source[0];
        analysis_params_.dc_start     = vstart[0];
        analysis_params_.dc_stop      = vstop[0];
        analysis_params_.dc_step      = vstep[0];
        analysis_params_.dc_nested      = (level > 1);
        analysis_params_.dc_variation2  = (level > 1) ? variation[1] : "LIN";
        analysis_params_.dc_source2     = (level > 1) ? source[1]    : "";
        analysis_params_.dc_start2      = (level > 1) ? vstart[1]    : 0.0;
        analysis_params_.dc_stop2       = (level > 1) ? vstop[1]     : 0.0;
        analysis_params_.dc_step2       = (level > 1) ? vstep[1]     : 0.0;
        if (verbose_) {
            std::cout << "  Analysis: DC sweep " << variation[0] << " "
                      << source[0] << " from " << vstart[0] << " to "
                      << vstop[0];
            if (level > 1) {
                std::cout << ", nested " << variation[1] << " " << source[1]
                          << " from " << vstart[1] << " to " << vstop[1];
            }
            std::cout << std::endl;
        }
    } else if (cmd == ".TEMP") {
        // .TEMP <celsius> — circuit operating temperature (SPICE: applies to
        // every card that did not pin its own TEMP=).
        if (tokens.size() < 2) {
            setError(".TEMP requires a temperature in Celsius: .TEMP <t>");
            return false;
        }
        double t = 0.0;
        if (!parseTemperature(tokens[1], t)) {
            setError(".TEMP: cannot parse temperature '" + tokens[1] + "'");
            return false;
        }
        if (t < -273.15) {
            setError(".TEMP: below absolute zero");
            return false;
        }
        analysis_params_.temp_given = true;
        analysis_params_.temp_value = t;
        if (verbose_) {
            std::cout << "  Circuit temperature: " << t << " C" << std::endl;
        }
    } else if (cmd == ".OPTIONS" || cmd == ".OPTION") {
        // Solver/sweep options. Unknown tokens warn instead of failing, so a
        // netlist written for another simulator still runs here.
        for (size_t k = 1; k < tokens.size(); ++k) {
            const std::string opt = toUpper(tokens[k]);
            const std::string body = opt.substr(0, opt.find('='));
            const std::string value = (opt.find('=') != std::string::npos)
                                    ? opt.substr(opt.find('=') + 1) : "";
            if (body == "NOFAIL") {
                // Keep the sweep going past a non-converged point (the row
                // is flagged and holds the last converged solution).
                analysis_params_.opt_nofail = true;
            } else if (body == "RETRY" || body == "MAXRETRY") {
                if (value.empty()) {
                    setError(".OPTIONS RETRY=<n> requires a value");
                    return false;
                }
                analysis_params_.opt_max_retries = std::max(0, std::stoi(value));
            } else if (body == "MAXITERS" || body == "ITL4") {
                if (value.empty()) {
                    setError(".OPTIONS MAXITERS=<n> requires a value");
                    return false;
                }
                analysis_params_.opt_max_iterations = std::max(1, std::stoi(value));
            } else if (body == "TEMP") {
                // SPICE3 spelling of `.TEMP <t>`.
                double t = 0.0;
                if (value.empty() || !parseTemperature(value, t) || t < -273.15) {
                    setError(".OPTIONS TEMP=<celsius> has a bad value");
                    return false;
                }
                analysis_params_.temp_given = true;
                analysis_params_.temp_value = t;
            } else if (body == "SWEEPPAR") {
                // Distribute independent sweep points across MPI ranks
                // (each rank keeps its own continuation chain).
                analysis_params_.opt_sweep_parallel = true;
            } else if (body == "ACFREQPAR") {
                // Distribute the AC FREQUENCY sweep across ranks. Unlike the
                // DC sweep, frequencies are independent and each rank
                // factorizes the same replicated linearization, so the answer
                // is bit-identical to a serial run.
                analysis_params_.opt_ac_freq_parallel = true;
            } else {
                analysis_params_.option_warnings.push_back(tokens[k]);
            }
        }
        if (verbose_ && !analysis_params_.option_warnings.empty()) {
            std::cout << "  .OPTIONS: ignoring unsupported option(s):";
            for (const auto& w : analysis_params_.option_warnings) {
                std::cout << " " << w;
            }
            std::cout << std::endl;
        }
    } else if (cmd == ".AC") {
        if (tokens.size() < 5) {
            setError(".AC requires: .AC <DEC|OCT|LIN> <points> <fstart> <fstop>");
            return false;
        }
        analysis_params_.type = AnalysisType::AC;
        analysis_params_.want_ac = true;
        analysis_params_.ac_variation = toUpper(tokens[1]);
        if (analysis_params_.ac_variation != "DEC" &&
            analysis_params_.ac_variation != "OCT" &&
            analysis_params_.ac_variation != "LIN") {
            setError(".AC variation must be LIN, DEC or OCT (got '" +
                     tokens[1] + "')");
            return false;
        }
        analysis_params_.ac_points = std::stoi(tokens[2]);
        analysis_params_.ac_fstart = parseValue(tokens[3]);
        analysis_params_.ac_fstop = parseValue(tokens[4]);
        if (analysis_params_.ac_points < 1) {
            setError(".AC needs at least one point per decade/octave");
            return false;
        }
        if (!(analysis_params_.ac_fstart > 0.0) ||
            !(analysis_params_.ac_fstop > 0.0)) {
            setError(".AC frequencies must be positive");
            return false;
        }
        if (verbose_) {
            std::cout << "  Analysis: AC " << analysis_params_.ac_variation 
                      << " sweep" << std::endl;
        }
    } else if (cmd == ".NOISE") {
        //   .NOISE V(out[,ref]) <src> (DEC|OCT|LIN) <pts> <fstart> <fstop>
        //         [pts_per_summary]
        // ngspice/Xyce/SPICE3 order. `src` may also be written SRC=V1.
        if (tokens.size() < 6) {
            setError(".NOISE requires: .NOISE V(out[,ref]) <src> "
                     "<DEC|OCT|LIN> <points> <fstart> <fstop> [pts_per_summary]");
            return false;
        }
        // ---- V(out[,ref]) ----
        // The tokenizer emits parentheses as separate tokens, so V(out,ref)
        // arrives as V ( out,ref ) and the name is rebuilt from the tokens in
        // between (a space after the comma is tolerated).
        if (tokens.size() < 3 || toUpper(tokens[1]) != "V" || tokens[2] != "(") {
            setError(".NOISE first argument must be V(node) or V(node,ref)");
            return false;
        }
        size_t close = 3;
        while (close < tokens.size() && tokens[close] != ")") ++close;
        if (close >= tokens.size() || close == 3) {
            setError(".NOISE: malformed V(...) argument");
            return false;
        }
        std::string spec;
        for (size_t k = 3; k < close; ++k) spec += tokens[k];
        const size_t comma = spec.find(',');
        if (comma == std::string::npos) {
            analysis_params_.noise_out_node = spec;
            analysis_params_.noise_ref_node = "0";
        } else {
            analysis_params_.noise_out_node = spec.substr(0, comma);
            analysis_params_.noise_ref_node = spec.substr(comma + 1);
        }
        if (analysis_params_.noise_out_node.empty()) {
            setError(".NOISE: empty output node in V(...)");
            return false;
        }

        // ---- input source (optional, ngspice style) / SRC= (Xyce style) ----
        size_t idx = close + 1;
        std::string src_token = toUpper(tokens[idx]);
        if (src_token.rfind("SRC=", 0) == 0) {
            analysis_params_.noise_source = tokens[idx].substr(4);
            analysis_params_.noise_source_given = true;
            ++idx;
        } else if (src_token != "DEC" && src_token != "OCT" &&
                   src_token != "LIN") {
            analysis_params_.noise_source = tokens[idx];
            analysis_params_.noise_source_given = true;
            ++idx;
        }

        // ---- sweep ----
        if (tokens.size() < idx + 4) {
            setError(".NOISE requires <DEC|OCT|LIN> <points> <fstart> <fstop>");
            return false;
        }
        const std::string variation = toUpper(tokens[idx]);
        if (variation != "DEC" && variation != "OCT" && variation != "LIN") {
            setError(".NOISE variation must be LIN, DEC or OCT (got '" +
                     tokens[idx] + "')");
            return false;
        }
        analysis_params_.type = AnalysisType::NOISE;
        analysis_params_.want_noise = true;
        analysis_params_.ac_variation = variation;
        analysis_params_.ac_points = std::stoi(tokens[idx + 1]);
        analysis_params_.ac_fstart = parseValue(tokens[idx + 2]);
        analysis_params_.ac_fstop = parseValue(tokens[idx + 3]);
        if (analysis_params_.ac_points < 1) {
            setError(".NOISE needs at least one point per decade/octave");
            return false;
        }
        if (!(analysis_params_.ac_fstart > 0.0) ||
            !(analysis_params_.ac_fstop > 0.0)) {
            setError(".NOISE frequencies must be positive");
            return false;
        }
        if (analysis_params_.ac_fstop < analysis_params_.ac_fstart) {
            setError(".NOISE: stop frequency is below the start frequency");
            return false;
        }
        // ---- pts_per_summary: per-generator contribution breakdown ----
        if (tokens.size() > idx + 4) {
            const int summary = std::stoi(tokens[idx + 4]);
            analysis_params_.noise_detail = (summary > 0);
        }
        if (verbose_) {
            std::cout << "  Analysis: NOISE V("
                      << analysis_params_.noise_out_node;
            if (analysis_params_.noise_ref_node != "0") {
                std::cout << "," << analysis_params_.noise_ref_node;
            }
            std::cout << ") " << variation << " sweep, input source "
                      << (analysis_params_.noise_source.empty()
                              ? "(first V source)" : analysis_params_.noise_source)
                      << std::endl;
        }
    } else if (cmd == ".TF") {
        if (tokens.size() < 3) {
            setError(".TF requires: .TF <outvar> <insrc>");
            return false;
        }
        // tokens[1] = V ( ( ... )  or  I ( ( ... )  — reassemble the
        // parenthesized argument from the tokenized form.
        if (tokens.size() < 4 || tokens[2] != "(") {
            setError(".TF first argument must be V(...) or I(...)");
            return false;
        }
        size_t close = 3;
        while (close < tokens.size() && tokens[close] != ")") ++close;
        if (close >= tokens.size() || close == 3) {
            setError(".TF: malformed output variable");
            return false;
        }
        std::string spec;
        for (size_t k = 1; k <= close; ++k) spec += tokens[k];
        analysis_params_.tf_out_var = spec;
        analysis_params_.tf_input_source = tokens[close + 1];
        analysis_params_.type = AnalysisType::TF;
        analysis_params_.want_tf = true;
        if (verbose_) {
            std::cout << "  Analysis: TF " << analysis_params_.tf_out_var
                      << " / " << analysis_params_.tf_input_source << std::endl;
        }
    } else if (cmd == ".PZ") {
        if (tokens.size() < 7) {
            setError(".PZ requires: .PZ node1 node2 node3 node4 (vol|cur) (pol|zer|pz)");
            return false;
        }
        analysis_params_.pz_out1 = tokens[1];
        analysis_params_.pz_out2 = tokens[2];
        analysis_params_.pz_in1 = tokens[3];
        analysis_params_.pz_in2 = tokens[4];
        const std::string mode = toUpper(tokens[5]);
        if (mode == "CUR") {
            analysis_params_.pz_current_input = true;
        } else if (mode == "VOL") {
            analysis_params_.pz_current_input = false;
        } else {
            setError(".PZ mode must be vol or cur (got '" + tokens[5] + "')");
            return false;
        }
        const std::string kind = toUpper(tokens[6]);
        if (kind == "POL") {
            analysis_params_.pz_poles_only = true;
            analysis_params_.pz_zeros_only = false;
        } else if (kind == "ZER") {
            analysis_params_.pz_poles_only = false;
            analysis_params_.pz_zeros_only = true;
        } else if (kind == "PZ") {
            analysis_params_.pz_poles_only = false;
            analysis_params_.pz_zeros_only = false;
        } else {
            setError(".PZ kind must be pol, zer or pz (got '" + tokens[6] + "')");
            return false;
        }
        analysis_params_.type = AnalysisType::PZ;
        analysis_params_.want_pz = true;
        if (verbose_) {
            std::cout << "  Analysis: PZ V(" << analysis_params_.pz_out1
                      << "," << analysis_params_.pz_out2 << ") / "
                      << (analysis_params_.pz_current_input ? "I(" : "V(")
                      << analysis_params_.pz_in1 << "," << analysis_params_.pz_in2
                      << ") " << kind << std::endl;
        }
    } else if (cmd == ".SENS") {
        // .SENS V(out) [param1, param2, ...]        (ngspice style)
        // .SENS objfunc=V(out) param=R1,R2            (Xyce style)
        // .SENS AC V(out) R1,R2                       (AC sensitivity)
        // .SENS TRAN V(out) R1,R2                     (transient sensitivity)
        if (tokens.size() < 2) {
            setError(".SENS requires an output variable");
            return false;
        }
        size_t idx = 1;
        std::string outvar;
        // AC / transient sensitivity: .SENS AC V(out) R1,R2
        bool ac_sens = false;
        bool tran_sens = false;
        std::string tok = toUpper(tokens[idx]);
        if (tok == "AC") {
            ac_sens = true;
            analysis_params_.sens_ac = true;
            ++idx;
            tok = toUpper(tokens[idx]);
        } else if (tok == "TRAN" || tok == "TR") {
            tran_sens = true;
            analysis_params_.sens_tran = true;
            ++idx;
            tok = toUpper(tokens[idx]);
        }
        // Xyce style: objfunc=... or acobjfunc=... (and the transient
        // counterpart tranobjfunc=...).
        //
        // The tokenizer splits the parentheses, so `acobjfunc=V(out)` arrives
        // as "ACOBJFUNC=V" "(" "out" ")": the parenthesized part has to be
        // reassembled (same as the ngspice spelling below).
        const auto reassemble = [&](size_t& i, const std::string& head) {
            std::string spec = head;
            if (i + 1 < tokens.size() && tokens[i + 1] == "(") {
                size_t close = i + 2;
                while (close < tokens.size() && tokens[close] != ")") ++close;
                for (size_t k = i + 1; k <= close && k < tokens.size(); ++k) {
                    spec += tokens[k];
                }
                i = close + 1;
            } else {
                ++i;
            }
            return spec;
        };

        bool xyce_style = false;
        if (tok.rfind("OBJFUNC=", 0) == 0) {
            xyce_style = true;
            size_t i = idx;
            outvar = reassemble(i, tokens[idx].substr(8));
            idx = i;
        } else if (tok.rfind("ACOBJFUNC=", 0) == 0) {
            xyce_style = true;
            ac_sens = true;
            analysis_params_.sens_ac = true;
            size_t i = idx;
            outvar = reassemble(i, tokens[idx].substr(10));
            idx = i;
        } else if (tok.rfind("TRANOBJFUNC=", 0) == 0) {
            xyce_style = true;
            tran_sens = true;
            analysis_params_.sens_tran = true;
            size_t i = idx;
            outvar = reassemble(i, tokens[idx].substr(12));
            idx = i;
        } else {
            // ngspice style: first token is the output variable (V(...)/I(...)).
            // Reassemble parenthesized form.
            if (tokens[idx] != "V" && toUpper(tokens[idx]) != "I" &&
                tokens[idx] != "v" && tokens[idx] != "i") {
                // maybe it's V(...) already tokenized
            }
            // Handle tokenized V ( node ) form.
            if (tokens.size() > idx + 1 && tokens[idx + 1] == "(") {
                size_t close = idx + 2;
                while (close < tokens.size() && tokens[close] != ")") ++close;
                std::string spec;
                for (size_t k = idx; k <= close; ++k) spec += tokens[k];
                outvar = spec;
                idx = close + 1;
            } else {
                outvar = tokens[idx];
                ++idx;
            }
        }
        analysis_params_.sens_output_var = outvar;

        // Parameter list.
        std::vector<std::string> param_tokens;
        if (xyce_style) {
            // param=...
            while (idx < tokens.size()) {
                std::string t = tokens[idx];
                if (toUpper(t).rfind("PARAM=", 0) == 0) {
                    std::string val = t.substr(6);
                    // split on comma
                    size_t start = 0;
                    for (size_t i = 0; i <= val.size(); ++i) {
                        if (i == val.size() || val[i] == ',') {
                            if (i > start) param_tokens.push_back(val.substr(start, i - start));
                            start = i + 1;
                        }
                    }
                }
                ++idx;
            }
        } else {
            // ngspice style: remaining tokens are params (comma-separated or space).
            while (idx < tokens.size()) {
                std::string t = tokens[idx];
                size_t start = 0;
                for (size_t i = 0; i <= t.size(); ++i) {
                    if (i == t.size() || t[i] == ',') {
                        if (i > start) param_tokens.push_back(t.substr(start, i - start));
                        start = i + 1;
                    }
                }
                ++idx;
            }
        }

        // Parse each param: "R1" or "R1:R"
        for (const auto& pt : param_tokens) {
            size_t colon = pt.find(':');
            if (colon == std::string::npos) {
                analysis_params_.sens_params.emplace_back(pt, "");
            } else {
                analysis_params_.sens_params.emplace_back(pt.substr(0, colon),
                                                          pt.substr(colon + 1));
            }
        }

        analysis_params_.type = AnalysisType::SENSITIVITY;
        analysis_params_.want_sens = true;
        if (verbose_) {
            std::cout << "  Analysis: SENS "
                      << (ac_sens ? "AC " : (tran_sens ? "TRAN " : ""))
                      << analysis_params_.sens_output_var
                      << " (" << param_tokens.size() << " param token(s))"
                      << std::endl;
        }
    } else if (cmd == ".MODEL") {
        return parseModel(tokens);
    } else if (cmd == ".END") {
        // handled in parseString
    } else {
        if (verbose_) {
            std::cout << "  Ignoring command: " << cmd << std::endl;
        }
    }
    
    return true;
}

// ===============================================================
// Utilities
// ===============================================================

Node* NetlistParser::getOrCreateNode(Circuit& circuit, const std::string& name, int id) {
    Node* n = circuit.getNode(id);
    if (!n) {
        n = circuit.addNode(id, name);
    }
    return n;
}

int NetlistParser::getNodeId(const std::string& node_name) {
    std::string name = node_name;
    
    // Ground aliases
    if (name == "0" || toUpper(name) == "GND") {
        return 0;
    }
    
    // Check if it already exists
    auto it = node_name_to_id_.find(name);
    if (it != node_name_to_id_.end()) {
        return it->second;
    }
    
    // Create new ID
    int id = next_node_id_++;
    node_name_to_id_[name] = id;
    return id;
}

std::vector<std::string> NetlistParser::tokenize(const std::string& line) {
    std::vector<std::string> tokens;
    std::string spaced;
    spaced.reserve(line.size() + 8);

    // Parentheses become standalone tokens so ".MODEL D (IS=1e-12)" and
    // "PULSE(0 5 1m)" parse with the same tokenizer.
    for (char c : line) {
        if (c == '(') { spaced += " ( "; }
        else if (c == ')') { spaced += " ) "; }
        else { spaced.push_back(c); }
    }

    std::istringstream stream(spaced);
    std::string token;
    while (stream >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

bool NetlistParser::parseTemperature(const std::string& token, double& celsius) {
    // SPICE convention: `.TEMP 27` and `.TEMP 27C` are the same; a trailing
    // `K` marks an explicit Kelvin value (no offset applied). parseValue()
    // is NOT used here: its `K` suffix means kilo (1e3).
    std::string s = token;
    bool kelvin = false;
    if (!s.empty()) {
        const char last = static_cast<char>(std::toupper(
            static_cast<unsigned char>(s.back())));
        if (last == 'C') {
            s.pop_back();
        } else if (last == 'K') {
            s.pop_back();
            kelvin = true;
        }
    }
    if (s.empty()) return false;
    try {
        size_t idx = 0;
        const double v = std::stod(s, &idx);
        if (idx != s.size()) return false;
        celsius = kelvin ? (v - Const::CtoK) : v;
    } catch (const std::exception&) {
        return false;
    }
    return true;
}

double NetlistParser::parseValue(const std::string& str) {
    // Supports multipliers: T, G, MEG, k, m, u, n, p, f
    if (str.empty()) {
        throw std::invalid_argument("empty value");
    }
    std::string s = str;
    double multiplier = 1.0;
    
    if (s.length() > 1) {
        std::string upper = toUpper(s);
        
        if (upper.find("MEG") != std::string::npos) {
            multiplier = 1e6;
            s = s.substr(0, s.length() - 3);
        } else {
            char last = std::toupper(s.back());
            switch (last) {
                case 'T': multiplier = 1e12; s.pop_back(); break;
                case 'G': multiplier = 1e9; s.pop_back(); break;
                case 'K': multiplier = 1e3; s.pop_back(); break;
                case 'M': multiplier = 1e-3; s.pop_back(); break;
                case 'U': multiplier = 1e-6; s.pop_back(); break;
                case 'N': multiplier = 1e-9; s.pop_back(); break;
                case 'P': multiplier = 1e-12; s.pop_back(); break;
                case 'F': multiplier = 1e-15; s.pop_back(); break;
            }
        }
    }
    
    size_t idx = 0;
    double value = 0.0;
    try {
        value = std::stod(s, &idx);
    } catch (const std::exception&) {
        throw std::invalid_argument("cannot parse value '" + str + "'");
    }
    if (idx != s.size()) {
        throw std::invalid_argument("cannot parse value '" + str + "'");
    }
    return value * multiplier;
}

std::string NetlistParser::cleanLine(const std::string& line) {
    std::string cleaned;
    
    // Remove inline comments (after $ or ;)
    size_t comment_pos = line.find_first_of(";$");
    if (comment_pos != std::string::npos) {
        cleaned = line.substr(0, comment_pos);
    } else {
        cleaned = line;
    }
    
    // Remove leading and trailing spaces
    size_t start = cleaned.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    
    size_t end = cleaned.find_last_not_of(" \t\r\n");
    return cleaned.substr(start, end - start + 1);
}

std::string NetlistParser::toUpper(const std::string& str) {
    std::string upper = str;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    return upper;
}

void NetlistParser::setError(const std::string& msg) {
    error_message_ = "Line " + std::to_string(line_number_) + ": " + msg;
}

} // namespace e_XSpice

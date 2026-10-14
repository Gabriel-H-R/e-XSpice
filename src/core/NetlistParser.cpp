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

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>

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
    
    char first_char = std::toupper(tokens[0][0]);
    
    switch (first_char) {
        case 'R':
            return parseResistor(tokens, circuit);
        case 'C':
            return parseCapacitor(tokens, circuit);
        case 'L':
            return parseInductor(tokens, circuit);
        case 'V':
            return parseVoltageSource(tokens, circuit);
        case 'I':
            return parseCurrentSource(tokens, circuit);
        case 'D':
            return parseDiode(tokens, circuit);
        case '.':
            return parseCommand(tokens);
        default:
            setError("Unknown component type: " + tokens[0]);
            return false;
    }
}

bool NetlistParser::parseResistor(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 4) {
        setError("Resistor requires: R<name> <n+> <n-> <value>");
        return false;
    }
    
    std::string name = tokens[0];
    int node_plus = getNodeId(tokens[1]);
    int node_minus = getNodeId(tokens[2]);
    double value = parseValue(tokens[3]);
    
    if (value <= 0) {
        setError("Resistor value must be positive");
        return false;
    }
    
    Node* n_plus = getOrCreateNode(circuit, tokens[1], node_plus);
    Node* n_minus = getOrCreateNode(circuit, tokens[2], node_minus);
    
    circuit.addDevice(new Resistor(name, n_plus, n_minus, value));
    
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
    
    circuit.addDevice(new Capacitor(name, n_plus, n_minus, value,
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
    
    circuit.addDevice(new Inductor(name, n_plus, n_minus, value,
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
    
    circuit.addDevice(new VoltageSource(name, n_plus, n_minus, waveform));
    
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
    
    circuit.addDevice(new CurrentSource(name, n_plus, n_minus, waveform));
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2]
                  << " (" << tokens[3] << ")" << std::endl;
    }
    
    return true;
}

bool NetlistParser::parseDiode(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 4) {
        setError("Diode requires: D<name> <anode> <cathode> <model> [AREA=<a>]");
        return false;
    }
    
    std::string name = tokens[0];
    int node_anode = getNodeId(tokens[1]);
    int node_cathode = getNodeId(tokens[2]);
    const std::string& model_name = tokens[3];
    
    auto it = diode_models_.find(toUpper(model_name));
    if (it == diode_models_.end()) {
        setError("Unknown diode model '" + model_name + "' (missing .MODEL card)");
        return false;
    }
    
    DiodeModel model = it->second;
    for (size_t i = 4; i < tokens.size(); ++i) {
        const std::string t = toUpper(tokens[i]);
        if (t.find("AREA=") == 0) {
            model.AREA = parseValue(tokens[i].substr(5));
        }
    }
    
    Node* n_anode = getOrCreateNode(circuit, tokens[1], node_anode);
    Node* n_cathode = getOrCreateNode(circuit, tokens[2], node_cathode);
    
    circuit.addDevice(new Diode(name, n_anode, n_cathode, model));
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " -> " << tokens[2]
                  << " model " << model_name << std::endl;
    }
    
    return true;
}

// ===============================================================
// .MODEL + waveform helpers
// ===============================================================

bool NetlistParser::parseModel(const std::vector<std::string>& tokens) {
    // .MODEL <name> D ( P1=V1 P2=V2 ... )
    if (tokens.size() < 3) {
        setError(".MODEL requires: .MODEL <name> <type> (...)");
        return false;
    }
    const std::string name = toUpper(tokens[1]);
    const std::string type = toUpper(tokens[2]);
    
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
    
    std::vector<std::string> params = tokenize(body);
    if (type == "D") {
        DiodeModel m;
        parseDiodeParams(params, m);
        diode_models_[name] = m;
    } else {
        if (verbose_) {
            std::cout << "  Ignoring .MODEL of type " << type << std::endl;
        }
    }
    
    if (verbose_) {
        std::cout << "  .MODEL " << name << " (" << type << ")" << std::endl;
    }
    return true;
}

void NetlistParser::parseDiodeParams(const std::vector<std::string>& params,
                                     DiodeModel& m) {
    for (const auto& p : params) {
        const size_t eq = p.find('=');
        if (eq == std::string::npos) continue;
        const std::string key = toUpper(p.substr(0, eq));
        const std::string val = p.substr(eq + 1);
        const double v = parseValue(val);
        
        if      (key == "IS" || key == "JS") m.IS = v;
        else if (key == "N")                m.N = v;
        else if (key == "RS")               m.RS = v;
        else if (key == "CJO" || key == "CJ" || key == "CJ0") m.CJO = v;
        else if (key == "VJ" || key == "PB") m.VJ = v;
        else if (key == "M")                m.M = v;
        else if (key == "FC")               m.FC = v;
        else if (key == "BV" || key == "VB") { m.BV = v; m.BVGiven = true; }
        else if (key == "IBV")              m.IBV = v;
        else if (key == "NBV")              m.NBV = v;
        else if (key == "TT")               m.TT = v;
        else if (key == "EG")               m.EG = v;
        else if (key == "XTI")              m.XTI = v;
        else if (key == "IKF")              m.IKF = v;
        else if (key == "ISR")              m.ISR = v;
        else if (key == "NR")               m.NR = v;
        else if (key == "TNOM")             m.TNOM = v + 273.15;
        else if (key == "TEMP")             m.TEMP = v + 273.15;
    }
}

bool NetlistParser::parseSourceSpec(const std::vector<std::string>& tokens,
                                    size_t first_param,
                                    SourceWaveform& waveform,
                                    Circuit& circuit) {
    (void)circuit;
    const std::string kw = toUpper(tokens[first_param]);
    std::vector<double> p = collectNumbers(tokens, first_param + 1);

    if (kw == "DC" || kw == "AC") {
        if (p.empty()) {
            setError("DC/AC source requires a value");
            return false;
        }
        waveform = SourceWaveform::dc(p[0]);
        return true;
    }
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
        return true;
    }
    if (kw == "SIN") {
        if (p.size() < 3) {
            setError("SIN requires VO VA FREQ");
            return false;
        }
        const double v0 = p[0], va = p[1], freq = p[2];
        const double td = p.size() > 3 ? p[3] : 0.0;
        const double theta = p.size() > 4 ? p[4] : 0.0;
        const double phase = p.size() > 5 ? p[5] : 0.0;
        waveform = SourceWaveform::sin(v0, va, freq, td, theta, phase);
        return true;
    }
    if (kw == "PWL") {
        if (p.size() < 2 || (p.size() % 2) != 0) {
            setError("PWL requires (t1 v1 t2 v2 ...)");
            return false;
        }
        std::vector<std::pair<double,double>> pts;
        for (size_t i = 0; i + 1 < p.size(); i += 2) {
            pts.emplace_back(p[i], p[i + 1]);
        }
        waveform = SourceWaveform::pwl(std::move(pts));
        return true;
    }
    
    setError("Unsupported source specification: " + kw);
    return false;
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
        if (verbose_) std::cout << "  Analysis: DC Operating Point" << std::endl;
    } else if (cmd == ".TRAN") {
        if (tokens.size() < 3) {
            setError(".TRAN requires: .TRAN <tstep> <tstop> [tstart] [tmax]");
            return false;
        }
        analysis_params_.type = AnalysisType::TRANSIENT;
        analysis_params_.tran_tstep = parseValue(tokens[1]);
        analysis_params_.tran_tstop = parseValue(tokens[2]);
        if (tokens.size() > 3) {
            analysis_params_.tran_tstart = parseValue(tokens[3]);
        }
        if (verbose_) {
            std::cout << "  Analysis: Transient (" << analysis_params_.tran_tstep 
                      << "s to " << analysis_params_.tran_tstop << "s)" << std::endl;
        }
    } else if (cmd == ".AC") {
        if (tokens.size() < 5) {
            setError(".AC requires: .AC <DEC|OCT|LIN> <points> <fstart> <fstop>");
            return false;
        }
        analysis_params_.type = AnalysisType::AC;
        analysis_params_.ac_variation = toUpper(tokens[1]);
        analysis_params_.ac_points = std::stoi(tokens[2]);
        analysis_params_.ac_fstart = parseValue(tokens[3]);
        analysis_params_.ac_fstop = parseValue(tokens[4]);
        if (verbose_) {
            std::cout << "  Analysis: AC " << analysis_params_.ac_variation 
                      << " sweep" << std::endl;
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

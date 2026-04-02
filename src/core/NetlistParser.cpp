/**
 * @author Gabriel Henrique Silva
 * @file NetlistParser.cpp
 * @brief SPICE netlist parser implementation
 */

#include "NetlistParser.h"
#include "../devices/Resistor.h"
#include "../devices/VoltageSource.h"
#include "../devices/CurrentSource.h"
#include "../devices/Capacitor.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>

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
        
        if (!parseLine(line, circuit)) {
            return false;
        }
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
    
    // Create nodes if they do not exist
    Node* n_plus = circuit.getNode(node_plus);
    if (!n_plus) {
        n_plus = circuit.addNode(node_plus, tokens[1]);
    }
    
    Node* n_minus = circuit.getNode(node_minus);
    if (!n_minus) {
        n_minus = circuit.addNode(node_minus, tokens[2]);
    }
    
    Resistor* r = new Resistor(name, n_plus, n_minus, value);
    circuit.addDevice(r);
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2] 
                 << " = " << value << " Ω" << std::endl;
    }
    
    return true;
}

bool NetlistParser::parseCapacitor(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 4) {
        setError("Capacitor requires: C<name> <n+> <n-> <value>");
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
    
    // Create nodes if they do not exist
    Node* n_plus = circuit.getNode(node_plus);
    if (!n_plus) {
        n_plus = circuit.addNode(node_plus, tokens[1]);
    }
    
    Node* n_minus = circuit.getNode(node_minus);
    if (!n_minus) {
        n_minus = circuit.addNode(node_minus, tokens[2]);
    }
    
    Capacitor* c = new Capacitor(name, n_plus, n_minus, value);
    
    // Initial condition handling (Ex: IC=0)
    for (size_t i = 4; i < tokens.size(); ++i) {
        std::string token_upper = toUpper(tokens[i]);
        if (token_upper.find("IC=") == 0) {
            // Extract value after "IC=" using its own parseValue
            double ic = parseValue(tokens[i].substr(3));
            
            // Assuming your Capacitor class has a method to set initial voltage.
            // If the method name is different, adjust here:
            // c->setInitialVoltage(ic); 
        }
    }
    
    circuit.addDevice(c);
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2] 
                  << " = " << value << " F" << std::endl;
    }
    
    return true;
}

/*
bool NetlistParser::parseCapacitor(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 4) {
        setError("Capacitor requires: C<name> <n+> <n-> <value>");
        return false;
    }
    
    // Similar implementation to resistor
    // For now, only a placeholder (will be implemented later)
    if (verbose_) {
        std::cout << "  " << tokens[0] << " (Capacitor) - not yet implemented" << std::endl;
    }
    
    return true;  // Temporarily ignored
}
*/

bool NetlistParser::parseInductor(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 4) {
        setError("Inductor requires: L<name> <n+> <n-> <value>");
        return false;
    }
    
    // Placeholder
    if (verbose_) {
        std::cout << "  " << tokens[0] << " (Inductor) - not yet implemented" << std::endl;
    }
    
    return true;  // Temporarily ignored
}

bool NetlistParser::parseVoltageSource(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 5) {
        setError("Voltage source requires: V<name> <n+> <n-> DC <value>");
        return false;
    }
    
    std::string name = tokens[0];
    int node_plus = getNodeId(tokens[1]);
    int node_minus = getNodeId(tokens[2]);
    
    // Check source type
    std::string type = toUpper(tokens[3]);
    if (type != "DC" && type != "AC") {
        setError("Only DC and AC voltage sources supported");
        return false;
    }
    
    double value = parseValue(tokens[4]);
    
    // Create nodes
    Node* n_plus = circuit.getNode(node_plus);
    if (!n_plus) {
        n_plus = circuit.addNode(node_plus, tokens[1]);
    }
    
    Node* n_minus = circuit.getNode(node_minus);
    if (!n_minus) {
        n_minus = circuit.addNode(node_minus, tokens[2]);
    }
    
    VoltageSource* v = new VoltageSource(name, n_plus, n_minus, value);
    circuit.addDevice(v);
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2] 
                 << " = " << value << " V" << std::endl;
    }
    
    return true;
}

bool NetlistParser::parseCurrentSource(const std::vector<std::string>& tokens, Circuit& circuit) {
    if (tokens.size() < 5) {
        setError("Current source requires: I<name> <n+> <n-> DC <value>");
        return false;
    }
    
    std::string name = tokens[0];
    int node_plus = getNodeId(tokens[1]);
    int node_minus = getNodeId(tokens[2]);
    
    std::string type = toUpper(tokens[3]);
    if (type != "DC") {
        setError("Only DC current sources supported");
        return false;
    }
    
    double value = parseValue(tokens[4]);
    
    // Create nodes
    Node* n_plus = circuit.getNode(node_plus);
    if (!n_plus) {
        n_plus = circuit.addNode(node_plus, tokens[1]);
    }
    
    Node* n_minus = circuit.getNode(node_minus);
    if (!n_minus) {
        n_minus = circuit.addNode(node_minus, tokens[2]);
    }
    
    CurrentSource* i = new CurrentSource(name, n_plus, n_minus, value);
    circuit.addDevice(i);
    
    if (verbose_) {
        std::cout << "  " << name << ": " << tokens[1] << " - " << tokens[2] 
                 << " = " << value << " A" << std::endl;
    }
    
    return true;
}

bool NetlistParser::parseCommand(const std::vector<std::string>& tokens) {
    std::string cmd = toUpper(tokens[0]);
    
    if (cmd == ".OP") {
        analysis_params_.type = AnalysisType::DC_OP;
        if (verbose_) {
            std::cout << "  Analysis: DC Operating Point" << std::endl;
        }
    } else if (cmd == ".TRAN") {
        if (tokens.size() < 3) {
            setError(".TRAN requires: .TRAN <tstep> <tstop> [tstart]");
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
    } else if (cmd == ".END") {
        // Already handled in parseString
    } else {
        if (verbose_) {
            std::cout << "  Ignoring command: " << cmd << std::endl;
        }
    }
    
    return true;
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
    std::istringstream stream(line);
    std::string token;
    
    while (stream >> token) {
        tokens.push_back(token);
    }
    
    return tokens;
}

double NetlistParser::parseValue(const std::string& str) {
    // Supports multipliers: T, G, MEG, k, m, u, n, p, f
    std::string s = str;
    double multiplier = 1.0;
    
    // Check suffix
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
    
    double value = std::stod(s);
    return value * multiplier;
}

std::string NetlistParser::cleanLine(const std::string& line) {
    std::string cleaned;
    
    // Remove inline comments (after ;)
    size_t comment_pos = line.find(';');
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

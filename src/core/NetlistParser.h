/**
 * @author Gabriel Henrique Silva
 * @file NetlistParser.h
 * @brief SPICE netlist parser
 * 
 * Supports basic SPICE syntax:
 * - Comments: * or .comment
 * - Resistors: R<name> <n+> <n-> <value>
 * - Capacitors: C<name> <n+> <n-> <value> [IC=<initial_value>]
 * - Inductors: L<name> <n+> <n-> <value> [IC=<initial_value>]
 * - DC Sources: V<name> <n+> <n-> DC <value>
 * - AC Sources: V<name> <n+> <n-> AC <mag> [phase]
 * - Current Sources: I<name> <n+> <n-> DC <value>
 * - Commands: .OP, .TRAN, .AC, .END
 * 
 * Example netlist:
 * ````
 * * Voltage Divider Circuit
 * V1 1 0 DC 10
 * R1 1 2 1k
 * R2 2 0 2k
 * .OP
 * .END
 * ````
 */

#ifndef NETLIST_PARSER_H
#define NETLIST_PARSER_H

#include "../core/Circuit.h"
#include "../devices/Diode.h"
#include <string>
#include <vector>
#include <map>

namespace e_XSpice {

class SourceWaveform;  // devices/SourceWaveform.h

/**
 * @brief Analysis type requested in netlist
 */
enum class AnalysisType {
    NONE,
    DC_OP,        // .OP
    TRANSIENT,    // .TRAN
    AC,           // .AC
    DC_SWEEP      // .DC
};

/**
 * @brief Analysis parameters extracted from netlist
 */
struct AnalysisParameters {
    AnalysisType type = AnalysisType::NONE;
    
    // Transient analysis
    double tran_tstep = 0.0;    // Time step
    double tran_tstop = 0.0;    // Final time
    double tran_tstart = 0.0;   // Initial time (optional)
    
    // AC analysis
    std::string ac_variation = "DEC";  // DEC, OCT, LIN
    int ac_points = 10;
    double ac_fstart = 1.0;
    double ac_fstop = 1e6;
    
    // DC sweep
    std::string dc_source;
    double dc_start = 0.0;
    double dc_stop = 0.0;
    double dc_step = 0.0;
};

/**
 * @brief SPICE netlist parser
 */
class NetlistParser {
public:
    /**
     * @brief Constructor
     */
    NetlistParser();
    
    /**
     * @brief Destructor
     */
    ~NetlistParser() = default;
    
    /**
     * @brief Parses a netlist file
     * @param filename Path to .cir file
     * @param circuit Circuit where to add components (output)
     * @return true if success, false if error
     */
    bool parse(const std::string& filename, Circuit& circuit);
    
    /**
     * @brief Parses a netlist from string
     * @param netlist_str String containing the netlist
     * @param circuit Circuit where to add components (output)
     * @return true if success, false if error
     */
    bool parseString(const std::string& netlist_str, Circuit& circuit);
    
    /**
     * @brief Gets analysis parameters found
     */
    const AnalysisParameters& getAnalysisParameters() const {
        return analysis_params_;
    }
    
    /**
     * @brief Gets error message (if parse failed)
     */
    std::string getErrorMessage() const { return error_message_; }
    
    /**
     * @brief Set verbose mode
     */
    void setVerbose(bool verbose) { verbose_ = verbose; }

private:
    bool verbose_;
    std::string error_message_;
    int line_number_;
    AnalysisParameters analysis_params_;
    
    // Mapping of node names to IDs
    std::map<std::string, int> node_name_to_id_;
    int next_node_id_;
    
    // Model cards (.MODEL name D (...))
    std::map<std::string, DiodeModel> diode_models_;
    
    // Accumulates a card (with its '+' continuation lines) before parsing
    std::string continuation_buffer_;
    
    /**
     * @brief Processes a netlist line
     */
    bool parseLine(const std::string& line, Circuit& circuit);
    
    /**
     * @brief Parses resistor: R<name> <n+> <n-> <value>
     */
    bool parseResistor(const std::vector<std::string>& tokens, Circuit& circuit);
    
    /**
     * @brief Parses capacitor: C<name> <n+> <n-> <value> [IC=<ic>]
     */
    bool parseCapacitor(const std::vector<std::string>& tokens, Circuit& circuit);
    
    /**
     * @brief Parses inductor: L<name> <n+> <n-> <value> [IC=<ic>]
     */
    bool parseInductor(const std::vector<std::string>& tokens, Circuit& circuit);
    
    /**
     * @brief Parses voltage source: V<name> <n+> <n-> DC <value>
     */
    bool parseVoltageSource(const std::vector<std::string>& tokens, Circuit& circuit);
    
    /**
     * @brief Parses current source: I<name> <n+> <n-> DC <value>
     */
    bool parseCurrentSource(const std::vector<std::string>& tokens, Circuit& circuit);
    
    /**
     * @brief Parses diode: D<name> <anode> <cathode> <model> [AREA=<a>]
     */
    bool parseDiode(const std::vector<std::string>& tokens, Circuit& circuit);
    
    /**
     * @brief Parses commands (.OP, .TRAN, .AC, .MODEL, .END)
     */
    bool parseCommand(const std::vector<std::string>& tokens);
    
    /**
     * @brief Parses .MODEL <name> D (params...)
     */
    bool parseModel(const std::vector<std::string>& tokens);
    
    /**
     * @brief Parses diode model parameters (IS=, N=, RS=, CJO=...)
     */
    void parseDiodeParams(const std::vector<std::string>& params, DiodeModel& model);
    
    /**
     * @brief Parses a source waveform: DC/AC <v> | PULSE(...) | SIN(...) | PWL(...)
     */
    bool parseSourceSpec(const std::vector<std::string>& tokens,
                         size_t first_param,
                         SourceWaveform& waveform,
                         Circuit& circuit);
    
    /**
     * @brief Collects the numeric tokens of a waveform spec.
     */
    std::vector<double> collectNumbers(const std::vector<std::string>& tokens,
                                       size_t from);
    
    /**
     * @brief Gets or creates node ID from name
     * Names "0", "GND", "gnd" map to 0 (ground)
     */
    int getNodeId(const std::string& node_name);                       // name → id
    
    /**
     * @brief Returns the node object for an id, creating it when needed.
     */
    Node* getOrCreateNode(Circuit& circuit, const std::string& name, int id);
    
    /**
     * @brief Tokenizes a line (separates by spaces)
     */
    std::vector<std::string> tokenize(const std::string& line);
    
    /**
     * @brief Converts string to double with multipliers (k, m, u, n, p, f, meg)
     */
    double parseValue(const std::string& str);
    
    /**
     * @brief Removes comments and extra spaces from a line
     */
    std::string cleanLine(const std::string& line);
    
    /**
     * @brief Converts string to uppercase
     */
    std::string toUpper(const std::string& str);
    
    /**
     * @brief Sets error message
     */
    void setError(const std::string& msg);
};

} // namespace e_XSpice

#endif // NETLIST_PARSER_H

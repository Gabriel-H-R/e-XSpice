/**
 * @author Gabriel Henrique Silva
 * @file NetlistParser.h
 * @brief SPICE netlist parser (model-agnostic device dispatch)
 *
 * Supports SPICE syntax:
 * - Comments: * or .comment
 * - Resistors: R<name> <n+> <n-> [model] <value>
 * - Capacitors: C<name> <n+> <n-> <value> [IC=<initial_value>]
 * - Inductors: L<name> <n+> <n-> <value> [IC=<initial_value>]
 * - Diodes: D<name> <anode> <cathode> <model> [AREA=<a>] [M=<m>]
 * - Transistors: Q<name> <C> <B> <E> [S] <model> [AREA=<a>] [TEMP=<c>]
 *   (4 nodes = substrate terminal, needs a LEVEL=2 BJT card)
 * - MOSFETs: M<name> <D> <G> <S> <B> <model>
 *   [W=<w>] [L=<l>] [M=<m>] [AD=<a>] [AS=<a>] [TEMP=<c>]
 * - Sources: V/I<name> <n+> <n-> DC|AC|PULSE|SIN|PWL ...
 * - Models: .MODEL <name> <TYPE> (P1=V1 P2=V2 ...)  (any registered TYPE)
 * - Commands: .OP, .TRAN, .AC, .END
 *
 * The parser knows NO device physics: device letters dispatch through a
 * creator map and model cards through the ModelRegistry's factories.
 * Supporting a new device (Q, M, ...) = registering one creator +
 * one model factory. Zero parser-core changes.
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
#include "../models/core/ModelRegistry.h"
#include <functional>
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

    // Explicit requests, independent of order: a netlist may carry both
    // .OP and .TRAN (each card sets its own flag; `type` keeps the last
    // one for backward compatibility).
    bool want_op = false;     // .OP seen
    bool want_tran = false;   // .TRAN seen
    bool want_dc_sweep = false;  // .DC sweep seen
    bool want_ac = false;       // .AC seen
    
    // Transient analysis
    double tran_tstep = 0.0;    // Time step
    double tran_tstop = 0.0;    // Final time
    double tran_tstart = 0.0;   // Initial time (optional)
    
    // AC analysis
    std::string ac_variation = "DEC";  // DEC, OCT, LIN
    int ac_points = 10;
    double ac_fstart = 1.0;
    double ac_fstop = 1e6;
    
    // DC sweep (.DC [LIN|DEC|OCT] <source> <start> <stop> <step|points>)
    // LIN is the default when the keyword is omitted (.DC VIN -10 15 1).
    // For DEC/OCT, <step> carries the points-per-decade/octave count.
    std::string dc_variation = "LIN";
    std::string dc_source;
    double dc_start = 0.0;
    double dc_stop = 0.0;
    double dc_step = 0.0;

    // Optional SECOND (inner) sweep: `.DC <spec1> [<spec2>]`. Xyce runs
    // spec2 as the inner loop and restarts it from its start value at every
    // point of spec1 (the outer loop).
    bool dc_nested = false;
    std::string dc_variation2 = "LIN";
    std::string dc_source2;
    double dc_start2 = 0.0;
    double dc_stop2 = 0.0;
    double dc_step2 = 0.0;

    // ---- .TEMP ----
    // Circuit temperature, SPICE convention: degrees Celsius (a trailing
    // `C` is accepted and stripped; a trailing `K` marks Kelvin). Cards
    // that pin their own TEMP= keep it (instance > card > circuit).
    bool temp_given = false;
    double temp_value = 27.0;        // Celsius

    // ---- .OPTIONS ----
    // Recognized: NOFAIL (keep sweeping past a non-converged point),
    // RETRY=<n> (source-step refinement depth per failed point),
    // MAXITERS=<n> (Newton iteration cap), TEMP=<t> (same as `.TEMP`),
    // SWEEPPAR (distribute independent sweep points across ranks),
    // ACFREQPAR (distribute the AC frequency sweep across ranks — results
    // stay bit-identical, frequencies are independent).
    // Unknown tokens are collected in `option_warnings` instead of aborting
    // the netlist (forward compatible).
    bool opt_nofail = false;
    int  opt_max_retries = 3;
    int  opt_max_iterations = 0;     // 0 = engine default
    bool opt_sweep_parallel = false;
    bool opt_ac_freq_parallel = false;
    std::vector<std::string> option_warnings;
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

    /**
     * @brief Registers a device-letter creator ("Q", "M", ...).
     *        Creators receive the tokenized line and build devices through
     *        the circuit + this parser's helpers (nodes, values, models).
     */
    using DeviceCreator =
        std::function<bool(const std::vector<std::string>& tokens,
                           Circuit& circuit)>;
    void registerDeviceCreator(char letter, DeviceCreator creator);

    /**
     * @brief Read access to the model registry (cards defined so far).
     */
    const ModelRegistry& models() const { return models_; }
    ModelRegistry& models() { return models_; }

private:
    bool verbose_;
    std::string error_message_;
    int line_number_;
    AnalysisParameters analysis_params_;

    // Mapping of node names to IDs
    std::map<std::string, int> node_name_to_id_;
    int next_node_id_;

    // Shared model cards (.MODEL), filled by parseModel, read by devices.
    ModelRegistry models_;

    // Device-letter dispatch: "R" -> parseResistor, ... (extensible).
    std::map<char, DeviceCreator> creators_;
    
    // Accumulates a card (with its '+' continuation lines) before parsing
    std::string continuation_buffer_;
    
    /**
     * @brief Processes a netlist line
     */
    bool parseLine(const std::string& line, Circuit& circuit);
    
    /**
     * @brief Parses resistor: R<name> <n+> <n-> [model] <value>
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
     * @brief Parses diode: D<name> <anode> <cathode> <model> [AREA=<a>] [M=<m>]
     */
    bool parseDiode(const std::vector<std::string>& tokens, Circuit& circuit);

    /**
     * @brief Parses transistor: Q<name> <C> <B> <E> [S] <model>
     *        [AREA=<a>] [TEMP=<c>] (4 nodes need a LEVEL=2 card)
     */
    bool parseBjt(const std::vector<std::string>& tokens, Circuit& circuit);

    /**
     * @brief Parses MOSFET: M<name> <D> <G> <S> <B> <model>
     *        [W=<w>] [L=<l>] [M=<m>] [AD=<a>] [AS=<a>] [TEMP=<c>]
     */
    bool parseMosfet(const std::vector<std::string>& tokens, Circuit& circuit);
    
    /**
     * @brief Parses commands (.OP, .TRAN, .AC, .MODEL, .END)
     */
    bool parseCommand(const std::vector<std::string>& tokens);
    
    /**
     * @brief Parses .MODEL <name> <TYPE> (params...) via the registry.
     *        Works for every registered TYPE (D today, Q/M tomorrow).
     */
    bool parseModel(const std::vector<std::string>& tokens);
    
    /**
     * @brief Parses a source line's specifications (SPICE dual spec):
     *        `<time-spec> [AC <mag> [phase]]` where time-spec is
     *        DC <v> | PULSE(...) | SIN(...) | PWL(...) | omitted when only
     *        `AC` is given (which implies a DC value of 0).
     */
    bool parseSourceSpec(const std::vector<std::string>& tokens,
                         size_t first_param,
                         SourceWaveform& waveform,
                         Circuit& circuit);

    /**
     * @brief Parses the `AC <mag> [phase]` small-signal part (phase in
     *        degrees) and stores it on the waveform.
     */
    bool parseACSpec(const std::vector<std::string>& tokens,
                     size_t from,
                     SourceWaveform& waveform);
    
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
     * @brief Fills a per-instance option struct from trailing KEY=value
     *        tokens, generically from the struct's own parameter table
     *        (same machinery as model cards: aliases, °C keys, given-hooks,
     *        unknown-key errors with the valid list). New instance options
     *        never touch parser code — only the struct's table.
     * @return false + setError() on malformed/unknown options
     */
    template <class O>
    bool parseInstanceTail(const std::vector<std::string>& tokens,
                           size_t from, O& inst,
                           const std::vector<TypedParamDesc<O>>& table,
                           const char* what) {
        ParamMap opts;
        for (size_t i = from; i < tokens.size(); ++i) {
            const size_t eq = tokens[i].find('=');
            if (eq == std::string::npos) {
                setError(std::string(what) + ": expected KEY=value, got '" +
                         tokens[i] + "'");
                return false;
            }
            const std::string key = toUpper(tokens[i].substr(0, eq));
            try {
                opts[key] = parseValue(tokens[i].substr(eq + 1));
            } catch (const std::exception& e) {
                setError(std::string(what) + " option '" + key + "': " +
                         e.what());
                return false;
            }
        }
        const std::vector<std::string> unknown =
            applyParamTable(inst, table, opts);
        if (!unknown.empty()) {
            std::string list;
            for (const auto& u : unknown) {
                if (!list.empty()) list += " ";
                list += u;
            }
            setError(std::string(what) + ": unknown option(s) '" + list +
                     "' (valid options: " + validParamList(table) + ")");
            return false;
        }
        return true;
    }
    
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
     * @brief Parses a `.TEMP`/`.OPTIONS TEMP=` temperature.
     *
     * SPICE convention: degrees Celsius, with an optional trailing `C`;
     * a trailing `K` marks an explicit Kelvin value. Deliberately does NOT
     * go through parseValue(), whose `K` suffix means kilo.
     *
     * @param token   netlist token ("27", "27C", "300.15K", "-40")
     * @param celsius output, degrees Celsius
     * @return false when the token is not a temperature
     */
    static bool parseTemperature(const std::string& token, double& celsius);
    
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

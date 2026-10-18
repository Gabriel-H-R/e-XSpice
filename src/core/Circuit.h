/**
 * @author Gabriel Henrique Silva
 * @file Circuit.h
 * @brief Represents the complete circuit with nodes and devices
 *
 * Ownership model (RCP-based):
 * - Circuit OWNS all nodes and devices via RCP<T> (shared, refcounted).
 *   The lookup maps share ownership with the vectors, so erasing/lookup
 *   can never dangle and no manual delete loop exists.
 * - Everybody else (Device terminals, DaeSystem, analyses, parser) holds
 *   RAW OBSERVER pointers from .get(). Hot paths never copy an RCP: a copy
 *   costs an atomic RMW (~13ns measured) vs ~1ns for a raw copy, and under
 *   OpenMP the shared refcount line would bounce between cores.
 * - addDevice takes RCP<Device> BY VALUE so call sites move
 *   make_rcp<T>(...) temporaries in with zero atomic overhead.
 */

#ifndef CIRCUIT_H
#define CIRCUIT_H

#include "Node.h"
#include "Device.h"
#include <vector>
#include <map>
#include <memory>
#include <string>

namespace e_XSpice {

class Circuit {
public:
    /**
     * @brief Constructor
     * @param name Circuit name
     */
    Circuit(const std::string& name = "Circuit");

    /**
     * @brief Destructor (default: RCP members release automatically)
     */
    ~Circuit();

    /**
     * @brief Adds a node to the circuit
     * @param id Node ID
     * @param name Node name (optional)
     * @return Raw observer pointer to the created node (owned by Circuit)
     */
    Node* addNode(int id, const std::string& name = "");

    /**
     * @brief Gets node by ID (raw observer, owned by Circuit)
     */
    Node* getNode(int id);

    /**
     * @brief Gets node by name (raw observer, owned by Circuit)
     */
    Node* getNodeByName(const std::string& name);

    /**
     * @brief Adds device to the circuit (Circuit shares ownership)
     * @param device RCP handle; pass make_rcp<T>(...) temporaries (moved,
     *        zero atomic overhead)
     */
    void addDevice(RCP<Device> device);

    /**
     * @brief Gets device by name (raw observer, owned by Circuit)
     */
    Device* getDevice(const std::string& name);

    /**
     * @brief Returns all nodes (except ground)
     */
    const std::vector<RCP<Node>>& getNodes() const { return nodes_; }

    /**
     * @brief Returns all devices
     */
    const std::vector<RCP<Device>>& getDevices() const { return devices_; }

    /**
     * @brief Returns ground node (always ID 0; raw observer)
     */
    Node* getGroundNode() { return ground_node_.get(); }
    
    /**
     * @brief Total number of nodes (including ground)
     */
    int getNumNodes() const { return nodes_.size() + 1; }
    
    /**
     * @brief Number of devices
     */
    int getNumDevices() const { return devices_.size(); }
    
    /**
     * @brief Counts voltage sources (to size expanded matrix)
     */
    int countVoltageSources() const;
    
    /**
     * @brief Resets state of all nodes
     */
    void reset();
    
    /**
     * @brief Prints circuit topology
     */
    void printTopology() const;
    
    /**
     * @brief Prints solution (voltages and currents)
     */
    void printSolution() const;
    
    /**
     * @brief Validates circuit (floating nodes, etc)
     */
    bool validate(std::string& error_msg) const;
    
    std::string getName() const { return name_; }

    // ---------------------------------------------------------------
    // Circuit temperature (SPICE `.TEMP`)
    // ---------------------------------------------------------------

    /**
     * @brief Operating temperature in Kelvin (300.15 K = 27 °C by default).
     *
     * Applied by setTemperature() to every model card that did not pin its
     * own TEMP, and to every device that precomputes temperature-dependent
     * values. Idempotent: calling it with the same temperature twice does no
     * work, which is what makes a temperature sweep cheap.
     *
     * Cost per DISTINCT temperature is O(cards + devices): cards are
     * collected (deduplicated) once and re-scaled by shared pointer, while
     * devices refresh only the scalars they precomputed.
     *
     * Apply it after the netlist is built (main does), so every device sees
     * the temperature from its constructor onwards. Devices added later
     * still get their card re-scaled (the card list is invalidated), but
     * their construction-time scalars would start at TNOM.
     */
    void setTemperature(double kelvin);

    /**
     * @brief Operating temperature in Kelvin.
     */
    double temperature() const { return temperature_; }

    /**
     * @brief True once setTemperature() (or `.TEMP`) pinned the circuit
     *        temperature, i.e. cards were re-scaled to it.
     */
    bool temperatureGiven() const { return temperature_given_; }

private:
    std::string name_;
    RCP<Node> ground_node_;                        // Node 0 (reference)
    std::vector<RCP<Node>> nodes_;                 // All other nodes
    std::vector<RCP<Device>> devices_;             // All devices

    std::map<int, RCP<Node>> node_map_;             // ID -> Node
    std::map<std::string, RCP<Node>> node_name_map_; // Name -> Node
    std::map<std::string, RCP<Device>> device_map_;  // Name -> Device

    double temperature_ = 300.15;      // Kelvin
    bool   temperature_given_ = false;
    // Distinct cards reached through Device::temperatureCard(), rebuilt
    // lazily after any addDevice() so a temperature change is O(cards).
    // Handles are const (the const expresses "parameters are read-only"),
    // but temperature re-scaling is an explicit, controlled mutation of the
    // shared snapshot — exactly like Xyce's Model::temperatureStep.
    std::vector<const DeviceModel*> temperature_cards_;
    bool temperature_cards_dirty_ = true;
    bool temperature_cards_scanned_ = false;

    int next_node_id_;
};

} // namespace e_XSpice

#endif // CIRCUIT_H

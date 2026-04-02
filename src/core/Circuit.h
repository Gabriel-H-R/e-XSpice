/**
 * @author Gabriel Henrique Silva
 * @file Circuit.h
 * @brief Represents the complete circuit with nodes and devices
 * 
 * The Circuit class maintains:
 * - List of all nodes
 * - List of all devices
 * - Mapping of names to objects
 * - Interface for adding components
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
     * @brief Destructor
     */
    ~Circuit();
    
    /**
     * @brief Adds a node to the circuit
     * @param id Node ID
     * @param name Node name (optional)
     * @return Pointer to the created node
     */
    Node* addNode(int id, const std::string& name = "");
    
    /**
     * @brief Gets node by ID
     */
    Node* getNode(int id);
    
    /**
     * @brief Gets node by name
     */
    Node* getNodeByName(const std::string& name);
    
    /**
     * @brief Adds device to the circuit
     * @param device Pointer to the device (Circuit assumes ownership)
     */
    void addDevice(Device* device);
    
    /**
     * @brief Gets device by name
     */
    Device* getDevice(const std::string& name);
    
    /**
     * @brief Returns all nodes (except ground)
     */
    const std::vector<Node*>& getNodes() const { return nodes_; }
    
    /**
     * @brief Returns all devices
     */
    const std::vector<Device*>& getDevices() const { return devices_; }
    
    /**
     * @brief Returns ground node (always ID 0)
     */
    Node* getGroundNode() { return ground_node_; }
    
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

private:
    std::string name_;
    Node* ground_node_;                          // Node 0 (reference)
    std::vector<Node*> nodes_;                   // All other nodes
    std::vector<Device*> devices_;               // All devices
    
    std::map<int, Node*> node_map_;             // ID -> Node
    std::map<std::string, Node*> node_name_map_; // Name -> Node
    std::map<std::string, Device*> device_map_;  // Name -> Device
    
    int next_node_id_;
};

} // namespace e_XSpice

#endif // CIRCUIT_H

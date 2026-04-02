/**
 * @author Gabriel Henrique Silva
 * @file Node.h
 * @brief Represents a node (connection point) in the circuit
 * 
 * Each node has a unique ID and stores its voltage.
 * Node 0 is always the ground (reference).
 */

#ifndef NODE_H
#define NODE_H

#include <string>
#include <vector>

namespace e_XSpice {

class Device; // Forward declaration

class Node {
public:
    /**
     * @brief Default constructor
     * @param id Unique node identifier
     * @param name Node name (optional)
     */
    Node(int id, const std::string& name = "");
    
    /**
     * @brief Destructor
     */
    ~Node() = default;
    
    // Getters
    int getId() const { return id_; }
    std::string getName() const { return name_; }
    double getVoltage() const { return voltage_; }
    bool isGround() const { return id_ == 0; }
    
    // Setters
    void setVoltage(double v) { voltage_ = v; }
    void setName(const std::string& name) { name_ = name; }
    
    /**
     * @brief Adds a device connected to this node
     * @param device Pointer to the device
     */
    void addConnectedDevice(Device* device);
    
    /**
     * @brief Returns all connected devices
     */
    const std::vector<Device*>& getConnectedDevices() const {
        return connected_devices_;
    }
    
    /**
     * @brief Reset node for new calculation
     */
    void reset();

private:
    int id_;                              // Unique node ID
    std::string name_;                    // Node name (e.g., "VDD", "GND")
    double voltage_;                      // Voltage at the node (simulation result)
    std::vector<Device*> connected_devices_; // Connected devices
};

} // namespace e_XSpice

#endif // NODE_H

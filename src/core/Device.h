/**
 * @author Gabriel Henrique Silva
 * @file Device.h
 * @brief Abstract base class for all circuit devices
 * 
 * This class defines the interface that all devices must implement.
 * Each device contributes to the conductance matrix (G) and current vector (b)
 * using the Modified Nodal Analysis (MNA) technique.
 */

#ifndef DEVICE_H
#define DEVICE_H

#include <string>
#include <vector>
#include <memory>

namespace e_XSpice {

class Node;
class SparseMatrix;

/**
 * @brief Supported device types
 */
enum class DeviceType {
    RESISTOR,
    CAPACITOR,
    INDUCTOR,
    VOLTAGE_SOURCE,
    CURRENT_SOURCE,
    DIODE,
    MOSFET,
    BJT,
    UNKNOWN
};

/**
 * @brief Base class for all devices
 */
class Device {
public:
    /**
     * @brief Constructor
     * @param name Device name (e.g., "R1", "C2")
     * @param type Device type
     */
    Device(const std::string& name, DeviceType type);
    
    /**
     * @brief Virtual destructor
     */
    virtual ~Device() = default;
    
    // Getters
    std::string getName() const { return name_; }
    DeviceType getType() const { return type_; }
    const std::vector<Node*>& getNodes() const { return nodes_; }
    
    /**
     * @brief Returns string with device type
     */
    std::string getTypeString() const;
    
    /**
     * @brief Stamping in the MNA system
     * 
     * Each device contributes to the matrix G and vector b:
     * Gv = b
     * 
     * @param G Conductance matrix (sparse)
     * @param b Vector of known currents/voltages
     * @param node_map Mapping of nodes to matrix indices
     */
    virtual void stamp(SparseMatrix& G, std::vector<double>& b,
                      const std::vector<int>& node_map) = 0;
    
    /**
     * @brief Updates device state based on node voltages
     * 
     * Used for non-linear devices that need to update
     * their contributions based on the current solution.
     */
    virtual void updateState() {}
    
    /**
     * @brief Returns if the device is linear
     */
    virtual bool isLinear() const = 0;
    
    /**
     * @brief Calculates current through the device
     */
    virtual double getCurrent() const = 0;
    
    /**
     * @brief Calculates dissipated power
     */
    virtual double getPower() const = 0;

protected:
    std::string name_;           // Unique device name
    DeviceType type_;           // Device type
    std::vector<Node*> nodes_;  // Nodes connected to the device
};

} // namespace e_XSpice

#endif // DEVICE_H

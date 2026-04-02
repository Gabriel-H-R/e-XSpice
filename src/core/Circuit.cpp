/**
 * @author Gabriel Henrique Silva
 * @file Circuit.cpp
 * @brief Circuit class implementation
 */

#include "Circuit.h"
#include "../devices/VoltageSource.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace e_XSpice {

Circuit::Circuit(const std::string& name)
    : name_(name)
    , next_node_id_(1)  // 0 is reserved for ground
{
    // Creates ground node (always exists)
    ground_node_ = new Node(0, "GND");
    node_map_[0] = ground_node_;
    node_name_map_["GND"] = ground_node_;
    node_name_map_["0"] = ground_node_;
}

Circuit::~Circuit() {
    // Frees memory of all nodes
    delete ground_node_;
    for (auto node : nodes_) {
        delete node;
    }
    
    // Frees memory of all devices
    for (auto device : devices_) {
        delete device;
    }
}

Node* Circuit::addNode(int id, const std::string& name) {
    // Checks if ID already exists
    if (node_map_.find(id) != node_map_.end()) {
        throw std::runtime_error("Node ID " + std::to_string(id) + " already exists");
    }
    
    if (id == 0) {
        return ground_node_;  // Ground already exists
    }
    
    Node* node = new Node(id, name);
    nodes_.push_back(node);
    node_map_[id] = node;
    
    if (!name.empty()) {
        node_name_map_[name] = node;
    }
    
    // Updates next available ID
    if (id >= next_node_id_) {
        next_node_id_ = id + 1;
    }
    
    return node;
}

Node* Circuit::getNode(int id) {
    auto it = node_map_.find(id);
    if (it != node_map_.end()) {
        return it->second;
    }
    return nullptr;
}

Node* Circuit::getNodeByName(const std::string& name) {
    auto it = node_name_map_.find(name);
    if (it != node_name_map_.end()) {
        return it->second;
    }
    return nullptr;
}

void Circuit::addDevice(Device* device) {
    if (device == nullptr) {
        throw std::invalid_argument("Cannot add null device");
    }
    
    // Checks if name already exists
    if (device_map_.find(device->getName()) != device_map_.end()) {
        throw std::runtime_error("Device " + device->getName() + " already exists");
    }
    
    devices_.push_back(device);
    device_map_[device->getName()] = device;
}

Device* Circuit::getDevice(const std::string& name) {
    auto it = device_map_.find(name);
    if (it != device_map_.end()) {
        return it->second;
    }
    return nullptr;
}

int Circuit::countVoltageSources() const {
    int count = 0;
    for (const auto& device : devices_) {
        if (device->getType() == DeviceType::VOLTAGE_SOURCE) {
            count++;
        }
    }
    return count;
}

void Circuit::reset() {
    ground_node_->reset();
    for (auto node : nodes_) {
        node->reset();
    }
}

void Circuit::printTopology() const {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Circuit: " << name_ << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::cout << "\nNodes (" << getNumNodes() << "):" << std::endl;
    std::cout << "  0: GND (ground)" << std::endl;
    for (const auto& node : nodes_) {
        std::cout << "  " << node->getId() << ": " << node->getName() << std::endl;
    }
    
    std::cout << "\nDevices (" << getNumDevices() << "):" << std::endl;
    for (const auto& device : devices_) {
        const auto& nodes = device->getNodes();
        std::cout << "  " << device->getName() << " (" 
                 << device->getTypeString() << "): ";
        
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (i > 0) std::cout << " - ";
            std::cout << nodes[i]->getName();
        }
        std::cout << std::endl;
    }
    std::cout << "========================================\n" << std::endl;
}

void Circuit::printSolution() const {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Solution for: " << name_ << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::cout << "\nNode Voltages:" << std::endl;
    std::cout << std::setw(10) << "Node" 
             << std::setw(15) << "Voltage (V)" << std::endl;
    std::cout << std::string(25, '-') << std::endl;
    
    std::cout << std::setw(10) << "GND" 
             << std::setw(15) << std::fixed << std::setprecision(6) 
             << 0.0 << std::endl;
    
    for (const auto& node : nodes_) {
        std::cout << std::setw(10) << node->getName() 
                 << std::setw(15) << std::fixed << std::setprecision(6) 
                 << node->getVoltage() << std::endl;
    }
    
    std::cout << "\nDevice Currents and Power:" << std::endl;
    std::cout << std::setw(15) << "Device" 
             << std::setw(15) << "Current (A)" 
             << std::setw(15) << "Power (W)" << std::endl;
    std::cout << std::string(45, '-') << std::endl;
    
    for (const auto& device : devices_) {
        std::cout << std::setw(15) << device->getName()
                 << std::setw(15) << std::scientific << std::setprecision(6)
                 << device->getCurrent()
                 << std::setw(15) << std::scientific << std::setprecision(6)
                 << device->getPower() << std::endl;
    }
    
    std::cout << "========================================\n" << std::endl;
}

bool Circuit::validate(std::string& error_msg) const {
    // 1. Checks if there is at least one node besides ground
    if (nodes_.empty()) {
        error_msg = "Circuit has no nodes besides ground";
        return false;
    }
    
    // 2. Checks if there are devices
    if (devices_.empty()) {
        error_msg = "Circuit has no devices";
        return false;
    }
    
    // 3. Checks for floating nodes (no connected devices)
    for (const auto& node : nodes_) {
        if (node->getConnectedDevices().empty()) {
            error_msg = "Node " + node->getName() + " is floating (no devices connected)";
            return false;
        }
    }
    
    // 4. Checks if there is at least one path to ground
    // (simplified - only checks if any device connects to ground)
    bool has_ground_connection = false;
    for (const auto& device : devices_) {
        const auto& nodes = device->getNodes();
        for (const auto& node : nodes) {
            if (node->isGround()) {
                has_ground_connection = true;
                break;
            }
        }
        if (has_ground_connection) break;
    }
    
    if (!has_ground_connection) {
        error_msg = "No device connects to ground";
        return false;
    }
    
    return true;
}

} // namespace e_XSpice

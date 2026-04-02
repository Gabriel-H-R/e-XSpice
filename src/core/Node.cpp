/**
 * @author Gabriel Henrique Silva
 * @file Node.cpp
 * @brief Node class implementation
 */

#include "Node.h"
#include <algorithm>

namespace e_XSpice {

Node::Node(int id, const std::string& name)
    : id_(id)
    , name_(name.empty() ? std::to_string(id) : name)
    , voltage_(0.0)
{
    // Ground (node 0) always has voltage 0
    if (isGround()) {
        voltage_ = 0.0;
    }
}

void Node::addConnectedDevice(Device* device) {
    if (device != nullptr) {
        // Avoid duplicates
        auto it = std::find(connected_devices_.begin(), 
                           connected_devices_.end(), 
                           device);
        if (it == connected_devices_.end()) {
            connected_devices_.push_back(device);
        }
    }
}

void Node::reset() {
    if (!isGround()) {
        voltage_ = 0.0;
    }
}

} // namespace e_XSpice

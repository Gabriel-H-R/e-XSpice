/**
 * @author Gabriel Henrique Silva
 * @file Device.cpp
 * @brief Device base class implementation.
 */

#include "Device.h"
#include "../solver/DaeMatrix.h"
#include "DaeSolution.h"
#include "Node.h"

namespace e_XSpice {

Device::Device(const std::string& name, DeviceType type)
    : name_(name), type_(type) {}

void Device::emit(DaeMatrix& J, int entry, double value) const {
    assert(entry >= 0 && entry < static_cast<int>(jac_slots_.size()));
    const int slot = jac_slots_[entry];
    if (slot >= 0) J.accumulateSlot(slot, value);
}

std::string Device::getTypeString() const {    switch (type_) {
        case DeviceType::RESISTOR:        return "Resistor";
        case DeviceType::CAPACITOR:       return "Capacitor";
        case DeviceType::INDUCTOR:        return "Inductor";
        case DeviceType::VOLTAGE_SOURCE:  return "Voltage Source";
        case DeviceType::CURRENT_SOURCE:  return "Current Source";
        case DeviceType::DIODE:           return "Diode";
        case DeviceType::MOSFET:          return "MOSFET";
        case DeviceType::BJT:             return "BJT";
        default:                          return "Unknown";
    }
}

} // namespace e_XSpice

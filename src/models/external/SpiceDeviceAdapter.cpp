/**
 * @author Gabriel Henrique Silva
 * @file SpiceDeviceAdapter.cpp
 * @brief Base class implementation for wrapping external SPICE device models.
 *
 * See SpiceDeviceAdapter.h and src/models/external/README.md for the
 * architecture and the three-step recipe for adding a new model (EKV,
 * BSIM-CMG, BSIM-SOI, ...).
 */

#include "SpiceDeviceAdapter.h"

#include "core/Device.h"
#include "core/DaeSolution.h"
#include "solver/DaeMatrix.h"

#include <algorithm>
#include <stdexcept>

namespace e_XSpice {

SpiceDeviceAdapter::SpiceDeviceAdapter(const std::string& name, DeviceType type,
                                       const std::vector<Node*>& nodes)
    : Device(name, type) {
    nodes_ = nodes;
    for (Node* n : nodes_) {
        if (n) n->addConnectedDevice(this);
    }
}

void SpiceDeviceAdapter::registerJacobianPattern(std::vector<int>& rows,
                                                std::vector<int>& cols) const {
    // The pattern is discovered by the Berkeley device's setup function.
    // Subclasses implement discoverPattern() to run that registration and
    // declare the resulting (row, col) pairs in the e-XSpice pattern.
    std::vector<int> r, c;
    discoverPattern(r, c);
    for (size_t i = 0; i < r.size(); ++i) {
        if (r[i] < 0 || c[i] < 0) continue;
        rows.push_back(r[i]);
        cols.push_back(c[i]);
    }
}

void SpiceDeviceAdapter::setJacobianSlots(const std::vector<int>& slots) {
    slots_ = slots;
    Device::setJacobianSlots(slots);
}

void SpiceDeviceAdapter::loadF(const DaeSolution& x, double time,
                               std::vector<double>& F) const {
    (void)x; (void)time; (void)F;
    // The actual evaluation happens in evaluate(), which stamps F, Q, B
    // and the Jacobians through the SPICE shim. Subclasses call evaluate()
    // with the shimmed CKTcircuit.
    // The default is a no-op: subclasses that need the e-XSpice path
    // override this directly.
}

void SpiceDeviceAdapter::loadQ(const DaeSolution& x, double time,
                               std::vector<double>& Q) const {
    (void)x; (void)time; (void)Q;
    // See loadF().
}

void SpiceDeviceAdapter::loadB(double time, std::vector<double>& B) const {
    (void)time; (void)B;
    // See loadF().
}

void SpiceDeviceAdapter::loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                                       double time) const {
    (void)J; (void)x; (void)time;
    // See loadF().
}

void SpiceDeviceAdapter::loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                                       double time) const {
    (void)J; (void)x; (void)time;
    // See loadF().
}

} // namespace e_XSpice

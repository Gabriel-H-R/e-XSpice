/**
 * @author Gabriel Henrique Silva
 * @file StateIndex.h
 * @brief Per-device state handle (charge / flux history slots).
 *
 * Single definition shared by Device.h (device interface) and StateStore.h
 * (history storage): each device that stores charge/flux gets a contiguous
 * range of state slots.
 */

#ifndef STATE_INDEX_H
#define STATE_INDEX_H

namespace e_XSpice {

struct StateIndex {
    int base  = 0;   // global slot of this device's first state
    int count = 0;   // number of states owned by this device

    bool isValid() const { return count > 0; }
    int  sliceBase() const { return base; }
    int  sliceEnd()  const { return base + count; }
};

} // namespace e_XSpice

#endif // STATE_INDEX_H

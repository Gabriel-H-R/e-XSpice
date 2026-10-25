/**
 * @author Gabriel Henrique Silva
 * @file SpiceContext.h
 * @brief Minimal SPICE circuit context (CKTcircuit) and sparse matrix
 *        (SMPmatrix) shim for wrapping external device models written against
 *        the Berkeley SPICE device interface.
 *
 * WHY THIS EXISTS
 * ---------------
 * The Berkeley device models (BSIM4, BSIM-CMG, BSIM-SOI, EKV, ...) are
 * written against the SPICE device interface: their load functions take a
 * `CKTcircuit *` and a `GENinstance *`, read the solution from
 * `ckt->CKTstate0`, stamp the MNA matrix through `SMPaddElt()` and write
 * the RHS into `ckt->CKTrhs`.
 *
 * e-XSpice uses the DAE device contract (loadF / loadQ / loadJacobianF /
 * loadJacobianQ with pre-resolved CSR slots). This shim provides a thin
 * translation layer:
 *
 *   1. A minimal CKTcircuit whose CKTstate0/CKTrhs/CKTmatrix borrow
 *      e-XSpice's own vectors/matrix (no copying, no re-implementation).
 *   2. An SMPmatrix that maps SMPaddElt(row, col) to the pre-resolved slot
 *      of the e-XSpice DaeMatrix.
 *
 * The device's "pattern registration" (BSIM4setup) is replaced by manual
 * setup: the adapter populates BSIM4model/BSIM4instance directly and sets
 * the BSIM4*Ptr matrix-element pointers to point into a conductance cache.
 *
 * SCALABILITY
 * -----------
 * Adding a new model (EKV, BSIM-CMG, ...) is three steps:
 *   1. Vendor the Berkeley C sources into src/models/<model>/code/.
 *   2. Write a `<Model>Adapter : SpiceDeviceAdapter` subclass that:
 *        - populates the model/instance structs from the e-XSpice card,
 *        - connects BSIM4*Ptr pointers to the conductance cache,
 *        - calls the Berkeley load function through the shim.
 *   3. Register the adapter in the device factory.
 *
 * The physics is NEVER re-implemented: the Berkeley code is the single
 * source of truth, and the adapter only translates interfaces.
 */

#ifndef SPICE_CONTEXT_H
#define SPICE_CONTEXT_H

// Pull in the full SPICE type declarations (CKTcircuit, SMPmatrix,
// CKTnode, GENmodel, GENinstance, IFparm, IFvalue, mode flags, ...).
// These are defined in global namespace inside an extern "C" block in
// SpiceHeaders.h so that the Berkeley C code can use them unchanged.
#include "SpiceHeaders.h"

#include <cstddef>
#include <vector>

namespace e_XSpice {

// Bring the global SPICE types into the e_XSpice namespace so that
// e_XSpice::CKTcircuit IS ::CKTcircuit (same type, no cast needed when
// passing &ctx_ to Berkeley load functions).
using ::CKTcircuit;
using ::SMPmatrix;
using ::CKTnode;
using ::GENinstance;
using ::GENmodel;

/**
 * @brief SMPaddElt shim: returns a pointer to the (row, col) slot.
 *
 * The slot must have been registered during setup (the device's
 * registerJacobianPattern). If it was not, returns nullptr and the device
 * must handle that gracefully.
 */
inline double* SMPaddElt(SMPmatrix* matrix, int row, int col) {
    if (!matrix || row < 0 || col < 0 || row >= matrix->size ||
        col >= matrix->size) {
        return nullptr;
    }
    const int idx = row * matrix->size + col;
    if (idx < 0 || idx >= static_cast<int>(matrix->slot_of.size())) {
        return nullptr;
    }
    const int slot = matrix->slot_of[static_cast<size_t>(idx)];
    if (slot < 0) return nullptr;
    return &matrix->values[static_cast<size_t>(slot)];
}

/**
 * @brief SMPmakeElt: alias for SMPaddElt (some Berkeley code uses this).
 */
inline double* SMPmakeElt(SMPmatrix* matrix, int row, int col) {
    return SMPaddElt(matrix, row, col);
}

} // namespace e_XSpice

#endif // SPICE_CONTEXT_H

/**
 * @author Gabriel Henrique Silva
 * @file SpiceDeviceAdapter.h
 * @brief Base class for wrapping external device models written against the
 *        Berkeley SPICE device interface (BSIM4, EKV, BSIM-CMG, ...) into
 *        e-XSpice's DAE device contract.
 *
 * WHY THIS EXISTS
 * ---------------
 * e-XSpice's DAE contract (loadF / loadQ / loadJacobianF / loadJacobianQ /
 * loadB / loadNoiseSources) is the canonical interface: every device writes
 * into pre-resolved CSR slots, and the analysis engines are device-agnostic.
 *
 * The Berkeley device models (BSIM4, ...) are written against a DIFFERENT
 * interface: their load functions take a CKTcircuit* and a GENinstance*,
 * read the solution from ckt->CKTstate0, and stamp the matrix through
 * SMPaddElt().
 *
 * The SpiceDeviceAdapter bridges the two: it implements the DAE contract by
 * (a) building a minimal CKTcircuit that points at the e-XSpice solution
 *     and matrix, (b) calling the Berkeley load function, and (c) letting
 *     the SMPaddElt shim translate the stamps into the pre-resolved slots.
 *
 * ARCHITECTURE
 * ------------
 *   SpiceDeviceAdapter (abstract, implements Device)
 *       │
*       ├── Bsim4Adapter   (wraps BSIM4load / BSIM4setup / ...)
*       ├── EkVAdapter    (wraps EKV ...)
*       └── ...
 *
 * The adapter is responsible for:
 *   - declaring the Jacobian pattern (registerJacobianPattern): it runs the
 *     Berkeley device's pattern registration (or, for models that stamp
 *     dynamically, a representative evaluation) and declares the resulting
 *     (row, col) pairs in the e-XSpice pattern.
 *   - evaluating the device (loadF, loadQ, loadJacobianF, loadJacobianQ,
 *     loadB): it builds the CKTcircuit, calls the Berkeley load function,
*     and the stamps land in the pre-resolved slots through the shim.
 *
 * The physics is NEVER re-implemented: the Berkeley code is the single
 * source of truth, and the adapter only translates between the two
 * interfaces.
 *
 * SCALABILITY
 * -----------
 * Adding a new model (EKV, BSIM-CMG, ...) is three steps:
 *   1. Vendor the Berkeley C sources into src/models/<model>/code/.
 *   2. Write a <Model>Adapter : SpiceDeviceAdapter subclass that:
 *        - declares the model's parameter table (reusing the existing
 *          DeviceModel / ModelCard infrastructure),
 *        - implements the device-specific hooks (evaluate, noise, ...),
 *        - calls the Berkeley load function through the shim.
 *   3. Register the adapter in the device factory.
 *
 * The physics is NEVER re-implemented: the Berkeley code is the single
 * source of truth, and the adapter only translates between the two
 * interfaces.
 */

#ifndef SPICE_DEVICE_ADAPTER_H
#define SPICE_DEVICE_ADAPTER_H

#include "../../core/Device.h"
#include "SpiceContext.h"

#include <memory>
#include <string>
#include <vector>

namespace e_XSpice {

/**
 * @brief Base class for wrapping external SPICE device models.
 *
 * Subclasses implement the model-specific hooks (parameter table,
 * instance setup, evaluation). The base class provides the DAE contract
 * implementation that delegates to the Berkeley load function through the
 * SPICE shim.
 */
class SpiceDeviceAdapter : public Device {
public:
    SpiceDeviceAdapter(const std::string& name, DeviceType type,
                       const std::vector<Node*>& nodes);
    ~SpiceDeviceAdapter() override = default;

    // ---- DAE contract (implemented by the base, delegating to the shim) ----
    void registerJacobianPattern(std::vector<int>& rows,
                                std::vector<int>& cols) const override;
    void setJacobianSlots(const std::vector<int>& slots) override;
    void loadF(const DaeSolution& x, double time,
               std::vector<double>& F) const override;
    void loadQ(const DaeSolution& x, double time,
               std::vector<double>& Q) const override;
    void loadB(double time, std::vector<double>& B) const override;
    void loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;
    void loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;

    // ---- model-specific hooks (implemented by subclasses) ------------------
    /**
     * @brief Runs the Berkeley device's pattern registration (e.g.
     *        BSIM4setup) to discover the (row, col) pairs the device will
     *        stamp. The pairs are returned in `rows`/`cols` (with -1 for
     *        ground/absent entries, matching the DAE contract).
     */
    virtual void discoverPattern(std::vector<int>& rows,
                                 std::vector<int>& cols) const = 0;

    /**
     * @brief Evaluates the device at the given solution and time, stamping
     *        F, Q, B and the Jacobians through the SPICE shim.
     *
     * @param ctx   the shimmed SPICE circuit context
     * @param inst  the device instance (Berkeley GENinstance*)
     */
    virtual void evaluate(const CKTcircuit& ctx, void* inst) const = 0;

    // ---- the SPICE shim ---------------------------------------------------
    CKTcircuit& spiceContext() { return ctx_; }
    const CKTcircuit& spiceContext() const { return ctx_; }

    /**
     * @brief Sets the pointer to the DaeMatrix value array so that the
     *        SMPaddElt shim writes conductances directly into the matrix.
     *        Called by the analysis/engine before assembly.
     */
    void setMatrixValues(double* values) { smp_.values = values; }

    SMPmatrix& smpRef() { return smp_; }
    const SMPmatrix& smpRef() const { return smp_; }

protected:
    // The shimmed SPICE context (borrowed from the e-XSpice state).
    mutable CKTcircuit ctx_;
    // The borrowed solution/RHS/matrix pointers.
    mutable std::vector<double> sol_;
    mutable std::vector<double> rhs_;
    // The SMPmatrix shim (borrows the DaeMatrix's value array).
    mutable SMPmatrix smp_;
    // The resolved slots (from setJacobianSlots).
    std::vector<int> slots_;
};

} // namespace e_XSpice

#endif // SPICE_DEVICE_ADAPTER_H

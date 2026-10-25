/**
 * @author Gabriel Henrique Silva
 * @file BSim4Adapter.h
 * @brief Adapter wrapping the Berkeley BSIM4 4.8.3 C code (src/models/bsim4/code)
 *        into e-XSpice's DAE device contract, via the SPICE shim (SpiceContext.h).
 *
 * WHY THIS EXISTS
 * ---------------
 * The Berkeley BSIM4 C source files (b4ld.c, b4temp.c, ...) are written against
 * the SPICE3 device interface: their load functions take a `CKTcircuit *` and a
 * `GENinstance *`, read terminal voltages from `ckt->CKTstate0`, stamp
 * conductances through `SMPaddElt()`, and write current contributions to
 * `ckt->CKTrhs`.
 *
 * e-XSpice uses the DAE contract (loadF / loadJacobianF / ... with pre-resolved
 * CSR slots). This adapter:
 *
 *   1. Holds a BSIM4model + BSIM4instance allocated from the Berkeley C code,
 *      populated from the e-XSpice model card (MosfetBSIM4Model).
 *   2. Keeps a private state buffer (`bsim4_state_`) that the Berkeley code
 *      reads/writes through CKTstate0 — scratch storage for the terminal
 *      voltages (vbd, vbs, vgs, vds, ...) that BSIM4load needs.
 *   3. Builds a minimal CKTcircuit whose CKTstate0/CKTrhs/CKTmatrix borrow
 *      e-XSpice's own vectors/matrix (no copying, no re-implementation).
 *   4. On each evaluation, copies terminal voltages from x into bsim4_state_
 *      and calls BSIM4load() — the upstream Berkeley function, unchanged.
 *   5. The SMPaddElt shim translates the Berkeley matrix stamps into the
 *      pre-resolved DaeMatrix slots (O(1) per stamp).
 *
 * The physics is NEVER re-implemented: every equation in b4ld.c / b4temp.c
 * runs verbatim. The adapter only translates between the two interfaces.
 *
 * SCALABILITY
 * -----------
 * This is the template for every other Berkeley model (EKV, BSIM-CMG, ...):
 * vendor the C sources, write an adapter that populates the model/instance
 * structs and calls the Berkeley load function, register in the factory.
 */

#ifndef BSIM4_ADAPTER_H
#define BSIM4_ADAPTER_H

#include "SpiceDeviceAdapter.h"
#include "SpiceHeaders.h"
#include "bsim4def.h"

#include "../../core_models/MosfetBSIM4Model.h"

#include <memory>
#include <string>
#include <vector>

namespace e_XSpice {

/**
 * @brief BSIM4 4.8.3 adapter: wraps the Berkeley C load functions.
 *
 * The instance uses NSTATVARS state variables (matching the Berkeley state
 * indexing relative to state_base_):
 *
 *   state[0] = vbd,  state[1] = vbs,  state[2] = vgs,  state[3] = vds,
 *   state[4] = vdbs, ... up to state[28] = qs
 *
 * BSIM4load reads terminal voltages from CKTstate0[state_base_ + offset].
 */
class BSim4Adapter : public SpiceDeviceAdapter {
public:
    /**
     * @param model   Shared BSIM4 model card (from the registry).
     * @param nodes   Terminal nodes [drain, gate, source, bulk] (4 nodes).
     * @param W       Channel width (m).
     * @param L       Channel length (m).
     * @param M       Multiplicity (fingers).
     * @param area    Drain/source area multiplier.
     */
    BSim4Adapter(RCP<const MosfetBSIM4Model> model,
                 const std::vector<Node*>& nodes,
                 double W, double L, double M, double area);
    ~BSim4Adapter() override;

    // ---- Device contract --------------------------------------------------
    int numExtraUnknowns() const override { return 0; }
    std::string getTypeString() const { return "BSIM4"; }
    bool isLinear() const override { return false; }
    double getCurrent() const override { return id_; }
    double getPower() const override { return 0.0; }

    // ---- DAE contract -----------------------------------------------------
    void registerJacobianPattern(std::vector<int>& rows,
                                 std::vector<int>& cols) const override;

    void loadF(const DaeSolution& x, double time,
               std::vector<double>& F) const override;

    void loadQ(const DaeSolution& x, double time,
               std::vector<double>& Q) const override;

    void loadB(double time, std::vector<double>& B) const override;

    void loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;

    void loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;

    // ---- Temperature notification ----------------------------------------
    void setCircuitTemperature(double kelvin) override;

    // ---- Outputs ----------------------------------------------------------
    double getVth() const { return vth_; }
    double getGm() const { return gm_; }
    double getGds() const { return gds_; }
    double getGmbs() const { return gmbs_; }
    double getDrainCurrent() const { return id_; }
    double getSourceCurrent() const { return is_; }
    double getBulkCurrent() const { return ib_; }

    // ---- Setup hooks (called by the analysis/engine) ----------------------
    /**
     * @brief Sets the base index into CKTstate0 where this instance's state
     *        values live. Called once at setup.
     */
    void setStateBase(int base) { state_base_ = base; }

    /**
     * @brief Sets the pointer where BSIM4load writes its ceq (current)
     *        contributions. Called before loadF with F.data().
     */
    void setRhsPtr(double* rhs) { ctx_.CKTrhs = rhs; }

    /**
     * @brief Sets the pointer to the DaeMatrix value array so that the
     *        SMPaddElt shim writes conductances into the matrix during
     *        BSIM4load. Called before loadJacobianF.
     */
    void setMatrixValues(double* values) { smp_.values = values; }

    /**
     * @brief Sets the matrix row numbers of the four terminals. Called after
     *        topology is assigned. These are the values BSIM4load uses as
     *        indices into CKTrhs/CKTmatrix.
     */
    void setNodeRows(int rd, int rg, int rs, int rb);

    /**
     * @brief Connects all BSIM4*Ptr matrix-element pointers to the
     *        corresponding g_cache_ slots. Called after setJacobianSlots.
     */
    void setJacobianSlots(const std::vector<int>& slots) override;

    // ---- evaluate (called internally by loadF / loadJacobianF) -----------
    using SpiceDeviceAdapter::evaluate;  // avoid hiding base-class overload

    /**
     * @brief Runs BSIM4load() at the given solution, stamping ceq into the
     *        RHS buffer and conductances into the Jacobian cache.
     */
    void evaluate(const DaeSolution& x, double time) const;

    // ---- SpiceDeviceAdapter pure virtuals --------------------------------
    void discoverPattern(std::vector<int>& rows,
                         std::vector<int>& cols) const override;
    void evaluate(const CKTcircuit& ctx, void* inst) const override;

private:
    RCP<const MosfetBSIM4Model> model_;
    double W_, L_, M_, area_;
    int p_;  // polarity: +1 NMOS, -1 PMOS

    // ---- Berkeley structures (owned by this adapter) ---------------------
    BSIM4model*  bsim4_model_ = nullptr;  // the model card (C struct)
    BSIM4instance* bsim4_inst_ = nullptr; // the device instance (C struct)
    bsim4SizeDependParam* pParam_ = nullptr; // size-dependent (per L/W) params

    // ---- State base index (into CKTstate0 / bsim4_state_) ---------------
    int state_base_ = 0;  // CKTstate0[state_base_ + 0..NSTATVARS-1]

    // ---- Node matrix rows (set via setNodeRows after topology) -----------
    int dNode_ = -1, gNode_ = -1, sNode_ = -1, bNode_ = -1;

    // ---- Private state buffer (CKTstate0 points here) ---------------------
    // BSIM4load reads/writes terminal voltages through CKTstate0. NSTATVARS
    // entries per instance (indices 0..28 cover vbd..qs).
    mutable std::vector<double> bsim4_state_;

    // ---- Cached Jacobian values from last BSIM4load call ------------------
    mutable std::vector<double> g_cache_;   // conductance stamps (size = #slots)

    // ---- Cached results from last BSIM4load call ---------------------------
    mutable double id_ = 0.0, is_ = 0.0, ib_ = 0.0;
    mutable double vth_ = 0.0, gm_ = 0.0, gds_ = 0.0, gmbs_ = 0.0;
    mutable bool   cache_valid_ = false;
    mutable const DaeSolution* cached_x_ = nullptr;

    // ---- Setup flag --------------------------------------------------------
    bool setup_done_ = false;

    // ---- Internal helpers --------------------------------------------------
    void buildModelCard();   // populate bsim4_model_ from MosfetBSIM4Model
    void buildInstance();    // populate bsim4_inst_ with nodes/params
    void setupStates();      // assign state_base_ and instance state indices
    void setup();            // full BSIM4setup(): defaults + pParam + instance fields
    void copyVoltagesToState(const DaeSolution& x) const;
    void readResults() const;   // extract id_, gm_, etc. from the instance
};

} // namespace e_XSpice

#endif // BSIM4_ADAPTER_H

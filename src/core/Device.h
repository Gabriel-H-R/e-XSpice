/**
 * @author Gabriel Henrique Silva
 * @file Device.h
 * @brief Abstract base class + DAE device contract for all circuit devices.
 *
 * ================================================================
 * THE DAE CONTRACT (formulation used by Xyce, modernized here)
 * ================================================================
 *
 * The circuit is solved as a nonlinear differential-algebraic system:
 *
 *      C(x)·ẋ + F(x) = b(t)          (equivalently: F(x) + dQ/dt = 0)
 *
 * where x is the MNA unknown vector:
 *
 *      x = [ v_1 ... v_{n-1} | i_{b1} ... i_{bk} ]
 *            node voltages        branch currents of voltage-defined
 *                                devices (voltage sources, inductors)
 *
 * Each device declares its own:
 *
 *  F(x, t)   static residual:  Kirchhoff laws (KCL rows) + constitutive law
 *                              (branch rows of voltage-defined devices and
 *                              nonlinear devices like the diode)
 *  Q(x, t)   dynamic state:    charges (capacitors) and fluxes (inductors)
 *  J_F = dF/dx, J_Q = dQ/dx    the two Jacobian contributions
 *
 * The TIME INTEGRATOR lives in the analysis engine, NOT in the device
 * (this is the key Xyce design decision we adopt): the engine computes the
 * effective system
 *
 *      R = F + s·Q − s·Q_n − ṡ·Q_n_history
 *      J = J_F + s·J_Q                (s: integration coefficient)
 *
 * with s = 2/dt (trapezoidal) or s = 1/dt (backward Euler), then runs the
 * Newton-Raphson loop J·Δx = −R. Because devices never see integration
 * coefficients, one device implementation serves DC, transient and
 * future AC analyses.
 *
 * Performance design (improvement over both Xyce and the legacy code):
 *  - Devices declare their Jacobian SPARSITY PATTERN once
 *    (registerJacobianPattern) and receive pre-resolved NNZ slots
 *    (setJacobianSlots). Evaluation writes values directly into a
 *    fixed-pattern CSR array (DaeMatrix) with no maps, no trees, no
 *    re-sorting per Newton iteration — and assembly is parallelizable
 *    because every slot belongs to exactly one device.
 *  - Devices own their state history (StateStore slots), so history
 *    rotation and save/restore are O(#states) with no per-device
 *    allocations.
 */

#ifndef DEVICE_H
#define DEVICE_H

#include <string>
#include <vector>
#include <cassert>

#include "e_XSpice_rcp.h"   // RCP<T>: intrusive refcounted ownership
#include "Node.h"
#include "StateIndex.h"

namespace e_XSpice {

class Circuit;
class SparseMatrix;     // legacy signature compatibility (stampAC)
class DaeMatrix;
class DaeSolution;
class StateStore;

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
 * @brief Base class for all devices (DAE formulation).
 *
 * Inherits RefCounted so devices live in RCP<T> handles owned by Circuit.
 * Hot paths (assembly loops) keep using raw observer pointers obtained via
 * .get() — an RCP copy costs an atomic RMW (~13ns measured) while a raw
 * copy is ~1ns, so RCPs never cross the per-iteration boundary.
 */
class Device : public RefCounted {
public:
    Device(const std::string& name, DeviceType type);
    virtual ~Device() = default;

    // ---------------------------------------------------------------
    // Identity / topology
    // ---------------------------------------------------------------
    std::string getName() const { return name_; }
    DeviceType  getType() const { return type_; }
    const std::vector<Node*>& getNodes() const { return nodes_; }
    std::string getTypeString() const;

    /**
     * @brief Number of extra unknowns (branch currents) this device adds
     *        to the solution vector. R, C, I, D: 0. V, L: 1.
     */
    virtual int numExtraUnknowns() const = 0;

    /**
     * @brief Total width allocated to this device in the solution vector
     *        (2 external nodes + extra unknowns + internal nodes).
     *        Used by the analysis to sanity-check index assignment.
     */
    virtual int numSolutionVariables() const { return 2 + numExtraUnknowns(); }

    /**
     * @brief Jacobian stamp: (row, col) pairs this device can fill in the
     *        combined F/J_Q Jacobian. Called once after indices are known.
     *        The SAME pattern is used for dFdx and dQdx (Xyce constraint),
     *        so devices must register the union of both contributions.
     */
    virtual void registerJacobianPattern(std::vector<int>& rows,
                                         std::vector<int>& cols) const = 0;

    /**
     * @brief Receives the resolved NNZ slot positions, ALIGNED to the
     *        entries registered in registerJacobianPattern (same length and
     *        order). Slots may be -1 when the entry references the ground
     *        row/column (no matrix entry exists); emit() skips those.
     */
    virtual void setJacobianSlots(const std::vector<int>& slots) { jac_slots_ = slots; }

    /**
     * @brief Assigns the absolute matrix index of this device's extra
     *        unknown (branch current), or -1 if it has none.
     */
    virtual void setBranchIndex(int idx) { branch_index_ = idx; }
    int getBranchIndex() const { return branch_index_; }

    /**
     * @brief SETUP: informs the device about the matrix layout.
     *
     * Equivalent of Xyce's registerLIDs phase. Called by the analysis once,
     * after the node map is built and before the pattern is registered.
     * The base implementation stores the map so rowOf() works; devices may
     * override to derive additional local indices.
     */
    void setTopology(const std::vector<int>& node_map) {
        node_map_ = &node_map;
        assignIndices();
    }

    /**
     * @brief Registers the state (charge/flux) slots owned by this device.
     *        Returns how many slots were requested. Called once at setup;
     *        the handle remains valid for the simulation lifetime.
     */
    virtual int registerStates(StateStore& store) { (void)store; return 0; }

    /**
     * @brief Number of charge/flux states this device owns (0 for purely
     *        algebraic devices).
     */
    virtual int numStateVariables() const { return 0; }

    /**
     * @brief Gives the device its state handle (after registerStates).
     */
    void setStateIndex(const StateIndex& idx) { state_index_ = idx; }
    const StateIndex& getStateIndex() const { return state_index_; }

    // ---------------------------------------------------------------
    // EVALUATION (called on every Newton iteration / time step)
    // ---------------------------------------------------------------

    /**
     * @brief Loads F(x, t): static residual of this device.
     *        Sign convention: contribution to the DAE residual F(x)+dQ/dt=0,
     *        i.e., for a resistor: F[pos] += (v+ − v−)/R, F[neg] −= same.
     */
    virtual void loadF(const DaeSolution& x, double time,
                       std::vector<double>& F) const { (void)x; (void)time; (void)F; }

    /**
     * @brief Loads Q(x, t): charges/fluxes of this device.
     *        Capacitor: Q[pos] += C·(v+ − v−); Q[neg] −= same.
     *        Inductor: Q[branch] += L·i_L.
     */
    virtual void loadQ(const DaeSolution& x, double time,
                       std::vector<double>& Q) const { (void)x; (void)time; (void)Q; }

    /**
     * @brief Loads B(t): independent sources (Xyce keeps the source values in
     *        their own vector so the residual is F + dQ/dt − B = 0 and the
     *        waveform is evaluated once per time point, not per iteration).
     *          Voltage source: B[branch] += V(t)
     *          Current source: B[pos] −= I(t); B[neg] += I(t)
     */
    virtual void loadB(double time, std::vector<double>& B) const {
        (void)time; (void)B;
    }

    /**
     * @brief Loads the static Jacobian J_F = dF/dx.
     */
    virtual void loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                               double time) const { (void)J; (void)x; (void)time; }

    /**
     * @brief Loads the charge Jacobian J_Q = dQ/dx.
     */
    virtual void loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                               double time) const { (void)J; (void)x; (void)time; }

    /**
     * @brief Adds the charge-history terms to the RHS vector (transient).
     *
     * After the analysis has assembled  rhs = B − F − s·Q(x), this adds
     *   + s·Q_n(t)   and (trapezoidal)   + dQ/dt(t)
     * using the ACCEPTED history of the device's own states. The device
     * writes exactly the rows (with exactly the signs) that loadQ uses.
     *
     * @param store      accepted state history
     * @param q_coeff    s  (integration coefficient; 0 → DC, no-op)
     * @param qdot_coeff 1 for trapezoidal, 0 for backward Euler
     * @param rhs        right-hand side vector (output)
     */
    virtual void loadChargeHistory(const StateStore& store, double q_coeff,
                                   double qdot_coeff,
                                   std::vector<double>& rhs) const {
        (void)store; (void)q_coeff; (void)qdot_coeff; (void)rhs;
    }

    /**
     * @brief Stores Q(x) into the device's state history slot (next state).
     *        Called after every Newton iteration so the integrator can use
     *        the converged charge; the analysis then computes dQ/dt.
     */
    virtual void storeState(const DaeSolution& x, StateStore& store) const {
        (void)x; (void)store;
    }

    /**
     * @brief Called by the analysis after each ACCEPTED time step, so the
     *        device can cache its accepted dQ/dt (capacitor/diode current)
     *        for output purposes.
     */
    virtual void commitStates(const StateStore& store) { (void)store; }

    /**
     * @brief Reads the device's accepted dQ/dt from history (for output).
     *        Devices whose "current" is dQ/dt (capacitors) read it through
     *        the store by their own state index.
     */
    virtual double stateCurrent() const { return 0.0; }

    /**
     * @brief Recomputes the device's cached output quantities (current,
     *        power helpers) from a fully assembled solution vector.
     *
     * Needed by distributed simulation: each rank evaluates loadF/loadQ
     * only for the devices it owns, so output caches of unowned devices
     * would be stale. After the solve, the analysis calls this on EVERY
     * device with the identical global solution, making getCurrent() /
     * getPower() exact on all ranks. Must be side-effect free w.r.t. the
     * system (no matrix/vector writes) and deterministic.
     */
    virtual void refreshOutputs(const DaeSolution& x, double time) const {
        (void)x; (void)time;
    }

    // ---------------------------------------------------------------
    // OPERATING MODES
    // ---------------------------------------------------------------

    /**
     * @brief Sets the "initial-condition operating point" mode.
     *
     * When true (DC operating point run before a transient with IC=
     * present — the SPICE/Xyce convention), voltage-defined state devices
     * are enforced at their initial conditions:
     *   capacitor → behaves as a voltage source of value IC
     *   inductor  → behaves as a current source of value IC
     */
    virtual void setICOperatingPoint(bool on) { ic_operating_point_ = on; }
    bool isICOperatingPoint() const { return ic_operating_point_; }

    // ---------------------------------------------------------------
    // OUTPUTS / INTERROGATION
    // ---------------------------------------------------------------

    /**
     * @brief Terminal current (A), positive from node[0] to node[1].
     *        For reactive devices: dQ/dt of the stored state.
     */
    virtual double getCurrent() const = 0;

    /**
     * @brief Instantaneous power (W), positive = dissipating.
     */
    virtual double getPower() const = 0;

    /**
     * @brief True when the device contributes a linear (constant) stamp.
     */
    virtual bool isLinear() const = 0;

    /**
     * @brief Writes the LIMITED unknowns back into the solution vector.
     *
     * SPICE/Xyce limiting semantics: a nonlinear device that rejected the
     * iteration (acceptSolution() == false) has already computed the limited
     * junction voltages used to evaluate F/Q/J. Without a write-back, the
     * next Newton iteration starts again from the unlimited unknown — and
     * on a current-source-driven junction the raw unknown can be millions
     * of volts, so the walk-up (vt per iteration) never converges. Writing
     * the limited value back restarts the iteration AT the limited point.
     *
     * Called by the solver only after an iteration that failed convergence,
     * and only when `x` is replicated (serial path) — the distributed path
     * must keep x identical on every rank, which per-device flags cannot
     * guarantee. Default no-op (linear devices never limit).
     */
    virtual void applySolutionLimit(std::vector<double>& x) const { (void)x; }

    /**
     * @brief QUERY: is the current solution acceptable by this device?
     *        Nonlinear devices use it to reject iterations where limiting
     *        was active (Xyce: origFlag=false → isConverged()==false).
     *        Must be side-effect free.
     */
    virtual bool acceptSolution() const { return true; }

    /**
     * @brief Notifies the device that its solution was accepted, so it may
     *        clear per-iteration bookkeeping (e.g. limiting flags).
     */
    virtual void confirmAccept() {}

    /**
     * @brief Called at the start of each nonlinear solve. Devices that keep
     *        limiting history (e.g. the diode's PNJLIM) reset their history
     *        when initial_solve is true (first solve of the simulation).
     */
    virtual void beginSolve(bool initial_solve) { (void)initial_solve; }

    /**
     * @brief Seeds the initial guess x0 for Newton-Raphson.
     *        Devices may write their initial-condition values here
     *        (e.g., inductor i_L = IC) before the first solve.
     */
    virtual void seedInitialGuess(std::vector<double>& x) const { (void)x; }

    /**
     * @brief Analysis breakpoint requests (source switching times, etc).
     *        Called once at analysis start; devices may append times
     *        smaller than tstop.
     */
    virtual void requestBreakpoints(std::vector<double>& times,
                                    double tstop) const {
        (void)times; (void)tstop;
    }

    /**
     * @brief Upper bound for the next time step requested by this device
     *        (e.g. 0.1·PER for a PULSE source). 0 = no constraint.
     */
    virtual double maxTimeStepHint(double t) const { (void)t; return 0.0; }

protected:
    /**
     * @brief Fast Jacobian write: adds `value` to the device's NNZ slot for
     *        stamp entry #entry (order = registerJacobianPattern). Entries
     *        referencing the ground (slot < 0) are silently skipped.
     */
    void emit(DaeMatrix& J, int entry, double value) const;

    /**
     * @brief Matrix row of a node (−1 for ground), from the node map given
     *        by setTopology(). Available during evaluation and setup.
     */
    int rowOf(const Node* node) const {
        if (!node || node->isGround()) return -1;
        const int id = node->getId();
        if (!node_map_ || id < 0 || id >= static_cast<int>(node_map_->size())) {
            return -1;
        }
        return (*node_map_)[id];
    }

    /**
     * @brief Hook for devices that need extra local indices (rare; most
     *        devices only need rowOf() + branch_index_).
     */
    virtual void assignIndices() {}

    std::string name_;             // Unique device name
    DeviceType  type_;             // Device type
    std::vector<Node*> nodes_;     // Terminal nodes (0: +, 1: −)

    const std::vector<int>* node_map_ = nullptr;  // node id -> matrix row

    int        branch_index_ = -1; // Matrix index of extra unknown (or −1)
    StateIndex state_index_;       // Charge/flux history handle

    std::vector<int> jac_slots_;   // NNZ slots for the registered pattern

    bool ic_operating_point_ = false;
};

} // namespace e_XSpice

#endif // DEVICE_H

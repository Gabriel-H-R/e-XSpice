/**
 * @author Gabriel Henrique Silva
 * @file Capacitor.h
 * @brief Capacitor model — charge-based DAE formulation.
 *
 * Constitutive law:  I = dQ/dt,   Q = C·V   (constant capacitance)
 *
 * DAE formulation (Xyce-style):
 *
 *   NORMAL (DC transient, and DC when no IC=):
 *     Q[pos] += +C·(v+ − v−)        Q[neg] += −C·(v+ − v−)
 *     J_Q:   [ +C  −C ]
 *            [ −C  +C ]
 *     F  and J_F are ZERO (a pure charge device contributes nothing static).
 *     The time integrator turns J_Q, Q and the charge history into the
 *     companion model, so this device is integration-method agnostic.
 *
 *   INITIAL-CONDITION OPERATING POINT (DC run before .TRAN when IC= given):
 *     the capacitor behaves as a VOLTAGE SOURCE of value IC using a branch
 *     current unknown (the same trick Xyce uses: numIntVars = 1 when ICGiven):
 *       F[pos] += i_b        F[neg] -= i_b
 *       F[bra] += (v+ − v− − IC)
 *       J_F:  [pos][bra]=+1  [neg][bra]=-1  [bra][pos]=+1  [bra][neg]=-1
 *
 * The branch unknown stays allocated for the whole simulation (fixed MNA
 * layout); in transient it is deactivated with J_F[bra][bra] += 1 so the
 * row stays non-singular.
 */

#ifndef CAPACITOR_H
#define CAPACITOR_H

#include "../core/Device.h"
#include "../core/Node.h"
#include "../solver/StateStore.h"

namespace e_XSpice {

class Capacitor : public Device {
public:
    Capacitor(const std::string& name,
              Node* node_plus,
              Node* node_minus,
              double capacitance,
              double initial_voltage = 0.0,
              bool   ic_given = false);

    ~Capacitor() override = default;

    // ---------------------------------------------------------------
    // Setup
    // ---------------------------------------------------------------
    int numExtraUnknowns() const override { return ic_given_ ? 1 : 0; }
    int numStateVariables() const override { return 1; }
    int registerStates(StateStore& store) override;

    void registerJacobianPattern(std::vector<int>& rows,
                                 std::vector<int>& cols) const override;

    // ---------------------------------------------------------------
    // Evaluation
    // ---------------------------------------------------------------
    void loadF(const DaeSolution& x, double time,
               std::vector<double>& F) const override;

    void loadQ(const DaeSolution& x, double time,
               std::vector<double>& Q) const override;

    void loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;

    void loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;

    void storeState(const DaeSolution& x, StateStore& store) const override;

    // Pulls the accepted dQ/dt (capacitor current) out of the store after
    // every accepted time step.
    void commitStates(const StateStore& store) override;

    // Seeds the branch current unknown of the IC operating point.
    void seedInitialGuess(std::vector<double>& x) const override;

    // ---------------------------------------------------------------
    // Outputs
    // ---------------------------------------------------------------
    bool   isLinear() const override { return true; }
    double getCurrent() const override { return last_current_; }  // = dQ/dt
    double getPower()   const override;
    double getEnergy()  const { return 0.5 * capacitance_ * vcap() * vcap(); }

    // Getters
    double getCapacitance() const { return capacitance_; }
    double getInitialVoltage() const { return initial_voltage_; }
    double getVoltage() const { return vcap(); }
    bool   hasIC() const { return ic_given_; }

    void loadChargeHistory(const StateStore& store, double q_coeff,
                           double qdot_coeff,
                           std::vector<double>& rhs) const override;

private:
    double capacitance_;     // Capacitance (F)
    double initial_voltage_; // IC (V); valid when ic_given_
    bool   ic_given_;

    double last_current_ = 0.0;   // last accepted dQ/dt (A)
    mutable double vcap_cache_ = 0.0;  // v+ − v− of the last evaluation

    double vcap() const { return vcap_cache_; }

    /**
     * @brief True only when THIS capacitor enforces its initial condition
     *        as a voltage source (it owns a branch unknown only then).
     *        The global flag is set whenever any device in the circuit has
     *        IC=, so the per-device ic_given_ check is mandatory.
     */
    bool inICMode() const { return ic_given_ && isICOperatingPoint(); }

    void refreshVcap(const DaeSolution& x) const;
};

} // namespace e_XSpice

#endif // CAPACITOR_H

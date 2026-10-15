/**
 * @author Gabriel Henrique Silva
 * @file Inductor.h
 * @brief Inductor model — flux-based DAE formulation with branch current.
 *
 * Constitutive law:  V = L·dI/dt,   Φ = L·I (flux linkage)
 *
 * DAE formulation (Xyce-style): the inductor current IS an MNA unknown
 * (like a voltage source), so the device owns an extra row/column:
 *
 *   Q[branch] += L·i_L                          (flux linkage)
 *   dQdx[branch][branch] += L                    → dΦ/dt = L·di_L/dt
 *   F[pos] += i_L ; F[neg] -= i_L                (KCL rows)
 *   F[branch] += −(v+ − v−)                      (KVL row)
 *   dFdx: [pos][bra] += +1 ; [neg][bra] -= 1
 *         [bra][pos] -= 1 ; [bra][neg] += 1
 *
 * With the integrator applying −(s/h) to the Q terms, the branch equation
 * becomes  v+ − v− − (s/h)·L·i_L = 0, i.e. the trapezoidal/backward-Euler
 * companion model of the inductor — obtained WITHOUT any companion-model
 * code inside the device.
 *
 * INITIAL-CONDITION OPERATING POINT (DC with IC=): the inductor behaves as
 * a current source of value IC: F[pos] += IC; F[neg] −= IC; the branch row
 * is pinned (dFdx[bra][bra] = 1) and the branch current is seeded with IC,
 * exactly like Xyce does at the operating point.
 */

#ifndef INDUCTOR_H
#define INDUCTOR_H

#include "../core/Device.h"
#include "../core/Node.h"
#include "../solver/StateStore.h"

namespace e_XSpice {

class Inductor : public Device {
public:
    Inductor(const std::string& name,
             Node* node_plus,
             Node* node_minus,
             double inductance,
             double initial_current = 0.0,
             bool   ic_given = false);

    ~Inductor() override = default;

    // ---------------------------------------------------------------
    // Setup
    // ---------------------------------------------------------------
    int numExtraUnknowns() const override { return 1; }
    int numStateVariables() const override { return 1; }
    int registerStates(StateStore& store) override;

    void registerJacobianPattern(std::vector<int>& rows,
                                 std::vector<int>& cols) const override;

    void seedInitialGuess(std::vector<double>& x) const override;

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

    void commitStates(const StateStore& store) override;

    /**
     * @brief Recomputes cached outputs from the global solution (exact on
     *        every rank, including ranks that do not own this device).
     */
    void refreshOutputs(const DaeSolution& x, double time) const override;

    // ---------------------------------------------------------------
    // Outputs
    // ---------------------------------------------------------------
    bool   isLinear() const override { return true; }
    double getCurrent() const override { return last_current_; }
    double getPower()   const override;
    double getEnergy()  const { return 0.5 * inductance_ * last_current_ * last_current_; }

    // Getters
    double getInductance() const { return inductance_; }
    double getInitialCurrent() const { return initial_current_; }
    bool   hasIC() const { return ic_given_; }

    void loadChargeHistory(const StateStore& store, double q_coeff,
                           double qdot_coeff,
                           std::vector<double>& rhs) const override;

private:
    double inductance_;      // Inductance (H)
    double initial_current_; // IC (A)
    bool   ic_given_;

    mutable double last_current_ = 0.0;  // last accepted i_L (A)
    mutable double vcap_cache_   = 0.0;  // v+ − v− at last evaluation
    mutable double flux_cache_   = 0.0;  // L·i_L at last evaluation
};

} // namespace e_XSpice

#endif // INDUCTOR_H

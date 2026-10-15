/**
 * @author Gabriel Henrique Silva
 * @file Diode.h
 * @brief PN-junction diode — charge-based, robust formulation (Xyce level 1).
 *
 * The device is AGNOSTIC to the model contents: it holds a shared
 * `RCP<const DiodeModel>` (physics + precomputed temperature scaling, see
 * models/DiodeModel.h) plus its own DiodeInstanceParams (AREA, M, TEMP).
 * Any number of diodes can share one card with zero per-instance copies
 * and zero per-iteration temperature math — the thermal snapshot is read
 * through a pointer.
 *
 * ================================================================
 * CURRENT MODEL (three regions, Xyce N_DEV_Diode.C)
 * ================================================================
 * Let  Vt = kT/q,  Vte = N·Vt,  VteR = NR·Vt,  VteBRK = NBV·Vt.
 *
 *  (a) forward:  Vd ≥ −3·Vte
 *        I_norm = IS·(exp(Vd/Vte) − 1) + GMIN·Vd
 *        I_rec  = ISR·(exp(Vd/VteR) − 1)
 *        K_hi   = sqrt(IKF/(IKF + I_norm))          (high-injection knee)
 *        I_D    = K_hi·I_norm + K_gen·I_rec
 *  (b) reverse:   I_D = −IS·(1 + (3·Vte/(e·Vd))³) + GMIN·Vd
 *  (c) breakdown: I_D = −IS·exp(−(BV + Vd)/VteBRK) + GMIN·Vd
 * G_d = dI_D/dVd in every region; the series resistance RS is carried by
 * an internal node (Pri) with conductance Gspr = 1/RS. Everything is
 * scaled by AREA and multiplicity M.
 *
 * ================================================================
 * JUNCTION CHARGE (function of V, not an ODE state)
 * ================================================================
 * With C0 = CJO·AREA, u = 1 − Vd/VJ, sarg = u^(−M), tDep = FC·VJ:
 *   Vd < tDep:  Q = TT·I_D + VJ·C0·(1 − u·sarg)/(1 − M)
 *               C = TT·G_d + C0·sarg
 *   Vd ≥ tDep:  Q = TT·I_D + C0·tF1 + (C0/F2)·[F3·(Vd − tDep)
 *                           + ½·(M/VJ)·(Vd² − tDep²)]
 *               C = TT·G_d + (C0/F2)·[F3 + (M/VJ)·Vd]
 * with F2 = (1−FC)^(1+M), F3 = 1−FC·(1+M) — chosen so that Q and C are
 * BOTH continuous at Vd = FC·VJ (the SPICE3f5/Xyce construction; C⁰ in
 * charge and in capacitance).
 *
 * Because Q is an explicit function of V, dQ/dx = C·dV/dx is loaded into
 * dQdx and the integrator supplies the transient companion terms — the
 * device itself contains NO integration code.
 *
 * CONVERGENCE: SPICE3f5 PNJLIM on Vd (and on the mapped breakdown
 * variable), plus GMIN in parallel with every junction and an exponent
 * clamp at 100.
 */

#ifndef DIODE_H
#define DIODE_H

#include "../core/Device.h"
#include "../core/Node.h"
#include "../solver/StateStore.h"
#include "../models/core_models/DiodeModel.h"
#include <cmath>

namespace e_XSpice {

class Diode : public Device {
public:
    /**
     * @brief Full constructor: shared model + per-instance parameters.
     * @param model  Shared card (from ModelRegistry or make_rcp<DiodeModel>).
     *               Must already be processParams()'d.
     * @param inst   Instance overrides (AREA, multiplicity, TEMP).
     */
    explicit Diode(const std::string& name,
                   Node* anode,
                   Node* cathode,
                   RCP<const DiodeModel> model,
                   const DiodeInstanceParams& inst = DiodeInstanceParams());

    /**
     * @brief Convenience constructor with the process-wide default model.
     *        Same physics as before, zero setup.
     */
    explicit Diode(const std::string& name,
                   Node* anode,
                   Node* cathode);

    ~Diode() override = default;

    // ---------------------------------------------------------------
    // Setup
    // ---------------------------------------------------------------
    int numExtraUnknowns() const override { return (model_->RS > 0.0) ? 1 : 0; }
    int numStateVariables() const override { return 1; }
    int registerStates(StateStore& store) override;

    void registerJacobianPattern(std::vector<int>& rows,
                                 std::vector<int>& cols) const override;

    void beginSolve(bool initial_solve) override;

    // ---------------------------------------------------------------
    // Evaluation
    // ---------------------------------------------------------------
    void loadF(const DaeSolution& x, double time,
               std::vector<double>& F) const override;

    void loadQ(const DaeSolution& x, double time,
               std::vector<double>& Q) const override;

    void applySolutionLimit(std::vector<double>& x) const override;

    void loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;

    void loadJacobianQ(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;

    void storeState(const DaeSolution& x, StateStore& store) const override;

    void loadChargeHistory(const StateStore& store, double q_coeff,
                           double qdot_coeff,
                           std::vector<double>& rhs) const override;

    bool acceptSolution() const override;

    void confirmAccept() override;

    // ---------------------------------------------------------------
    // Outputs
    // ---------------------------------------------------------------
    bool   isLinear() const override { return false; }
    double getCurrent() const override;
    double getPower()   const override;
    double getVoltage() const { return junctionVoltage(); }

    /**
     * @brief Shared model card (observer; owned by registry/creator).
     */
    RCP<const DiodeModel> model() const { return model_; }

    /**
     * @brief Per-instance parameters.
     */
    const DiodeInstanceParams& instanceParams() const { return inst_; }

private:
    RCP<const DiodeModel> model_;   // shared card (never null)
    DiodeInstanceParams   inst_;    // per-instance overrides

    // Thermal snapshot in force: the model's, or the instance override's.
    // Plain pointer, setup-time binding, zero per-iteration cost.
    const DiodeThermal* thermal_ = nullptr;
    DiodeThermal ownThermal_;       // used only when inst_.TEMPGiven

    // junction operating state (last limited evaluation)
    mutable double vd_limited_ = 0.0;
    mutable double qd_limited_ = 0.0;
    mutable double cd_limited_ = 0.0;
    mutable double id_limited_ = 0.0;
    mutable double gd_limited_ = 0.0;
    mutable bool   limiting_active_ = false;

    mutable double vd_prev_ = 0.0;  // last limited junction voltage (PNJLIM history)
    mutable bool   first_solve_ = true;

public:
    // Debug switch (diagnostics for convergence studies)
    static void setDebug(bool on) { debug_ = on; }

private:
    static bool debug_;

    int priRow() const {
        return (model_->RS > 0.0) ? branch_index_ : rowOf(nodes_[0]);
    }

    /**
     * @brief Evaluates current/conductance and junction charge/capacitance
     *        at the (limited) voltage Vd; fills the cached members.
     *        Result includes AREA and multiplicity scaling.
     */
    void evaluateJunction(double vd) const;

    /**
     * @brief SPICE3f5 pnjlim — see N_DEV_DeviceSupport.C in Xyce.
     */
    static double pnjlim(double vnew, double vold, double vt, double vcrit);

    /**
     * @brief Reads the raw junction voltage from the solution and applies
     *        the limiters (returns the limited voltage).
     */
    double limitedJunctionVoltage(const DaeSolution& x) const;

    /**
     * @brief Conduction current at voltage vd (no cache mutation).
     */
    double currentForVoltage(double vd) const;

    double junctionVoltageFromNodes() const;
    double junctionVoltage() const { return junctionVoltageFromNodes(); }
};

} // namespace e_XSpice

#endif // DIODE_H

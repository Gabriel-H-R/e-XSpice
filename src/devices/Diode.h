/**
 * @author Gabriel Henrique Silva
 * @file Diode.h
 * @brief PN-junction diode — charge-based, robust formulation (Xyce level 1).
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
 *        I_D    = K_hi·I_norm + K_hi... (see .cpp for I_rec·K_gen)
 *  (b) reverse:   I_D = −IS·(1 + (3·Vte/(e·Vd))³) + GMIN·Vd
 *  (c) breakdown: I_D = −IS·exp(−(BV + Vd)/VteBRK) + GMIN·Vd
 * G_d = dI_D/dVd in every region; the series resistance RS is carried by
 * an internal node (Pri) with conductance Gspr = 1/RS.
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
#include <cmath>

namespace e_XSpice {

/**
 * @brief Diode model parameters (defaults per SPICE/Xyce level 1).
 */
struct DiodeModel {
    double IS  = 1e-14;   // Saturation current (A)
    double N   = 1.0;     // Emission coefficient
    double RS  = 0.0;     // Series resistance (Ω)
    double CJO = 0.0;     // Zero-bias junction capacitance (F)
    double VJ  = 1.0;     // Junction potential (V)
    double M   = 0.5;     // Grading coefficient
    double FC  = 0.5;     // Forward-bias depletion coefficient
    double BV  = 1e99;    // Reverse breakdown voltage (V); BVGiven=false → off
    double IBV = 1e-3;    // Current at breakdown (A)
    double NBV = 1.0;     // Breakdown emission coefficient
    double TT  = 0.0;     // Transit time (s)
    double EG  = 1.11;    // Bandgap (eV)
    double XTI = 3.0;     // IS temperature exponent
    double AREA = 1.0;    // Area multiplier
    double IKF = 0.0;     // High-injection knee current (A)
    double ISR = 0.0;     // Recombination saturation current (A)
    double NR  = 2.0;     // Recombination emission coefficient
    double TNOM = 300.15; // Reference temperature (K)
    double TEMP = 300.15; // Device temperature (K)
    double GMIN = 1e-12;  // Parallel conductance (S) per junction

    bool BVGiven = false;
};

class Diode : public Device {
public:
    using ModelParams = DiodeModel;

    explicit Diode(const std::string& name,
                   Node* anode,
                   Node* cathode,
                   const ModelParams& params = ModelParams());

    ~Diode() override = default;

    // ---------------------------------------------------------------
    // Setup
    // ---------------------------------------------------------------
    int numExtraUnknowns() const override { return (params_.RS > 0.0) ? 1 : 0; }
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

    const ModelParams& getParams() const { return params_; }

private:
    ModelParams params_;

    // temperature-scaled internals
    double vt_       = 0.0;   // kT/q
    double vte_      = 0.0;   // N·vt
    double vteR_     = 0.0;   // NR·vt
    double vteBRK_   = 0.0;   // NBV·vt
    double isat_     = 0.0;   // IS(T)
    double isatR_    = 0.0;   // ISR(T)
    double tJctPot_  = 0.0;   // VJ(T)
    double tJctCap_  = 0.0;   // CJO(T)
    double tDepCap_  = 0.0;   // FC·VJ(T)
    double F2_       = 0.0;   // (1−FC)^(1+M)
    double F3_       = 0.0;   // 1−FC·(1+M)
    double tF1_      = 0.0;   // depletion charge at the transition point
    double tBrkdwnV_ = 0.0;   // breakdown voltage fitted from BV/IBV
    double tVcrit_   = 0.0;   // critical voltage for PNJLIM
    double gspr_     = 0.0;   // 1/RS
    double ikf_      = 0.0;   // IKF(T)

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

    int priRow() const { return (params_.RS > 0.0) ? branch_index_ : rowOf(nodes_[0]); }

    void updateTemperature();
    void fitBreakdownVoltage();

    /**
     * @brief Evaluates current/conductance and junction charge/capacitance
     *        at the (limited) voltage Vd; fills the cached members.
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
    double expClamp(double arg) const { return std::exp(std::min(100.0, arg)); }
};

} // namespace e_XSpice

#endif // DIODE_H

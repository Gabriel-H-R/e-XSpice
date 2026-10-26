/**
 * @author Gabriel Henrique Silva
 * @file MosfetBSIM4Device.h
 * @brief BSIM4 4.8.3 device evaluation (level 2 MOSFET).
 *
 * Implements the BSIM4 DC current equations directly in C++, following the
 * Berkeley equations (b4ld.c / b4temp.c). The device is polarity-agnostic:
 * it receives N-oriented voltages (vgs, vds, vbs) and returns N-oriented
 * currents/charges/Jacobians. The Mosfet device wrapper handles the p=±1
 * polarity remapping and stamping.
 *
 * Equations implemented (BSIM4 4.8.3):
 *   - Threshold voltage with body effect (Vth = Vth0 + ...)
 *   - Effective mobility (surface & vertical field degradation)
 *   - Saturation velocity (velocity saturation)
 *   - Drain current: Ids = beta * Vgsteff^2 / ... (simplified)
 *   - Output conductance (LAMBDA)
 *   - Bulk diodes (SPICE3)
 *   - Charge-thickness model (capMod 1/2, Xyce C-V block) with overlap and
 *     junction charges; capacitances are central differences of the charges
 *
 * Not implemented (roadmap): pocket implant, NQS, gate current,
 * substrate current, AC-only parameters, noise.
 */

#ifndef MOSFET_BSIM4_DEVICE_H
#define MOSFET_BSIM4_DEVICE_H

#include "../../solver/DaeMatrix.h"
#include "../../core/DaeSolution.h"
#include "../../core/Node.h"
#include "../../models/core/DeviceConstants.h"
#include "MosfetBSIM4Model.h"
#include "Junction.h"
#include <cmath>

namespace e_XSpice {

/**
 * @brief BSIM4 4.8.3 device evaluation (level 2).
 *
 * DAE-facing shell: limiters, stamping and state live here. Physics lives
 * behind IBsim4Kernel (reference subset by default; Berkeley C via the
 * adapter skeleton when the shim is vendored). All voltages N-oriented.
 */
class MosfetBSIM4Device : public RefCounted {
public:
    /**
     * @param model  Shared BSIM4 model card (never null)
     * @param inst   Per-instance parameters (W, L, M, etc.)
     */
    MosfetBSIM4Device(RCP<const MosfetBSIM4Model> model,
                      double W, double L, double M, double AREA,
                      double SA = 0.0, double SB = 0.0, double SD = 0.0,
                      double SCA = 0.0, double SCB = 0.0, double SCC = 0.0,
                      double SC = 0.0);
    // Set by the Mosfet wrapper after construction.
    void setPolarity(int p) { p_ = p; }
    void setThermal(const MosfetBSIM4Thermal* t) { thermal_ = t; }

    // ---------------------------------------------------------------
    // Evaluation
    // ---------------------------------------------------------------
    /**
     * @brief Evaluates currents, charges and Jacobian ingredients at the
     *        given N-oriented controlling voltages.
     * @param vgs  gate-source voltage (N-oriented)
     * @param vds  drain-source voltage (N-oriented)
     * @param vbs  bulk-source voltage (N-oriented)
     */
    void evaluate(double vgs, double vds, double vbs);
    void applyScale(double lengthScale);
    void setIC(double icVDS, double icVGS, double icVBS);
    double eval1ovFNoise(double vds, double freq, double temp) const;
    int getNumNoiseSources() const;
    void updatePrimaryState();
    void auxChargeCalculations();
    void setupCapacitors();

    /**
     * @brief Reads raw terminal voltages from the solution, applies
     *        limiters, and stores the limited N-oriented voltages.
     * @param x     solution vector
     * @param rd,rg,rs,rb  matrix rows of drain,gate,source,bulk
     */
    void limitedVoltages(const DaeSolution& x,
                         int rd, int rg, int rs, int rb,
                         bool allow_limiting);

    /**
     * @brief Resets limiter history at the start of a solve (mirrors
     *        Mosfet::beginSolve/Bjt::beginSolve Xyce initJctFlag behaviour).
     */
    void beginSolve(bool initial_solve);

    /**
     * @brief Re-anchors limiter history to raw voltages (called by the
     *        wrapper's resetIterationState after a rejected step restore).
     */
    void resetHistory(double vgs_raw, double vds_raw, double vbs_raw);

    /**
     * @brief Gate-anchored write-back of the limited voltages into x
     *        (mirrors Mosfet::applySolutionLimit).
     */
    void applySolutionLimit(std::vector<double>& x,
                            int rd, int rg, int rs, int rb, int p) const;

    /**
     * @brief Advances limiter history after an evaluation (mirrors the
     *        vbs_prev_/vgs_prev_ updates of the level-1 loaders) and
     *        clears the first-solve flag so later iterations can limit.
     */
    void commitLimited();

    bool acceptSolution() const { return !limiting_active_; }
    void confirmAccept() { limiting_active_ = false; }

    double vgsLim() const { return vgs_lim_; }
    double vdsLim() const { return vds_lim_; }
    double vbsLim() const { return vbs_lim_; }
    bool limitingActive() const { return limiting_active_; }

    // ---------------------------------------------------------------
    // Loaders — called by Mosfet::loadF/loadQ/loadJacobianF/loadJacobianQ
    //
    // The residual/charge vectors are the GLOBAL nodal vectors (rows =
    // nodes + branch unknowns), so the terminal rows must be passed in:
    // this device knows its currents and charges, not the MNA layout.
    // The N-oriented currents/charges are multiplied by the polarity here,
    // exactly like the level-1 wrapper does.
    // ---------------------------------------------------------------
    void loadF(std::vector<double>& F,
               int rd, int rg, int rs, int rb) const;
    void loadQ(std::vector<double>& Q,
               int rd, int rg, int rs, int rb) const;
    void loadJacobianF(DaeMatrix& J,
                       const std::vector<int>& slots) const;
    void loadJacobianQ(DaeMatrix& J,
                       const std::vector<int>& slots) const;

    // ---------------------------------------------------------------
    // Outputs
    // ---------------------------------------------------------------
    double getDrainCurrent() const { return id_; }
    double getSourceCurrent() const { return is_; }
    double getBulkCurrent() const { return ib_; }
    double getGateCurrent() const { return 0.0; }
    double getVth() const { return vth_cached_; }
    double getGm() const { return gm_; }
    double getChargeG() const { return qg_; }
    double getChargeD() const { return qd_; }
    double getChargeB() const { return qb_; }
    double getChargeDef() const { return qdef_; }
    bool inSaturation() const { return saturated_; }

private:
    RCP<const MosfetBSIM4Model> model_;
    double W_ = 0.0, L_ = 0.0, M_ = 1.0, AREA_ = 1.0;
    double SA_ = 0.0, SB_ = 0.0, SD_ = 0.0;
    double SCA_ = 0.0, SCB_ = 0.0, SCC_ = 0.0, SC_ = 0.0;

    // Polarity multiplier: +1 NMOS, −1 PMOS (set by Mosfet wrapper).
    int p_ = 1;

    // Thermal snapshot from the model card.
    const MosfetBSIM4Thermal* thermal_ = nullptr;

    // ---- last evaluation (N-oriented laws, physical reports) ----
    mutable double vgs_cached_ = 0.0;  // EXTERNAL terminal voltages
    mutable double vds_cached_ = 0.0;
    mutable double vbs_cached_ = 0.0;
    mutable double vth_cached_ = 0.0;  // N-oriented threshold
    mutable bool   saturated_ = false;
    mutable bool   swapped_ = false;   // vds_n < 0 frame in force

    // N-oriented limited controlling voltages (for the write-back).
    mutable double vgs_lim_ = 0.0;
    mutable double vds_lim_ = 0.0;
    mutable double vbs_lim_ = 0.0;

    // Terminal currents, N-oriented (positive INTO the device).
    mutable double id_ = 0.0, is_ = 0.0, ib_ = 0.0, ig_ = 0.0;

    // Nodal charges, N-oriented physical frame.
    mutable double qg_ = 0.0, qd_ = 0.0, qb_ = 0.0, qdef_ = 0.0;  // qdef for NQS

    // Channel derivatives (effective frame): gm/dId/dvgs, gds, gmbs.
    mutable double gm_ = 0.0, gds_ = 0.0, gmbs_ = 0.0;
    // Bulk-diode conductances (physical frame), for Jacobian consistency.
    mutable double gibs_ = 0.0, gibd_ = 0.0;

    // Charge Jacobian dQ/d(vgs,vds,vbs) of the N-oriented node charges
    // (rows: drain, gate, bulk), by central differences in evaluate().
    mutable double cq_[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0}};

    mutable bool limiting_active_ = false;
    mutable bool first_solve_ = true;
    mutable bool computing_derivatives_ = false;
    // PNJLIM history (N-oriented bulk voltages) + step-clamp history.
    mutable double vbs_prev_ = 0.0;
    mutable double vbd_prev_ = 0.0;
    mutable double vgs_prev_ = 0.0;
    mutable double vds_prev_ = 0.0;

    static bool debug_;
    static constexpr double MaxStep_ = 8.0;  // bounded-step clamp (V)

    static double pnjlim(double vnew, double vold, double vt, double vcrit);
    static double clampStep(double vnew, double vold, double maxStep);

    // Operating point and terminal charges at one bias (no Jacobian). The
    // flag suppresses the diagnostic dump for the probing evaluations.
    void evaluateBase(double vgs, double vds, double vbs, bool probe);
};

} // namespace e_XSpice

#endif // MOSFET_BSIM4_DEVICE_H

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
 *   - Meyer charges (simplified Meyer partition for BSIM4)
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
#include "../external/Bsim4Kernel.h"
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
                      double W, double L, double M, double AREA);

    ~MosfetBSIM4Device() = default;

    // Set by the Mosfet wrapper after construction.
    void setPolarity(int p) { p_ = p; }
    void setThermal(const MosfetBSIM4Thermal* t) { thermal_ = t; }
    void setKernel(RCP<IBsim4Kernel> k) { kernel_ = std::move(k); }

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
    void evaluate(double vgs, double vds, double vbs) const;

    /**
     * @brief Reads raw terminal voltages from the solution, applies
     *        limiters, and stores the limited N-oriented voltages.
     * @param x     solution vector
     * @param rd,rg,rs,rb  matrix rows of drain,gate,source,bulk
     */
    void limitedVoltages(const DaeSolution& x,
                         int rd, int rg, int rs, int rb,
                         bool allow_limiting) const;

    // ---------------------------------------------------------------
    // Loaders — called by Mosfet::loadF/loadQ/loadJacobianF/loadJacobianQ
    // ---------------------------------------------------------------
    void loadF(std::vector<double>& F) const;
    void loadQ(std::vector<double>& Q) const;
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
    bool inSaturation() const { return saturated_; }

private:
    RCP<const MosfetBSIM4Model> model_;
    RCP<IBsim4Kernel> kernel_;  // physics source (reference default)
    double W_ = 0.0, L_ = 0.0, M_ = 1.0, AREA_ = 1.0;

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
    mutable double id_ = 0.0, is_ = 0.0, ib_ = 0.0;

    // Nodal charges, N-oriented physical frame.
    mutable double qg_ = 0.0, qd_ = 0.0, qb_ = 0.0;

    // Channel derivatives (effective frame): gm/dId/dvgs, gds, gmbs.
    mutable double gm_ = 0.0, gds_ = 0.0, gmbs_ = 0.0;

    // Capacitances (dQ/dV)
    mutable double cgs_ = 0.0, cgd_ = 0.0, cgb_ = 0.0;
    mutable double cdg_ = 0.0, cdd_ = 0.0, cdb_ = 0.0;
    mutable double csg_ = 0.0, csd_ = 0.0, css_ = 0.0, csb_ = 0.0;
    mutable double cbg_ = 0.0, cbd_ = 0.0, cbs_ = 0.0, cbb_ = 0.0;

    mutable bool limiting_active_ = false;
    mutable bool first_solve_ = true;
    // PNJLIM history (N-oriented bulk voltages) + step-clamp history.
    mutable double vbs_prev_ = 0.0;
    mutable double vbd_prev_ = 0.0;
    mutable double vgs_prev_ = 0.0;
    mutable double vds_prev_ = 0.0;

    static bool debug_;
    static constexpr double MaxStep_ = 8.0;  // bounded-step clamp (V)

    static double pnjlim(double vnew, double vold, double vt, double vcrit);
    static double clampStep(double vnew, double vold, double maxStep);
};

} // namespace e_XSpice

#endif // MOSFET_BSIM4_DEVICE_H

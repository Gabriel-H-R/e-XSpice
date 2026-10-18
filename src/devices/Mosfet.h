/**
 * @author Gabriel Henrique Silva
 * @file Mosfet.h
 * @brief MOSFET (NMOS and PMOS) — level 1 (Shichman-Hodges) in the shared
 *        DAE formulation.
 *
 * ================================================================
 * TOPOLOGY AND POLARITY
 * ================================================================
 * Four terminals (SPICE order), no internal nodes:
 *
 *      M<name>  D  G  S  B  <model> [W=<w>] [L=<l>] [M=<m>] [TEMP=<degC>]
 *
 *      nodes_[0] = drain   nodes_[1] = gate
 *      nodes_[2] = source  nodes_[3] = bulk
 *
 * Polarity-agnostic like the BJT: the DEVICE remaps terminal voltages to
 * N-oriented laws  vgs_n = p·(v_g − v_s), vds_n = p·(v_d − v_s),
 * vbs_n = p·(v_b − v_s)  (p = +1 NMOS / −1 PMOS) and stamps currents with
 * the same p, so every Jacobian entry is IDENTICAL for NMOS and PMOS.
 * VTO follows the SPICE sign convention (negative for PMOS cards) and is
 * mapped to the N-oriented frame (Vth_n = p·VTO); GAMMA/PHI stay positive
 * for both (reverse bulk bias is negative in the N-oriented frame either
 * way). When vds_n < 0 source and drain swap roles internally (the device
 * evaluates the swapped frame and maps currents/derivatives/charges back;
 * the loaders see one unified stamp).
 *
 * ================================================================
 * DAE STAMPS (rows/cols d, g, s, b — full 4×4 pattern, row-major)
 * ================================================================
 * N-oriented laws (vds ≥ 0 branch; Ig = 0 always):
 *   id = Id + Ibd,   is = −Id + Ibs,   ib = −(Ibs + Ibd)
 * with Id/gm/gds/gmbs from the card and Ibs/Ibd/gibs/gibd the bulk diodes.
 * J_F entries (node derivatives, e.g. row d: d/dvg = gm,
 * d/dvd = gds − gibd, d/dvs = −(gm+gds+gmbs), d/dvb = gmbs + gibd;
 * every COLUMN sums to zero ✓).
 *
 * Charges (Meyer + overlap + bulk depletion, charge-conserving by KCL):
 *   Qg = Cgs·vgs + Cgd·vgd + Cgb·vgb,  Qs = −Cgs·vgs − Qbs,
 *   Qd = −Cgd·vgd − Qbd,               Qb = −Cgb·vgb + Qbs + Qbd
 * Three states: Qg, Qd, Qb (Qs = −Qg−Qd−Qb, never stored).
 *
 * ================================================================
 * CONVERGENCE
 * ================================================================
 * SPICE3f5 PNJLIM on both bulk diodes, a bounded-step clamp on the
 * controlling voltages (the channel is quadratic, not exponential, so a
 * step bound suffices where diodes need pnjlim), gmin on every junction,
 * and the limiting write-back shared with diode/BJT.
 */

#ifndef MOSFET_H
#define MOSFET_H

#include "../core/Device.h"
#include "../core/Node.h"
#include "../solver/StateStore.h"
#include "../models/core_models/MosfetModel.h"
#include "../models/core_models/MosfetBSIM4Model.h"
#include "../models/core_models/MosfetBSIM4Device.h"
#include <cmath>

namespace e_XSpice {

class Mosfet : public Device {
public:
    /**
     * @brief Full constructor: shared card + per-instance parameters.
     * @param model shared card (must already be processParams()'d)
     */
    Mosfet(const std::string& name,
           Node* drain,
           Node* gate,
           Node* source,
           Node* bulk,
           RCP<const DeviceModel> model,
           const MosfetInstanceParams& inst = MosfetInstanceParams());

    /**
     * @brief Convenience constructors with process-wide default cards.
     */
    static RCP<Mosfet> makeNmos(const std::string& name, Node* drain,
                                Node* gate, Node* source, Node* bulk);
    static RCP<Mosfet> makePmos(const std::string& name, Node* drain,
                                Node* gate, Node* source, Node* bulk);

    ~Mosfet() override;

    // ---------------------------------------------------------------
    // Setup
    // ---------------------------------------------------------------
    int numExtraUnknowns() const override { return 0; }
    int numSolutionVariables() const override { return 4; }
    int numStateVariables() const override { return 3; }
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

    void resetIterationState(const DaeSolution& x) override;

    void refreshOutputs(const DaeSolution& x, double time) const override;

    // ---------------------------------------------------------------
    // Outputs
    // ---------------------------------------------------------------
    bool   isLinear() const override { return false; }
    // Drain current, positive INTO the drain (terminal convention).
    double getCurrent() const override { return p_ * id_; }
    double getPower()   const override;

    double getVgs() const { return vgs_cached_; }
    double getVds() const { return vds_cached_; }
    double getVbs() const { return vbs_cached_; }
    double getVth() const { return vth_cached_; }
    double getDrainCurrent() const { return p_ * id_; }
    double getSourceCurrent() const { return p_ * is_; }
    double getBulkCurrent()  const { return p_ * ib_; }
    double getTransconductance() const { return gm_; }
    bool   inSaturation() const { return saturated_; }

    RCP<const DeviceModel> model() const { return model_; }
    const MosfetInstanceParams& instanceParams() const { return inst_; }

    /** @brief True when using a BSIM4 level-2 card. */
    bool isBSIM4() const { return (bool)bsim4Model_; }

    /** @brief Access the level-1 card when available (null when BSIM4). */
    RCP<const MosfetModel> level1Model() const {
        return rcp_dynamic_cast<const MosfetModel>(model_);
    }

    static void setDebug(bool on) { debug_ = on; }

    // ---- circuit temperature (`.TEMP`) ----
    RCP<const DeviceModel> temperatureCard() const override { return model_; }

    /**
     * @brief Re-derives the setup-time scalars that depend on the card's
     *        thermal snapshot (cox). The card itself is re-scaled by
     *        Circuit::setTemperature(); only the values this device
     *        precomputed for the hot path are refreshed here.
     */
    void setCircuitTemperature(double kelvin) override;

private:
    RCP<const DeviceModel> model_;
    MosfetInstanceParams   inst_;

    // BSIM4 level-2 support: shared card + device evaluator (null when
    // LEVEL=1 / Shichman-Hodges).
    RCP<const MosfetBSIM4Model> bsim4Model_;
    MosfetBSIM4Device* bsim4Device_ = nullptr;  // owned by this Mosfet

    const MosfetThermal* thermal_ = nullptr;
    MosfetThermal ownThermal_;
    MosfetBSIM4Thermal ownThermalBSIM4_;

    // Polarity multiplier: +1 NMOS, −1 PMOS.
    int p_ = 1;

    // Setup-time geometry (no per-iteration math).
    double leff_ = 0.0;   // L − 2·LD
    double beta_ = 0.0;   // KP·(W/Leff)·M
    double cox_ = 0.0;    // ε_ox/TOX·W·Leff·M (intrinsic Meyer total)
    double cgso_ = 0.0;   // overlap totals (×M already)
    double cgdo_ = 0.0;
    double cgbo_ = 0.0;

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

    // Unified physical-frame stamps, terminal order D(0) G(1) S(2) B(3),
    // row-major: ef_[r*4+c], eq_[r*4+c]. Filled by evaluate(); the loaders
    // only emit them — no branches in the hot path.
    mutable double ef_[16] = {};
    mutable double eq_[16] = {};

    mutable bool limiting_active_ = false;
    mutable bool first_solve_ = true;
    // PNJLIM history (N-oriented bulk voltages) + step-clamp history.
    mutable double vbs_prev_ = 0.0;
    mutable double vbd_prev_ = 0.0;
    mutable double vgs_prev_ = 0.0;
    mutable double vds_prev_ = 0.0;

    static bool debug_;
    static constexpr double MaxStep_ = 8.0;  // bounded-step clamp (V)

    void limitedVoltages(const DaeSolution& x, bool allow_limiting) const;

    /**
     * @brief Copies the BSIM4 device evaluation into this wrapper's
     *        members so storeState/getCurrent/getPower/refreshOutputs and
     *        the limiting vote all work unchanged (single source of truth
     *        stays the device; this is just a pull, no second evaluation).
     */
    void pullBSIM4() const;

    /**
     * @brief Evaluates channel + bulk + charge at the (already limited)
     *        N-oriented voltages; fills currents, charges and ef_/eq_.
     */
    void evaluate(double vgs, double vds, double vbs) const;

    static double pnjlim(double vnew, double vold, double vt, double vcrit);
    static double clampStep(double vnew, double vold, double maxStep);
};

} // namespace e_XSpice

#endif // MOSFET_H

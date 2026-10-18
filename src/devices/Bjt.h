/**
 * @author Gabriel Henrique Silva
 * @file Bjt.h
 * @brief Bipolar junction transistor (NPN and PNP) — level 1 (Ebers-Moll
 *        transport) and level 2 (Gummel-Poon) in the shared DAE formulation.
 *
 * ================================================================
 * TOPOLOGY AND POLARITY
 * ================================================================
 * Three external terminals (SPICE order), optional substrate:
 *
 *      Q<name>  C  B  E  [S]  <model name> [AREA=<a>] [TEMP=<degC>]
 *
 *      nodes_[0] = collector   nodes_[1] = base   nodes_[2] = emitter
 *      nodes_[3] = substrate (only when the netlist gives 4 nodes; needs a
 *                  LEVEL=2 card — LEVEL=1 with 4 nodes is rejected)
 *
 * The implementation is polarity-agnostic: the DEVICE remaps the terminal
 * voltages so that the CARD always evaluates the NPN-oriented laws
 *
 *      vje_B = p·(v_b − v_e),     vje_C = p·(v_b − v_c),    p = +1 NPN / −1 PNP
 *
 * and the terminal currents are stamped with the same p. Because the sign of
 * the junction voltage AND the sign of the current flip together, every
 * Jacobian entry is IDENTICAL for NPN and PNP (verified in the derivation
 * below) — there is exactly one code path for both polarities.
 *
 * ================================================================
 * LEVEL 2 INTERNAL NODES (ohmic RB/RC/RE)
 * ================================================================
 * Each nonzero ohmic resistance gets one internal unknown (a voltage, held
 * in the device's contiguous extra-unknown range, exactly like the diode's
 * Pri node): bi between external base and the intrinsic base, ci, ei ditto.
 * With RB=RC=RE=0 the intrinsic rows alias the external rows and level 2
 * costs exactly the qb evaluation — no extra unknowns, no extra stamp.
 * Fixed order of the internal range: base, collector, emitter (whichever
 * exist). The intrinsic junctions are evaluated on the INTRINSIC nodes;
 * the resistors contribute the usual 4-entry conductance quads.
 *
 * ================================================================
 * DAE STAMPS, level 1 (rows b, c, e; unknowns v_b, v_c, v_e)
 * ================================================================
 * Laws (courtesy of BjtModel, with β_F = IS/BF, β_R = IS/BR):
 *      i_b = I_BE(vbe) + I_BC(vbc)
 *      i_c = I_CT(vbe,vbc) − I_BC(vbc)
 *      i_e = −i_b − i_c
 *
 * F:  F[b] += p·i_b,  F[c] += p·i_c,  F[e] += p·i_e
 * Q:  Q[b] += p·(Q_BE + Q_BC),  Q[e] −= p·Q_BE,  Q[c] −= p·Q_BC
 *
 * J_F = dF/dx (writing gib=∂I_BE/∂vbe, gibc=∂I_BC/∂vbc, gm=∂I_CT/∂vbe,
 *               go=∂I_CT/∂vbc):
 *
 *        |  b            c             e      |
 *   -----+------------------------------------+ col b → v_b etc.
 *    F[b]| gib+gibc     −gibc         −gib     |
 *    F[c]| gm+go−gibc   −go+gibc      −gm      |
 *    F[e]| −(gib+gm+go)  go            gib+gm  |
 *
 *   (each COLUMN sums to zero — KCL consistency of the stamp ✓)
 *
 * J_Q = dQ/dx (writing Cbe_v=∂Q_BE/∂vbe, Cbe_w=∂Q_BE/∂vbc,
 *               Cbc_u=∂Q_BC/∂vbe, Cbc_v=∂Q_BC/∂vbc):
 *
 *        |  b                     c              e          |
 *    F[b]| Cbe_v+Cbe_w+Cbc_u+Cbc_v −Cbe_w−Cbc_v   −Cbe_v−Cbc_u|
 *    F[c]| −Cbc_u−Cbc_v            Cbc_v          Cbc_u      |
 *    F[e]| −Cbe_v−Cbe_w            Cbe_w          Cbe_v      |
 *
 *   (each COLUMN sums to zero — the device is charge neutral ✓)
 *
 * The cross terms (Cbe_w, Cbc_u) are the transit-charge derivatives: they
 * exist because the base charge is transported by I_CT, which depends on
 * both junction voltages. Dropping them would break Newton's quadratic
 * convergence in any switching circuit.
 *
 * Level 2 uses the SAME stamp shapes on the intrinsic nodes (plus the
 * resistor quads and the substrate quad); only the laws change (qb,
 * non-ideal base currents). See BjtModel.h for the Gummel-Poon equations.
 * Frame rule (applies to every stamp in this device): F/Q rows accumulate
 * N-oriented currents for the INTRINSIC laws (p² = 1 makes the block
 * polarity-free) but PHYSICAL currents for the symmetric elements
 * (RB/RC/RE ohmic drops — no p anywhere). Mixing the frames negates the
 * PNP equation and converges to a non-solution.
 *
 * ================================================================
 * CONVERGENCE
 * ================================================================
 * SPICE3f5 PNJLIM on both junction voltages (plus the Early/q1 clamp and
 * the exponent clamp inside the card evaluation), gmin on every junction,
 * and the limiting write-back (applySolutionLimit) shared with the diode.
 */

#ifndef BJT_H
#define BJT_H

#include "../core/Device.h"
#include "../core/Node.h"
#include "../solver/StateStore.h"
#include "../models/core_models/BjtModel.h"
#include <cmath>

namespace e_XSpice {

class Bjt : public Device {
public:
    /**
     * @brief Full constructor: shared card + per-instance parameters.
     * @param model shared card (must already be processParams()'d)
     * @param substrate optional 4th terminal (needs a LEVEL=2 card)
     */
    Bjt(const std::string& name,
        Node* collector,
        Node* base,
        Node* emitter,
        RCP<const BjtModel> model,
        const BjtInstanceParams& inst = BjtInstanceParams(),
        Node* substrate = nullptr);

    /**
     * @brief Convenience constructor with the process-wide default NPN card.
     */
    Bjt(const std::string& name,
        Node* collector,
        Node* base,
        Node* emitter);

    ~Bjt() override = default;

    // ---------------------------------------------------------------
    // Setup
    // ---------------------------------------------------------------
    int numExtraUnknowns() const override { return numInternal_; }
    int numSolutionVariables() const override {
        return 3 + numInternal_ + (hasSubstrate() ? 1 : 0);
    }
    int numStateVariables() const override {
        return 2 + (hasSubstrate() ? 1 : 0);
    }
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
    // Terminal convention: positive INTO the terminal. ib_/ic_/ie_ hold the
    // NPN-oriented card laws, so a PNP reports the sign-flipped values.
    double getCurrent() const override { return p_ * ic_; }
    double getPower()   const override;                    // dissipated power

    /** Base-emitter voltage of the last evaluation (terminal difference). */
    double getVbe() const { return vbe_cached_; }
    /** Base-collector voltage of the last evaluation (terminal difference). */
    double getVbc() const { return vbc_cached_; }
    double getBaseCurrent()      const { return p_ * ib_; }
    double getCollectorCurrent() const { return p_ * ic_; }
    double getEmitterCurrent()   const { return p_ * ie_; }
    /** Stored base-emitter charge of the last evaluation. */
    double getBaseCharge()       const { return p_ * qbe_; }
    /** Stored base-collector charge of the last evaluation. */
    double getCollectorCharge()  const { return p_ * qbc_; }
    /** Substrate current (into the substrate terminal, 0 without one). */
    double getSubstrateCurrent() const {
        // F[sub] −= p·isc_: the node→device current is −p·isc_.
        return hasSubstrate() ? -p_ * isc_ : 0.0;
    }
    /** Intrinsic (post-RB/RC/RE) base-emitter voltage, NPN-oriented. */
    double getIntrinsicVbe() const { return vbe_lim_; }
    /** Intrinsic base-collector voltage, NPN-oriented. */
    double getIntrinsicVbc() const { return vbc_lim_; }

    /**
     * @brief Shared model card (observer; owned by registry/creator).
     */
    RCP<const BjtModel> model() const { return model_; }

    // ---- circuit temperature (`.TEMP`) ----
    // Thermal snapshot is read through a pointer into the card (like the
    // diode), so re-scaling the card is enough.
    RCP<const DeviceModel> temperatureCard() const override {
        return RCP<const DeviceModel>(model_);
    }

    /**
     * @brief Per-instance parameters.
     */
    const BjtInstanceParams& instanceParams() const { return inst_; }

    bool hasSubstrate() const { return nodes_.size() > 3; }

    /**
     * @brief Debug switch (diagnostics for convergence studies).
     */
    static void setDebug(bool on) { debug_ = on; }

private:
    RCP<const BjtModel> model_;    // shared card (never null)
    BjtInstanceParams   inst_;     // per-instance overrides

    // Thermal snapshot in force: the card's, or the instance override's.
    const BjtThermal* thermal_ = nullptr;
    BjtThermal ownThermal_;

    // Polarity multiplier: +1 NPN, −1 PNP.
    int p_ = 1;

    // ---- level-2 internal-node layout (offsets into the extra range) ----
    int numInternal_ = 0;
    int offB_ = -1, offC_ = -1, offE_ = -1;   // -1 ⇒ aliased to external
    // Pattern entry bases (positional indices into the pushed pattern).
    int patRb_ = -1, patRc_ = -1, patRe_ = -1;
    int patInt_ = 0;    // intrinsic 3x3 base (always present)
    int patSub_ = -1;   // substrate quad base (-1 when no substrate)

    // Conductances of the ohmic resistors (0 when absent).
    double gb_ = 0.0, gc_ = 0.0, ge_ = 0.0;

    // ---- last evaluation (all AREA-scaled, NPN-oriented) ----
    mutable double vbe_cached_ = 0.0;   // EXTERNAL terminal voltages
    mutable double vbc_cached_ = 0.0;
    mutable double vbe_lim_ = 0.0;      // intrinsic, limited, NPN-oriented
    mutable double vbc_lim_ = 0.0;
    mutable double ib_ = 0.0, ic_ = 0.0, ie_ = 0.0;
    mutable double qbe_ = 0.0, qbc_ = 0.0;
    mutable double isc_ = 0.0, gsc_ = 0.0, qcs_ = 0.0, csc_ = 0.0;

    // Jacobian ingredients (see the header derivation).
    mutable double gib_ = 0.0;    // ∂i_b/∂vbe
    mutable double gibc_ = 0.0;   // ∂i_b/∂vbc
    mutable double gm_ = 0.0;     // ∂i_c/∂vbe
    mutable double go_ = 0.0;     // ∂i_c/∂vbc
    mutable double cbe_v_ = 0.0;  // ∂Q_BE/∂vbe
    mutable double cbe_w_ = 0.0;  // ∂Q_BE/∂vbc
    mutable double cbc_u_ = 0.0;  // ∂Q_BC/∂vbe
    mutable double cbc_v_ = 0.0;  // ∂Q_BC/∂vbc

    // limiting bookkeeping
    mutable bool limiting_active_ = false;
    mutable bool first_solve_ = true;
    mutable double vbe_prev_ = 0.0;   // last limited junction voltages (PNJLIM)
    mutable double vbc_prev_ = 0.0;

    static bool debug_;

    // ---- internal/external row accessors ----
    int rowBi() const {
        return (offB_ >= 0) ? branch_index_ + offB_ : rowOf(nodes_[1]);
    }
    int rowCi() const {
        return (offC_ >= 0) ? branch_index_ + offC_ : rowOf(nodes_[0]);
    }
    int rowEi() const {
        return (offE_ >= 0) ? branch_index_ + offE_ : rowOf(nodes_[2]);
    }
    int rowSub() const {
        return hasSubstrate() ? rowOf(nodes_[3]) : -1;
    }

    /**
     * @brief Reads the raw junction voltages from the solution, applies the
     *        limiters and returns them through the argument list.
     */
    void limitedJunctionVoltages(const DaeSolution& x, bool allow_limiting) const;

    /**
     * @brief Evaluates the card laws at the (already applied) junction
     *        voltages: currents, charges and every Jacobian ingredient.
     * @param vbe, vbc  NPN-ORIENTED intrinsic junction voltages
     */
    void evaluate(double vbe, double vbc) const;

    /**
     * @brief Level-1 (Ebers-Moll) laws. Bit-identical to the pre-level-2
     *        implementation for LEVEL=1 cards.
     */
    void evaluateLevel1(double vbe, double vbc) const;

    /**
     * @brief Level-2 (Gummel-Poon) laws: qb, non-ideal base currents,
     *        substrate diode. With all GP extras off (IKF=IKR=ISE=ISC=0,
     *        no RB/RC/RE, no substrate) this reduces to level 1 up to the
     *        qb→q1 identity (icc = (fwd−rev)·ke, same ke floor).
     */
    void evaluateLevel2(double vbe, double vbc) const;

    /**
     * @brief SPICE3f5 pnjlim — see N_DEV_DeviceSupport.C in Xyce.
     */
    static double pnjlim(double vnew, double vold, double vt, double vcrit);
};

} // namespace e_XSpice

#endif // BJT_H

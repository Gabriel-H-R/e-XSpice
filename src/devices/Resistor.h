/**
 * @author Gabriel Henrique Silva
 * @file Resistor.h
 * @brief Linear resistor model (DAE formulation).
 *
 * Constitutive law:  V = I·R,  I = G·V,  G = 1/R.
 *
 * DAE formulation (resistor contributes only to F, never to Q):
 *
 *   F rows (KCL):
 *     F[i] += g·(v_i − v_j)
 *     F[j] −= g·(v_i − v_j)
 *
 *   J_F = dF/dx  (constant stamp):
 *       [ +g  −g ]
 *       [ −g  +g ]
 *
 * This device carries NO reactive state and NO internal unknowns.
 */

#ifndef RESISTOR_H
#define RESISTOR_H

#include "../core/Device.h"
#include "../core/Node.h"
#include "../models/core_models/PassiveModels.h"

namespace e_XSpice {

class Resistor : public Device {
public:
    /**
     * @brief Constructor with optional shared model card.
     * @param model Shared R card (multiplier + TC); null = ideal resistor.
     *              Effective R = R·Rmult·(1+TC1·ΔT+TC2·ΔT²), evaluated at the
     *              model's TNOM until temperature stepping exists (so TC
     *              terms are dormant and only Rmult applies today).
     */
    Resistor(const std::string& name,
             Node* node_plus,
             Node* node_minus,
             double resistance,
             RCP<const ResistorModel> model = RCP<const ResistorModel>());

    ~Resistor() override = default;

    // ---- DAE contract ----
    int numExtraUnknowns() const override { return 0; }

    void registerJacobianPattern(std::vector<int>& rows,
                                 std::vector<int>& cols) const override;

    void setJacobianSlots(const std::vector<int>& slots) override { jac_slots_ = slots; }

    void loadF(const DaeSolution& x, double time,
               std::vector<double>& F) const override;

    void loadJacobianF(DaeMatrix& J, const DaeSolution& x,
                       double time) const override;

    // ---- Outputs ----
    bool     isLinear() const override { return true; }
    double   getCurrent() const override;
    double   getPower()   const override;

    // Getters
    double getResistance() const { return resistance_; }
    double getConductance() const { return conductance_; }
    Node*  getNodePlus()   const { return nodes_[0]; }
    Node*  getNodeMinus()  const { return nodes_[1]; }
    RCP<const ResistorModel> model() const { return model_; }

    /**
     * @brief Voltage across the resistor (last accepted solution).
     */
    double getVoltage() const;

private:
    double resistance_;   // Nominal instance resistance (Ω)
    double conductance_;  // Effective conductance (S) = 1/(R·factor)
    RCP<const ResistorModel> model_;  // shared card (null = ideal)
};

} // namespace e_XSpice

#endif // RESISTOR_H

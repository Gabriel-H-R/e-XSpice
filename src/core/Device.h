/**
 * @author Gabriel Henrique Silva
 * @file Device.h
 * @brief Abstract base class for all circuit devices
 * 
 * This class defines the interface that all devices must implement.
 * Each device contributes to the conductance matrix (G) and current vector (b)
 * using the Modified Nodal Analysis (MNA) technique.
 */

#ifndef DEVICE_H
#define DEVICE_H

#include <string>
#include <vector>
#include <memory>
#include <complex>

namespace e_XSpice {

class Node;
class SparseMatrix;

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
 * @brief Base class for all devices
 */
class Device {
public:
    /**
     * @brief Constructor
     * @param name Device name (e.g., "R1", "C2")
     * @param type Device type
     */
    Device(const std::string& name, DeviceType type);
    
    /**
     * @brief Virtual destructor
     */
    virtual ~Device() = default;
    
    // Getters
    std::string getName() const { return name_; }
    DeviceType getType() const { return type_; }
    const std::vector<Node*>& getNodes() const { return nodes_; }
    
    /**
     * @brief Returns string with device type
     */
    std::string getTypeString() const;
    
    /**
     * @brief Stamping in the MNA system (DC or linearized Newton-Raphson)
     *
     * For LINEAR devices: always the same stamp (G and b are constant).
     * For NON-LINEAR devices (Diode, MOSFET...): calculates the linearized Norton equivalent
     * at the CURRENT operating point (reads node voltages):
     *
     * G_eq = dI/dV|_{V_k} (Jacobian)
     * I_eq = I(V_k) - G_eq * V_k (equivalent source)
     *
     * So stamp() directly implements the Newton-Raphson iteration:
     * J(x_k) * x_{k+1} = J(x_k)*x_k - F(x_k)
     *
     * @param G Conductance / Jacobian matrix (sparse)
     * @param b RHS vector
     * @param node_map Mapping node_id → index in the matrix
     */
    virtual void stamp(SparseMatrix& G, std::vector<double>& b,
                      const std::vector<int>& node_map) = 0;

    /**
     * @brief Stamping for AC analysis (small-signal, frequency ω)
     *
     * Inspired by Xyce's DAE formulation: devices provide their
     * complex admittance Y(jω) = G + jωC at the DC bias point.
     *
     * Linear devices:
     * - Resistor: Y = 1/R (real)
     * - Capacitor: Y = jωC (pure imaginary)
     * - Inductor: Y = 1/(jωL) (pure imaginary, negative)
     * - V-source: enforces KVL constraint with complex AC amplitude
     *
     * Non-linear devices:
     * - Uses G_d = dI/dV|_{V_dc} (small-signal at DC point)
     *
     * @param Z Complex AC system matrix (n x n)
     * @param b_ac Complex RHS vector (AC sources)
     * @param node_map Node ID → index mapping
     * @param omega Angular frequency ω = 2πf
     */
    virtual void stampAC(std::vector<std::vector<std::complex<double>>>& Z,
                 std::vector<std::complex<double>>& b_ac,
                 const std::vector<int>& node_map,
                 double omega) = 0;

    // ================================================================
    // DAE FORMULATION (Differential-Algebraic Equation)
    // ================================================================
    //
    // Xyce solves: C(x)·ẋ + F(x) = b(t)
    //
    // F(x) = current vector (Kirchhoff laws + device equations)
    // Q(x) = charge/flux vector (capacitors, inductors)
    // C(x) = dQ/dx (dynamic Jacobian)
    // J(x) = dF/dx (static Jacobian — used in NR)
    // b(t) = independent sources (voltage/current)
    // ẋ = time derivative of state variables (voltages, currents)
    //
    // Each device is responsible for knowing its own contribution.
    // The solver (NR, Transient) calls the appropriate methods.
    // ================================================================
    /**
     * @brief Loads F(x) — static KCL residual
     * @brief Current/voltage contribution of the device to the DAE residual
     * @brief vector F. Linear devices: F = G*x - b. Nonlinear devices: F = I(x) (current as function of voltages).
     * @param F Residual vector (output, accumulated)
     * @param node_map Mapping node_id → index
     */
    virtual void loadDaeF(std::vector<double>& F,
                          const std::vector<int>& node_map) const {
                            (void)F; (void)node_map;
        }

    /**
     * @brief Loads Q(x) — charges and fluxes (dynamic part)
     * @brief Capacitors contribute: Q_i = C*(V_i - V_j)
     * @brief Inductors contribute: Q_i = L*I_i (through current variable)
     * @brief Resistors and DC sources: Q = 0
     * @param Q Charge/flux vector (output)
     * @param node_map Mapping node_id → index
     */
    virtual void loadDaeQ(std::vector<double>& Q,
                          const std::vector<int>& node_map) const {
                            (void)Q; (void)node_map;
        }

    /**
     * @brief Loads Jacobian J = dF/dx
     * @brief For Newton-Raphson: J * Δx = -F(x)
     * @brief Linear devices: J = G (constant)
     * @brief Nonlinear devices: J = ∂I/∂V (e.g. diode G_d)
     * @param J Sparse Jacobian (output, accumulated)
     * @param node_map Mapping node_id → index
     */
    virtual void loadJacobianF(SparseMatrix& J,
                               const std::vector<int>& node_map) {
                                (void)J; (void)node_map;
        // Devices that use stamp() as the Jacobian do not need
        // to reimplement this method separately.
        }

    /**
     * @brief Loads Jacobian dQ/dx (for transient integration)
     * @brief Required by implicit integration methods that need
     * @brief d(C·ẋ)/dx. Capacitor: dQ/dV = C. Inductor: dQ/dI = L.
     * @param J Sparse Jacobian (output)
     * @param node_map Mapping node_id → index
     */
    virtual void loadJacobianQ(SparseMatrix& J,
                               const std::vector<int>& node_map) const {
                                (void)J; (void)node_map;
        }

    /**
     * @brief Updates device state based on node voltages after each step (NR or transient)
     * 
     * Used for non-linear devices that need to update
     * their contributions based on the current solution.
     */
    virtual void updateState() {}
    
    /**
     * @brief Returns if the device is linear
     */
    virtual bool isLinear() const = 0;
    
    /**
     * @brief Calculates current through the device
     */
    virtual double getCurrent() const = 0;
    
    /**
     * @brief Calculates dissipated power
     */
    virtual double getPower() const = 0;

protected:
    std::string name_;           // Unique device name
    DeviceType type_;           // Device type
    std::vector<Node*> nodes_;  // Nodes connected to the device
};

} // namespace e_XSpice

#endif // DEVICE_H

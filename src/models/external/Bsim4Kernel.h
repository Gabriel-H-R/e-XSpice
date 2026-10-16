/**
 * @author Gabriel Henrique Silva
 * @file Bsim4Kernel.h
 * @brief Narrow adapter interface between e-XSpice DAE and BSIM4 physics.
 *
 * Direction of dependency (deliberate, and the key difference from Xyce):
 * the external C code adapts to e-XSpice, never the reverse. The kernel
 * receives plain doubles plus the e-XSpice card/thermal snapshots and
 * returns plain doubles. No CKTcircuit, no GENmodel, no IFuid, no
 * SMPmatrix, no global state — so the kernel is thread-safe and
 * MPI/OpenMP-safe by construction.
 *
 * Xyce reference (verified in-repo layout, Oct 2026):
 *   src/DeviceModelPKG/OpenModels/N_DEV_MOSFET_B4.{C,h} plus versioned
 *   variants B4p61/B4p70/B4p82, with ADMS/Verilog-A under
 *   src/DeviceModelPKG/ADMS. Xyce wraps versioned Berkeley-derived code
 *   inside per-version device classes. We keep the versioning idea but
 *   invert the coupling: one DAE evaluator (MosfetBSIM4Device) + one
 *   kernel per physics source, selected at setup time.
 */
#ifndef BSIM4_KERNEL_H
#define BSIM4_KERNEL_H

#include "../core_models/MosfetBSIM4Model.h"
#include <string>

namespace e_XSpice {

/**
 * @brief Plain inputs to any BSIM4 physics kernel (no solver types).
 */
struct Bsim4KernelInputs {
    const MosfetBSIM4Model*   card = nullptr;    // card parameters (setup)
    const MosfetBSIM4Thermal* thermal = nullptr; // temperature snapshot
    double W = 0.0, L = 0.0, M = 1.0, AREA = 1.0;
    double vgs = 0.0, vds = 0.0, vbs = 0.0;      // N-oriented voltages
};

/**
 * @brief Plain outputs of any BSIM4 physics kernel (N-oriented).
 */
struct Bsim4KernelOutputs {
    double id = 0.0, is = 0.0, ib = 0.0;
    double gm = 0.0, gds = 0.0, gmbs = 0.0;
    double qg = 0.0, qd = 0.0, qb = 0.0;
    double cgs = 0.0, cgd = 0.0, cgb = 0.0;
    double cdg = 0.0, cdd = 0.0, cdb = 0.0;
    double csg = 0.0, csd = 0.0, css = 0.0, csb = 0.0;
    double cbg = 0.0, cbd = 0.0, cbs = 0.0, cbb = 0.0;
    double vth = 0.0;
    bool saturated = false;
    bool swapped = false;
};

/**
 * @brief Physics source behind MosfetBSIM4Device::evaluate().
 */
class IBsim4Kernel : public RefCounted {
public:
    ~IBsim4Kernel() override = default;
    virtual const char* name() const = 0;
    virtual bool available(std::string& whyNot) const = 0;
    virtual bool eval(const Bsim4KernelInputs& in,
                      Bsim4KernelOutputs& out,
                      std::string& err) const = 0;
};

} // namespace e_XSpice

#endif // BSIM4_KERNEL_H

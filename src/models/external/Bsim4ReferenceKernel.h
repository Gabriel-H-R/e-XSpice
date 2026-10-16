/**
 * @author Gabriel Henrique Silva
 * @file Bsim4ReferenceKernel.h
 * @brief Simplified BSIM4-style kernel in pure C++ (no Berkeley headers).
 *
 * This is the physics previously inlined in MosfetBSIM4Device::evaluate(),
 * moved verbatim behind the IBsim4Kernel interface. It is a deliberately
 * SMALL subset (threshold + mobility degradation + Early + Meyer + bulk
 * diodes) — sufficient for regression and smoke tests, explicitly NOT a
 * claim of BSIM4 4.8.3 conformance (no NQS, gate leakage, RF, stress).
 */
#ifndef BSIM4_REFERENCE_KERNEL_H
#define BSIM4_REFERENCE_KERNEL_H

#include "Bsim4Kernel.h"

namespace e_XSpice {

class Bsim4ReferenceKernel : public IBsim4Kernel {
public:
    const char* name() const override { return "reference"; }
    bool available(std::string& whyNot) const override {
        (void)whyNot;
        return true;
    }
    bool eval(const Bsim4KernelInputs& in,
              Bsim4KernelOutputs& out,
              std::string& err) const override;
};

} // namespace e_XSpice

#endif // BSIM4_REFERENCE_KERNEL_H

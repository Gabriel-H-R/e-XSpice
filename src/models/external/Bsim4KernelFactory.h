/**
 * @author Gabriel Henrique Silva
 * @file Bsim4KernelFactory.h
 * @brief Setup-time kernel selector (no per-iteration cost).
 *
 * "reference" (default): pure-C++ subset, always available, used by the
 * committed regression netlist. "berkeley": upstream b4ld.c/b4temp.c via
 * the adapter skeleton; reports unavailable until the shim is vendored.
 * "auto": berkeley when compiled in, else reference.
 */
#ifndef BSIM4_KERNEL_FACTORY_H
#define BSIM4_KERNEL_FACTORY_H

#include "Bsim4Kernel.h"
#include "Bsim4ReferenceKernel.h"
#include "Bsim4BerkeleyKernel.h"

namespace e_XSpice {

class Bsim4KernelFactory {
public:
    static inline RCP<IBsim4Kernel> create(const std::string& which) {
        if (which == "berkeley") {
            return make_rcp<Bsim4BerkeleyKernel>();
        }
        if (which == "auto") {
            std::string why;
            RCP<IBsim4Kernel> b = make_rcp<Bsim4BerkeleyKernel>();
            if (b->available(why)) return b;
            return make_rcp<Bsim4ReferenceKernel>();
        }
        return make_rcp<Bsim4ReferenceKernel>();
    }

    static inline const char* activeDefault() {
#ifdef EXSPICE_USE_BERKELEY_BSIM4
        return "berkeley (compiled in)";
#else
        return "reference (Berkeley shim absent)";
#endif
    }
};

} // namespace e_XSpice

#endif // BSIM4_KERNEL_FACTORY_H

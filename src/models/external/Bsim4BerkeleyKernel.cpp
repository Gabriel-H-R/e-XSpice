/**
 * @author Gabriel Henrique Silva
 * @file Bsim4BerkeleyKernel.cpp
 * @brief Berkeley binding point (disabled until shim headers are vendored).
 */
#include "Bsim4BerkeleyKernel.h"

namespace e_XSpice {

bool Bsim4BerkeleyKernel::available(std::string& whyNot) const {
#ifdef EXSPICE_USE_BERKELEY_BSIM4
    (void)whyNot;
    return true;
#else
    whyNot = "Berkeley BSIM4 C kernel not compiled: "
             "src/models/bsim4/code needs spice.h/cktdefs.h/devdefs.h "
             "(absent). Vendor the ngspice/SPICE3 shim under "
             "src/models/external/berkeley_shim/ and rebuild with "
             "-DEXSPICE_USE_BERKELEY_BSIM4=ON. Binding points: "
             "BSIM4temp (b4temp.c) at setup, BSIM4load (b4ld.c) per "
             "iteration. Falling back to the reference kernel.";
    return false;
#endif
}

bool Bsim4BerkeleyKernel::eval(const Bsim4KernelInputs& in,
                               Bsim4KernelOutputs& out,
                               std::string& err) const {
    (void)in;
    (void)out;
#ifdef EXSPICE_USE_BERKELEY_BSIM4
    err = "Berkeley kernel selected but per-instance BSIM4load bridge is "
          "not yet wired to sBSIM4instance/sBSIM4model; see header contract.";
    return false;
#else
    err = "Berkeley kernel unavailable (see available()). "
          "Use kernel \"reference\".";
    return false;
#endif
}

} // namespace e_XSpice

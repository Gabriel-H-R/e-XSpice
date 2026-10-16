/**
 * @author Gabriel Henrique Silva
 * @file Bsim4BerkeleyKernel.h
 * @brief Adapter skeleton for the upstream Berkeley BSIM4 4.8.3 C code.
 *
 * Vendored source (untouched, ECL-2.0) lives in src/models/bsim4/code/:
 * b4ld.c (BSIM4load), b4temp.c (BSIM4temp), b4set.c, b4mpar.c, bsim4def.h.
 * Direct compilation fails by design — verified in-repo:
 *   b4ld.c:25        #include "spice.h"    (absent)
 *   b4ld.c:29-35     cktdefs.h, devdefs.h, trandefs.h, const.h, sperror.h
 *   bsim4def.h:28-32 ifsim.h, gendefs.h, cktdefs.h, complex.h, noisedef.h
 *   bsim4def.h:34     struct sBSIM4instance { BSIM4modPtr, nextInstance,
 *                     IFuid name, states, dNode/gNode/sNode/bNode, ... }
 * i.e. the C code is coupled to SPICE3 data structures (CKTcircuit,
 * GENmodel/GENinstance, IFuid, SMPmatrix), not to standalone equations.
 *
 * Binding contract (when EXSPICE_USE_BERKELEY_BSIM4=ON with shim headers):
 *   setup:  BSIM4mpar/b4set.c  (card params) + BSIM4temp/b4temp.c (T scaling)
 *   eval:   BSIM4load/b4ld.c    (per-iteration I/Q + gm/gds/gmbs/C)
 * The adapter translates e-XSpice POD <-> Berkeley structs at these two
 * points only. Until the shim exists, available() is false with a message.
 */
#ifndef BSIM4_BERKELEY_KERNEL_H
#define BSIM4_BERKELEY_KERNEL_H

#include "Bsim4Kernel.h"

namespace e_XSpice {

class Bsim4BerkeleyKernel : public IBsim4Kernel {
public:
    const char* name() const override { return "berkeley"; }
    bool available(std::string& whyNot) const override;
    bool eval(const Bsim4KernelInputs& in,
              Bsim4KernelOutputs& out,
              std::string& err) const override;
};

} // namespace e_XSpice

#endif // BSIM4_BERKELEY_KERNEL_H

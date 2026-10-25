/**
 * @author Gabriel Henrique Silva
 * @file SpiceHeaders.h
 * @brief Umbrella header — includes all SPICE compatibility stubs.
 *
 * Each individual wrapper (spice.h, cktdefs.h, devdefs.h, smpdefs.h, ...)
 * is self-contained. This umbrella header includes them all for convenience
 * in C++ code.
 */

#ifndef EXSPICE_SPICE_HEADERS_H
#define EXSPICE_SPICE_HEADERS_H

#ifdef __cplusplus
#include <vector>
extern "C" {
#endif

#include "spice.h"
#include "cktdefs.h"
#include "devdefs.h"
#include "smpdefs.h"
#include "ifsim.h"
#include "sperror.h"
#include "const.h"
#include "trandefs.h"
#include "noisedef.h"
#include "util.h"
#include "suffix.h"
#include "complex.h"
#include "gendefs.h"
#include "jobdefs.h"
#include "fteconst.h"
#include "ftedefs.h"
#include "fteext.h"
#include "acdefs.h"
#include "inpdefs.h"
#include "inpmacs.h"
#include "misc.h"
#include "strext.h"
#include "isrcdefs.h"
#include "vsrcdefs.h"

/* NIintegrate: SPICE integration routine (stub — level-2 BSIM4 subset
 * does not model charges, so it returns 0 without doing anything). */
int NIintegrate(void *ckt, double *geq, double *ceq,
                double cf, int state);

/* DEVpnjlim: SPICE PNJLIM diode limiting function (stub). */
int DEVpnjlim(double vnew, double vold, double vt, double vcrit,
              double *vlimited);

/* BSIM4PAeffGeo / BSIM4RdseffGeo: geometry effect functions (stubs —
 * level-2 subset doesn't model geometry/proximity effects). */
int BSIM4PAeffGeo(int nf, int geoMod, int min, double *pseff,
                  double *adeff, double *pdeff);
int BSIM4RdseffGeo(int nf, int geoMod, int min, double *wr,
                   double *dwg, double *dwb, double *prwb, double *prwg);

/* BSIM4checkModel: model consistency check (stub — always OK).
 * Declared in bsim4ext.h (included via bsim4def.h). */

/* SPfrontEnd: SPICE front-end global pointer (stub).
 * Defined in spice.h and BSim4Adapter.cpp. */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* EXSPICE_SPICE_HEADERS_H */

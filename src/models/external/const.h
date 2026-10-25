#ifndef EXSPICE_CONST_H
#define EXSPICE_CONST_H
#include "spice.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define Charge_q     1.60219e-19   /* elementary charge (C) */
#define Boltzmann    1.3806226e-23 /* Boltzmann constant (J/K) */
#define KboQ         8.617087e-5   /* k/q (V/K) */
#define EPS0         8.85418e-12   /* permittivity of vacuum (F/m) */
#define EPSSI        1.03594e-10   /* permittivity of silicon (F/m) */

#ifndef CHARGE
#define CHARGE       Charge_q
#endif

#endif /* EXSPICE_CONST_H */

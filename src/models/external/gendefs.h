#ifndef EXSPICE_GENDEFS_H
#define EXSPICE_GENDEFS_H
#include "spice.h"

#define NMOS   1
#define PMOS  -1

/* Number of BSIM4 state variables per instance */
#ifndef NSTATVARS
#define NSTATVARS 30
#endif

/* NUMELEMS: array size macro used by the Berkeley parameter tables. */
#ifndef NUMELEMS
#define NUMELEMS(a) (sizeof(a) / sizeof((a)[0]))
#endif

/* DEV_DEFAULT alias for compatibility */
#ifndef DEV_DEFAULT
#define DEV_DEFAULT 0
#endif

#endif /* EXSPICE_GENDEFS_H */

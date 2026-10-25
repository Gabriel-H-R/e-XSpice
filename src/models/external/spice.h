/**
 * @file spice.h
 * @brief SPICE core types — master header for the Berkeley BSIM4 C code.
 *
 * Defines ALL types the Berkeley BSIM4 C code needs in a single header,
 * avoiding any circular include issues:
 *
 *   - Basic types: REAL, BOOLEAN, IFuid, SPcomplex
 *   - Mode flags: MODEDCOP, MODETRAN, ...
 *   - IF types: IF_REAL, IF_INTEGER, ...
 *   - Device types: GENmodel, GENinstance
 *   - Circuit types: CKTnode, CKTcircuit (with extra fields BSIM4load reads)
 *   - Matrix type: SMPmatrix
 *   - Parameter types: IFparm, IFvalue
 *
 * Individual wrapper headers (devdefs.h, cktdefs.h, smpdefs.h, ...) simply
 * include spice.h for backward compatibility with the Berkeley #include
 * chain. Include guards prevent circular re-entry.
 */

#ifndef EXSPICE_SPICE_H
#define EXSPICE_SPICE_H

/* =====================================================================
 *  Basic types
 * ===================================================================== */

typedef double  REAL;
typedef int     BOOLEAN;
typedef char   *IFuid;

typedef struct  { double r, i; } SPcomplex;

/* =====================================================================
 *  Mode flags (CKTcircuit->CKTmode)
 * ===================================================================== */

#define MODEDCOP          0x0001
#define MODETRAN          0x0002
#define MODEAC            0x0004
#define MODEINITSMSIG     0x0008
#define MODEINITTRAN      0x0010
#define MODEINITJCT       0x0020
#define MODEINITFIX       0x0040
#define MODEUIC           0x0080
#define MODEDCTRANCURVE   0x0100
#define MODETRANOP        0x0200
#define MODEPREDICTOR     0x0400
#define MODEINITPRED      0x0800

#define DEV_DEFAULT       0

/* =====================================================================
 *  SPICE utility macros
 * ===================================================================== */

#ifndef MAX
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#endif
#ifndef FABS
#define FABS(x)   ((x) < 0.0 ? -(x) : (x))
#endif

/* =====================================================================
 *  IF parameter types
 * ===================================================================== */

#define IF_REAL      0
#define IF_INTEGER   1
#define IF_FLAG      2
#define IF_STRING    3
#define IF_NODE      4
#define IF_INSTANCE  5
#define IF_REALVEC   6
#define IF_COMPLEX   7

#define IF_SET       0x01
#define IF_ASK       0x02
#define IF_NEW       0x04
#define IF_NOCHECK   0x08

/* =====================================================================
 *  Device types (GENmodel, GENinstance)
 * ===================================================================== */

struct GENinstance {
    int             GENnumNodes;   /* number of nodes */
    int             GENstate;      /* base state index for this instance */
    int            *GENnode;       /* node numbers (into CKTnodes) */
    struct GENmodel *GENmodel;     /* back-pointer to the model card */
};

struct GENmodel {
    IFuid       GENmodName;   /* model name */
    int         GENmodType;   /* device type */
    void       *GENnextModel; /* next model in list */
    void       *GENinstances; /* head of instance list */
};

/* C requires typedef for GENmodel/GENinstance without 'struct' keyword */
typedef struct GENinstance GENinstance;
typedef struct GENmodel    GENmodel;

/* =====================================================================
 *  Circuit types (CKTnode, CKTcircuit)
 * ===================================================================== */

struct SMPmatrix;

struct CKTnode {
    int     number;     /* matrix row index (-1 for ground) */
    IFuid   name;       /* node name */
};

struct CKTcircuit {
    double    CKTtemp;          /* circuit temperature (Kelvin) */
    double    CKTgmin;          /* parallel conductance */
    int       CKTmode;          /* analysis mode flags */
    double   *CKTrhs;           /* right-hand side vector */
    double   *CKTrhsOld;        /* previous RHS vector */
    double   *CKTirhsOld;       /* previous RHS vector (implicit) */
    double   *CKTstate0;        /* solution vector (current iteration) */
    double   *CKTstate1;        /* previous solution / state */
    double   *CKTstate2;        /* second previous solution */
    double   *CKTstate3;        /* third previous solution */
    struct SMPmatrix *CKTmatrix; /* sparse matrix */
    int       CKTnumNodes;      /* number of non-ground nodes */
    struct CKTnode **CKTnodes;  /* node list */
    /* Integration / time-step fields: */
    double    CKTdelta;         /* current time step */
    double    CKTdeltaOld[8];   /* previous time steps */
    double    CKTag[4];         /* integration coefficients */
    /* Convergence fields: */
    int       CKTnoncon;        /* non-convergence counter */
    double    CKTreltol;        /* relative tolerance */
    double    CKTabstol;        /* absolute tolerance */
    /* Job / iteration pointers: */
    void     *CKTcurJob;        /* current job */
    void     *CKTcurPtr;        /* current pointer */
    void     *CKTorder;         /* ordering */
    double    CKTstep;          /* time step */
    int       CKTbypass;        /* bypass flag */
    double    CKTvoltTol;       /* voltage tolerance */
};

/* C requires typedef for CKTcircuit/CKTnode/SMPmatrix without 'struct' keyword */
typedef struct CKTnode      CKTnode;
typedef struct CKTcircuit   CKTcircuit;
typedef struct SMPmatrix    SMPmatrix;

/* =====================================================================
 *  Matrix type (SMPmatrix)
 * ===================================================================== */

struct SMPmatrix {
    int size;               /* dimension of the matrix */
#ifdef __cplusplus
    // C++-only fields used by the SMPaddElt shim:
    #include <vector>
    std::vector<int> slot_of;   // slot_of[row * size + col] (or -1)
    double* values = nullptr;   // borrowed DaeMatrix value array
#endif
};

/* =====================================================================
 *  IFparm / IFvalue structs
 * ===================================================================== */

typedef struct sIFparm   IFparm;
typedef struct sIFvalue  IFvalue;

struct sIFparm {
    IFuid   name;
    int     id;
    int     type;
    int     flags;
    double  defaultValue;
    IFuid   stringDefault;
};

struct sIFvalue {
    int     type;
    union {
        double  rValue;
        int     iValue;
        IFuid   sValue;
        struct { int nElements; double *rVector; } vValue;
        SPcomplex cValue;
        void   *pValue;
    };
};

/* =====================================================================
 *  Noise analysis types
 * ===================================================================== */

struct Ndata {
    int    dummy;
};
typedef struct Ndata Ndata;

/* =====================================================================
 *  SPICE error codes
 * ===================================================================== */

#define OK          0
#define E_BADPARM   1
#define E_NOMEM     2
#define E_BADCHAR   3
#define E_CYCLE     4
#define E_NODE      5
#define E_TRANS     6
#define E_UNSUPP    7
#define E_FATAL     8
#define E_TIME      9
#define E_CONVERG   10

#define ERR_FATAL   E_FATAL
#define ERR_WARNING 11
#define ERR_INFO    12

/* =====================================================================
 *  SPICE physical constants (from const.h)
 * ===================================================================== */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define Charge_q     1.60219e-19   /* elementary charge (C) */
#define Boltzmann    1.3806226e-23 /* Boltzmann constant (J/K) */
#define KboQ         8.617087e-5   /* k/q (V/K) */
#define EPS0         8.85418e-12   /* permittivity of vacuum (F/m) */
#define EPSSI        1.03594e-10   /* permittivity of silicon (F/m) */

#define CONSTvt0     0.025852      /* kT/q at 300.15K */
#define CONSTroot2   1.414213562   /* sqrt(2) */

/* =====================================================================
 *  SPICE front end (global context pointer)
 * ===================================================================== */

typedef int (*IFerrorFunc)(int, ...);

struct SPfrontEnd {
    IFerrorFunc *IFerror;
    void       *dummy;
};
extern struct SPfrontEnd *SPfrontEnd;

#endif /* EXSPICE_SPICE_H */

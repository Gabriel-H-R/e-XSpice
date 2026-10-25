#ifndef EXSPICE_JOBDEFS_H
#define EXSPICE_JOBDEFS_H
#include "spice.h"

struct JOB {
    char   *JOBname;
    struct JOB *JOBnextJob;
};

struct TSKtask {
    struct {
        struct JOB *jobs;
    } *ci_curTask;
};

#define ft_curckt ((struct { struct TSKtask *ci_curTask; } *)0)

#endif /* EXSPICE_JOBDEFS_H */

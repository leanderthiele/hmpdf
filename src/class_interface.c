#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <class.h>

#include "data.h"
#include "class_interface.h"

#include "hmpdf.h"

static
void alloc_class(all_data *d)
{//{{{
    d->cls->pr = malloc(sizeof(struct precision));
    d->cls->ba = malloc(sizeof(struct background));
    d->cls->th = malloc(sizeof(struct thermo));
    d->cls->pt = malloc(sizeof(struct perturbs));
    d->cls->tr = malloc(sizeof(struct transfers));
    d->cls->pm = malloc(sizeof(struct primordial));
    d->cls->sp = malloc(sizeof(struct spectra));
    d->cls->nl = malloc(sizeof(struct nonlinear));
    d->cls->le = malloc(sizeof(struct lensing));
    d->cls->op = malloc(sizeof(struct output));
}//}}}

static
void run_class(all_data *d)
{//{{{
    fprintf(stdout, "\trun_class\n");
    fflush(stdout);

    struct precision *pr = (struct precision *)d->cls->pr;
    struct background *ba = (struct background *)d->cls->ba;
    struct thermo *th = (struct thermo *)d->cls->th;
    struct perturbs *pt = (struct perturbs *)d->cls->pt;
    struct primordial *pm = (struct primordial *)d->cls->pm;
    struct nonlinear *nl = (struct nonlinear *)d->cls->nl;

    fprintf(stdout, "\t\tbackground\n");
    fflush(stdout);
    SAFECLASS(background_init(pr, ba), ba->error_message)
    fprintf(stdout, "\t\tthermodynamics\n");
    fflush(stdout);
    SAFECLASS(thermodynamics_init(pr, ba, th), th->error_message)
    fprintf(stdout, "\t\tperturbs\n");
    fflush(stdout);
    SAFECLASS(perturb_init(pr, ba, th, pt), pt->error_message)
    fprintf(stdout, "\t\tprimordial\n");
    fflush(stdout);
    SAFECLASS(primordial_init(pr, pt, pm), pm->error_message)
    fprintf(stdout, "\t\tnonlinear\n");
    fflush(stdout);
    SAFECLASS(nonlinear_init(pr, ba, th, pt, pm, nl), nl->error_message)
}//}}}

void init_class_interface(all_data *d)
{//{{{
    fprintf(stdout, "In class_interface.h -> init_class.\n");
    fflush(stdout);
    char **argv = (char **)malloc(3 * sizeof(char *));
    argv[1] = d->cls->class_ini;
    argv[2] = d->cls->class_pre;

    int argc = (strcmp(argv[2], "none")) ? 3 : 2;

    alloc_class(d);

    struct precision *pr = (struct precision *)d->cls->pr;
    struct background *ba = (struct background *)d->cls->ba;
    struct thermo *th = (struct thermo *)d->cls->th;
    struct perturbs *pt = (struct perturbs *)d->cls->pt;
    struct primordial *pm = (struct primordial *)d->cls->pm;
    struct nonlinear *nl = (struct nonlinear *)d->cls->nl;
    struct spectra *sp = (struct spectra *)d->cls->sp;
    struct lensing *le = (struct lensing *)d->cls->le;
    struct output *op = (struct output *)d->cls->op;
    struct transfers *tr = (struct transfers *)d->cls->tr;

    ErrorMsg errmsg;
    SAFECLASS(input_init_from_arguments(argc, argv, pr, ba, th,
                                       pt, tr, pm, sp,
                                       nl, le, op, errmsg),
              errmsg)

    run_class(d);

    free(argv);
}//}}}

void null_class_interface(all_data *d)
{//{{{
    d->cls->pr = NULL;
    d->cls->ba = NULL;
    d->cls->th = NULL;
    d->cls->pt = NULL;
    d->cls->pm = NULL;
    d->cls->nl = NULL;
    d->cls->sp = NULL;
    d->cls->le = NULL;
    d->cls->op = NULL;
    d->cls->tr = NULL;
}//}}}

void reset_class_interface(all_data *d)
{//{{{
    if (d->cls->op != NULL) { free(d->cls->op); }
    if (d->cls->le != NULL) { free(d->cls->le); }
    if (d->cls->sp != NULL) { free(d->cls->sp); }
    if (d->cls->tr != NULL) { free(d->cls->tr); }
    if (d->cls->nl != NULL) { nonlinear_free(d->cls->nl); free(d->cls->nl); }
    if (d->cls->pt != NULL) { perturb_free(d->cls->pt); free(d->cls->pt); }
    if (d->cls->pm != NULL) { primordial_free(d->cls->pm); free(d->cls->pm); }
    if (d->cls->th != NULL) { thermodynamics_free(d->cls->th); free(d->cls->th); }
    if (d->cls->ba != NULL) { background_free(d->cls->ba); free(d->cls->ba); }
    if (d->cls->pr != NULL) { free(d->cls->pr); }
}//}}}

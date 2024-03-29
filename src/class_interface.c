#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <class.h>

#include "utils.h"
#include "object.h"
#include "class_interface.h"

#include "hmpdf.h"

static int
alloc_class(hmpdf_obj *d)
{//{{{
    STARTFCT

    SAFEALLOC(d->cls->pr, malloc(sizeof(struct precision)));
    SAFEALLOC(d->cls->ba, malloc(sizeof(struct background)));
    SAFEALLOC(d->cls->th, malloc(sizeof(struct thermodynamics)));
    SAFEALLOC(d->cls->pt, malloc(sizeof(struct perturbations)));
    SAFEALLOC(d->cls->tr, malloc(sizeof(struct transfer)));
    SAFEALLOC(d->cls->pm, malloc(sizeof(struct primordial)));
    SAFEALLOC(d->cls->sp, malloc(sizeof(struct harmonic)));
    SAFEALLOC(d->cls->nl, malloc(sizeof(struct fourier)));
    SAFEALLOC(d->cls->le, malloc(sizeof(struct lensing)));
    SAFEALLOC(d->cls->sd, malloc(sizeof(struct distortions)));
    SAFEALLOC(d->cls->op, malloc(sizeof(struct output)));

    ENDFCT
}//}}}

static int
run_class(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\trun_class\n");

    struct precision *pr = (struct precision *)d->cls->pr;
    struct background *ba = (struct background *)d->cls->ba;
    struct thermodynamics *th = (struct thermodynamics *)d->cls->th;
    struct perturbations *pt = (struct perturbations *)d->cls->pt;
    struct primordial *pm = (struct primordial *)d->cls->pm;
    struct fourier *nl = (struct fourier *)d->cls->nl;

    HMPDFPRINT(3, "\t\tbackground\n");
    SAFECLASS(background_init(pr, ba), ba->error_message);

    // check if user passed curved geometry
    if (ba->sgnK)
    {
        HMPDFWARN("The code does not reliably support non-flat universes "
                  "at the moment.");
    }

    HMPDFPRINT(3, "\t\tthermodynamics\n");
    SAFECLASS(thermodynamics_init(pr, ba, th), th->error_message);
    
    HMPDFPRINT(3, "\t\tperturbs\n");
    SAFECLASS(perturbations_init(pr, ba, th, pt), pt->error_message);

    // check if user passed a correct CLASS .ini file
    if (pt->has_pk_matter != _TRUE_)
    {
        HMPDFWARN("You should set output=mPk in the CLASS .ini file.");
    }

    // check if user passed large enough k_max
    if (pt->k_max_for_pk < 1.0)
    {
        HMPDFWARN("You should set P_k_max_h/Mpc or P_k_max_1/Mpc "
                  "to > 1 h/Mpc in the CLASS .ini file.");
    }
    
    HMPDFPRINT(3, "\t\tprimordial\n");
    SAFECLASS(primordial_init(pr, pt, pm), pm->error_message);

    HMPDFPRINT(3, "\t\tnonlinear\n");
    SAFECLASS(fourier_init(pr, ba, th, pt, pm, nl), nl->error_message);

    ENDFCT
}//}}}

int
init_class_interface(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tinit_class_interface\n");

    char **argv;
    SAFEALLOC(argv, malloc(3 * sizeof(char *)));
    argv[1] = d->cls->class_ini;
    argv[2] = d->cls->class_pre;

    int argc = (strcmp(argv[2], "none")) ? 3 : 2;

    SAFEHMPDF(alloc_class(d));

    struct precision *pr = (struct precision *)d->cls->pr;
    struct background *ba = (struct background *)d->cls->ba;
    struct thermodynamics *th = (struct thermodynamics *)d->cls->th;
    struct perturbations *pt = (struct perturbations *)d->cls->pt;
    struct primordial *pm = (struct primordial *)d->cls->pm;
    struct fourier *nl = (struct fourier *)d->cls->nl;
    struct harmonic *sp = (struct harmonic *)d->cls->sp;
    struct lensing *le = (struct lensing *)d->cls->le;
    struct output *op = (struct output *)d->cls->op;
    struct transfer *tr = (struct transfer *)d->cls->tr;
    struct distortions *sd = (struct distortions *)d->cls->sd;

    ErrorMsg errmsg;
    SAFECLASS(input_init(argc, argv, pr, ba, th,
                         pt, tr, pm, sp,
                         nl, le, sd, op ,errmsg),
              errmsg);

    SAFEHMPDF(run_class(d));

    free(argv);

    ENDFCT
}//}}}

int
null_class_interface(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->cls->pr = NULL;
    d->cls->ba = NULL;
    d->cls->th = NULL;
    d->cls->pt = NULL;
    d->cls->pm = NULL;
    d->cls->nl = NULL;
    d->cls->sp = NULL;
    d->cls->le = NULL;
    d->cls->tr = NULL;
    d->cls->sd = NULL;
    d->cls->op = NULL;

    ENDFCT
}//}}}

int
reset_class_interface(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_class_interface\n");

    struct fourier *nl;
    struct perturbations *pt;
    struct primordial *pm;
    struct thermodynamics *th;
    struct background *ba;

    if (d->cls->op != NULL) { free(d->cls->op); }
    if (d->cls->le != NULL) { free(d->cls->le); }
    if (d->cls->sp != NULL) { free(d->cls->sp); }
    if (d->cls->tr != NULL) { free(d->cls->tr); }
    if (d->cls->sd != NULL) { free(d->cls->sd); }
    if (d->cls->nl != NULL)
    { 
        nl = (struct fourier *)d->cls->nl;
        SAFECLASS(fourier_free(nl),
                  nl->error_message);
        free(d->cls->nl);
    }
    if (d->cls->pt != NULL)
    {
        pt = (struct perturbations *)d->cls->pt;
        SAFECLASS(perturbations_free(pt),
                  pt->error_message);
        free(d->cls->pt);
    }
    if (d->cls->pm != NULL)
    {
        pm = (struct primordial *)d->cls->pm;
        SAFECLASS(primordial_free(pm),
                  pm->error_message);
        free(d->cls->pm);
    }
    if (d->cls->th != NULL)
    {
        th = (struct thermodynamics *)d->cls->th;
        SAFECLASS(thermodynamics_free(th),
                  th->error_message);
        free(d->cls->th);
    }
    if (d->cls->ba != NULL)
    {
        ba = (struct background *)d->cls->ba;
        SAFECLASS(background_free(ba),
                  ba->error_message);
        free(d->cls->ba); }
    if (d->cls->pr != NULL) { free(d->cls->pr); }

    ENDFCT
}//}}}

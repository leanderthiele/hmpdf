#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <class.h>

#include "data.h"
#include "class_interface.h"

struct class_interface_s
{//{{{
    char *class_ini;
    char *class_pre;

    struct precision *pr;
    struct background *ba;
    struct thermo *th;
    struct primordial *pm;
    struct perturbs *pt;
    struct nonlinear *nl;
    struct transfers *tr;
    struct spectra *sp;
    struct lensing *le;
    struct output *op;

    ErrorMsg errmsg;
};//}}}

typedef struct class_interface_s cls;


static
void alloc_class(all_data *d)
{//{{{
    d->cls = (cls *)malloc(sizeof(cls));
    cls *_c = (cls *)d->cls;

    _c->pr = malloc(sizeof(struct precision));
    _c->ba = malloc(sizeof(struct background));
    _c->th = malloc(sizeof(struct thermo));
    _c->pt = malloc(sizeof(struct perturbs));
    _c->tr = malloc(sizeof(struct transfers));
    _c->pm = malloc(sizeof(struct primordial));
    _c->sp = malloc(sizeof(struct spectra));
    _c->nl = malloc(sizeof(struct nonlinear));
    _c->le = malloc(sizeof(struct lensing));
    _c->op = malloc(sizeof(struct output));
}//}}}

static
void run_class(all_data *d)
{//{{{
    cls *_c = (cls *)d->cls;

    printf("\trun_class\n");
    printf("\t\tbackground\n");
    SAFECLASS(background_init(_c->pr, _c->ba),
              _c->ba->error_message)
    printf("\t\tthermodynamics\n");
    SAFECLASS(thermodynamics_init(_c->pr, _c->ba, _c->th),
              _c->th->error_message)
    printf("\t\tperturbs\n");
    SAFECLASS(perturb_init(_c->pr, _c->ba, _c->th, _c->pt),
              _c->pt->error_message)
    printf("\t\tprimordial\n");
    SAFECLASS(primordial_init(_c->pr, _c->pt, _c->pm),
              _c->pm->error_message)
    printf("\t\tnonlinear\n");
    SAFECLASS(nonlinear_init(_c->pr, _c->ba, _c->th, _c->pt, _c->pm, _c->nl),
              _c->nl->error_message);
}//}}}

void init_class(all_data *d)
{//{{{
    printf("In class_interface.h -> init_class.\n");
    char **argv = (char **)malloc(3 * sizeof(char *));
    argv[1] = d->class_ini;
    argv[2] = d->class_pre;

    int argc = (strcmp(argv[2], "none")) ? 3 : 2;

    alloc_class(d);

    cls *_c = (cls *)d->cls;
    SAFECLASS(input_init_from_arguments(argc, argv, _c->pr, _c->ba, _c->th,
                                       _c->pt, _c->tr, _c->pm, _c->sp,
                                       _c->nl, _c->le, _c->op, _c->errmsg),
              _c->errmsg)

    run_class(d);
}//}}}

void free_class(all_data *d)
{//{{{
    if (d->cls != NULL)
    {
        cls *_c = (cls *)d->cls;

        free(_c->op);
        lensing_free(_c->le);
        free(_c->le);
        spectra_free(_c->sp);
        free(_c->sp);
        transfer_free(_c->tr);
        free(_c->tr);
        nonlinear_free(_c->nl);
        free(_c->nl);
        perturb_free(_c->pt);
        free(_c->pt);
        primordial_free(_c->pm);
        free(_c->pm);
        thermodynamics_free(_c->th);
        free(_c->th);
        background_free(_c->ba);
        free(_c->ba);
        free(_c->pr);

        free(d->cls);
    }
}//}}}

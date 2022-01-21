#ifndef CLASS_INTERFACE_H
#define CLASS_INTERFACE_H

#include "object.h"

#include "hmpdf.h"

typedef struct//{{{
{
    char *class_ini;
    char *class_pre;

    void /*struct precision*/ *pr;
    void /*struct background*/ *ba;
    void /*struct thermodynamics*/ *th;
    void /*struct primordial*/ *pm;
    void /*struct perturbations*/ *pt;
    void /*struct fourier*/ *nl;
    void /*struct transfer*/ *tr;
    void /*struct harmonic*/ *sp;
    void /*struct lensing*/ *le;
    void /*struct distortions*/ *sd;
    void /*struct output*/ *op;
}//}}}
class_interface_t;

int null_class_interface(hmpdf_obj *d);
int reset_class_interface(hmpdf_obj *d);
int init_class_interface(hmpdf_obj *d);

#endif

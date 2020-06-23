/*! [compile] */
/* gcc -I../include -o example example.c -L.. -lhmpdf */
/*! [compile] */
#include "hmpdf.h"

/*! [example_kappa_onepoint] */
int example_kappa_onepoint(void)
{
    /* construct some binedges */
    int Nbins = 100; double kappamin = 0.0; double kappamax = 0.3;
    double binedges[Nbins+1];
    for (int ii=0; ii<=Nbins; ii++)
        binedges[ii] = kappamin + (double)(ii)*(kappamax-kappamin)/(double)(Nbins);

    /* get a new hmpdf_obj */
    hmpdf_obj *d = hmpdf_new();
    if (!(d))
        return -1;

    /* initialize with default settings */
    if (hmpdf_init(d, "example.ini", hmpdf_kappa,
                   1.0/* source redshift */,
                   hmpdf_end_configs/* always include! */))
        return -1;

    /* get the one-point PDF */
    double op[Nbins];
    if (hmpdf_get_op(d, Nbins, binedges, op,
                     1/* include two-halo term */,
                     0/* don't include noise */))
        return -1;

    /* free memory associated with the hmpdf_obj */
    if (hmpdf_delete(d))
        return -1;

    /* do something with the one-point PDF ... */

    return 0;
}
/*! [example_kappa_onepoint] */

/*! [example_ell_filter] */
double example_ell_filter(double ell, void *p)
{
    double ell_max = *(double *)p;
    return (double)(ell < ell_max);
}
/*! [example_ell_filter] */

/*! [example_ell_filter_use] */
int example_ell_filter_use(void)
{
    hmpdf_obj *d = hmpdf_new();
    if (!(d))
        return -1;

    /* initialize with an ell-space filter which cuts off all
     * modes above ell = 5000
     */
    double ell_max = 5000.0;
    if (hmpdf_init(d, "example.ini", hmpdf_kappa, 1.0,
                   hmpdf_custom_ell_filter, &example_ell_filter,
                   hmpdf_custom_ell_filter_params, &ell_max,
                   hmpdf_end_configs))
        return -1;

    /* get your results ... */

    if (hmpdf_delete(d))
        return -1;

    return 0;
}
/*! [example_ell_filter_use] */

int main(void)
{
    return example_ell_filter_use();
}

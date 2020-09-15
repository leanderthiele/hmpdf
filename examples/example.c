/*! [compile] */
/* gcc -I../include -o example example.c -L.. -lhmpdf */
/*! [compile] */
#include <stdio.h>

#include "utils.h"
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
                   1.0/* source redshift */))
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
                   hmpdf_custom_ell_filter_params, &ell_max))
        return -1;

    /* get your results ... */

    if (hmpdf_delete(d))
        return -1;

    return 0;
}
/*! [example_ell_filter_use] */

double noisepwr(double ell, void *p)
{
    // this is a good testing value for tSZ
    //return 1e-17;
    return 0.0;
}

double kfilter(double k, double z, void *p)
{
    return 1.0;
}

int example_tsz_map(void)
{
    int Nbins = 50; double ymin=0.0; double ymax=1e-4;
    double binedges[Nbins+1];
    for (int ii=0; ii<=Nbins; ii++)
        binedges[ii] = ymin + (double)(ii)*(ymax-ymin)/(double)(Nbins);

    hmpdf_obj *d = hmpdf_new();
    if (!(d))
        return -1;

    if (hmpdf_init(d, "example.ini", hmpdf_tsz,
                   hmpdf_N_threads, 4,
                   hmpdf_verbosity, 5,
                   hmpdf_M_min, 1e11,
                   hmpdf_N_signal, 2048,

                   hmpdf_N_z, 40,
                   hmpdf_N_M, 40,
                   hmpdf_N_theta, 200,

                   hmpdf_custom_k_filter, &kfilter,
                   
                   hmpdf_map_pixelgrid, 2,
                   hmpdf_map_poisson, 0,
                   
                   hmpdf_pixel_side, 1.0,
                   hmpdf_map_fsky, 1e-2))
        return -1;

    double op[Nbins];
    if (hmpdf_get_op(d, Nbins, binedges, op, 0, 0))
        return -1;

    savetxt("test_op", Nbins, 1, op);

    double map_op[Nbins];
    if (hmpdf_get_map_op(d, Nbins, binedges, map_op, 0))
        return -1;

    double *map;
    long Nside;
    if (hmpdf_get_map(d, &map, &Nside, 0))
        return -1;

    savetxt("test_map_op", Nbins, 1, map_op);
    tofile("test_map", Nside*Nside, 1, map);

    free(map);

    if (hmpdf_delete(d))
        return -1;

    return 0;
}

int main(void)
{
    if (example_tsz_map())
    {
        fprintf(stderr, "failed\n");
        return -1;
    }
    else
    {
        return 0;
    }
}

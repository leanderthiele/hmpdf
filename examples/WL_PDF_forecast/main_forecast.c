/* Driver code to reproduce Fisher matrix computation
 * for the WL convergence PDF, using analytic covariance matrix.
 *
 * This code computes PDFs and covariance matrix for the fiducial model
 * and varied cosmological and concentration model parameters.
 *
 * The executable takes one command line argument, integer IDX :
 *      IDX = 0     fid       -- fiducial model
 *            1..6  varcosmo  -- varied cosmological parameters
 *            7..12 varconc   -- varied concentration model parameters
 *            
 *            (we have 3 parameters on the cosmology and concentration side
 *             respectively, each with variations plus and minus)
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "hmpdf.h"
#include "utils.h"

#define KMIN -0.03 // minimum kappa 
#define KMAX 0.3   // maximum kappa
#define SN   0.3   // shape noise parameter
#define NGAL 45.0  // galaxy density per arcmin^2

// prefix to read .ini files from
#define INI_PATH "./cosmologies/"
// prefix to put output
#define OUT_PATH "./results/productionrun_"
// Wiener filter file
#define WIENER_FILE "./Wiener_filter.txt"
// CLASS precision file
#define PRE_FILE "./cl_ref.pre"

// we bin into different resolutions to check if it influences constraints
const int Nbinnings = 3;
const int binnings[] = {100, 200, 400};

// convenience macro for common settings
#define COMMON_INIT(inifile, ...)                 \
    hmpdf_init(d, inifile, hmpdf_kappa, 1.0,      \
               hmpdf_verbosity, 5,                \
               hmpdf_class_pre, PRE_FILE,         \
               hmpdf_N_threads, 4,                \
               hmpdf_pixel_side, 0.41,            \
               hmpdf_signal_max, 0.4,             \
               hmpdf_N_theta, 500,                \
               hmpdf_M_max, 1e17,                 \
               hmpdf_N_M, 80,                     \
               hmpdf_N_z, 80,                     \
               hmpdf_N_signal, 1024L,             \
               hmpdf_N_phi, 2000,                 \
               hmpdf_rout_scale, 1.6,             \
               hmpdf_custom_ell_filter, &Well,    \
               hmpdf_custom_ell_filter_params, p, \
               hmpdf_noise_pwr, &Nell,            \
               ##__VA_ARGS__)
// note : dealing with empty __VA_ARGS__ this way is an extension
//        not supported by every compiler (but gcc does)

// Wiener filter interpolator
double Well(double ell, void *p)
{//{{{
    interp1d *i = (interp1d *)p;
    double out;
    interp1d_eval(i, ell, &out);
    return out;
}//}}}

#define SQUARE(x) ((x)*(x))
// noise power spectrum
double Nell(double ell, void *p)
{//{{{
    return SQUARE(SN) / NGAL
           * SQUARE(M_PI/180.0/60.0); // convert to rad units
}//}}}
#undef SQUARE

int _run(hmpdf_obj *d, char *outname, int Nbins)
{//{{{
    double binedges[Nbins+1];
    linspace(Nbins+1, KMIN, KMAX, binedges);

    double *PDF      = malloc(Nbins * sizeof(double));
    double *PDFnoisy = malloc(Nbins * sizeof(double));
    double *COV      = malloc(Nbins * Nbins * sizeof(double));
    double *COVnoisy = malloc(Nbins * Nbins * sizeof(double));

    if (hmpdf_get_op(d, Nbins, binedges, PDF, 1, 0))
        return -1;
    if (hmpdf_get_op(d, Nbins, binedges, PDFnoisy, 1, 1))
        return -1;
    if (hmpdf_get_cov(d, Nbins, binedges, COV, 0))
        return -1;
    if (hmpdf_get_cov(d, Nbins, binedges, COVnoisy, 1))
        return -1;

    char pdfname[512];
    char covname[512];
    char pdfname_noisy[512];
    char covname_noisy[512];
    sprintf(pdfname, "%s_%d_pdf", outname, Nbins);
    sprintf(pdfname_noisy, "%s_%d_pdf_noisy", outname, Nbins);
    sprintf(covname, "%s_%d_cov", outname, Nbins);
    sprintf(covname_noisy, "%s_%d_cov_noisy", outname, Nbins);

    savetxt(pdfname, Nbins, 1, PDF);
    savetxt(pdfname_noisy, Nbins, 1, PDFnoisy);
    savetxt(covname, Nbins*Nbins, 1, COV);
    savetxt(covname_noisy, Nbins*Nbins, 1, COVnoisy);

    free(PDF); free(PDFnoisy); free(COV); free(COVnoisy);

    return 0;
}//}}}

int run(hmpdf_obj *d, char *outname)
{//{{{
    for (int ii=0; ii<Nbinnings; ii++)
    {
        if (_run(d, outname, binnings[ii]))
            return -1;
    }

    return 0;
}//}}}

int fid(hmpdf_obj *d, void *p)
{//{{{
    printf("-------------fid---------------\n");
    if (COMMON_INIT(INI_PATH"fid.ini"))
        return -1;
    if (run(d, OUT_PATH"fid"))
        return -1;

    return 0;
}//}}}

int varcosmo(hmpdf_obj *d, void *p, int idx /*0..5*/)
{//{{{
    printf("-----------varcosmo-------------\n");
    char names[][32] = {"As", "Om", "Mn"};
    char pm[] = {'m', 'p'};
    char inname[512];
    char outname[512];

    // ii = 0..2 (index of parameter to vary)
    int ii = idx/2;
    // jj = 0..1 (index of plus or minus variation)
    int jj = idx%2;

    printf("\t---%s%c----\n", names[ii], pm[jj]);
    sprintf(inname, INI_PATH"%s%c.ini", names[ii], pm[jj]);
    sprintf(outname, OUT_PATH"%s%c", names[ii], pm[jj]);
    if (COMMON_INIT(inname))
        return -1;

    if (run(d, outname))
        return -1;

    return 0;
}//}}}

void conc_at_peak(double *concp, double *out)
{//{{{
    double mref = pow(10.0, 14.5)/1e12; // where kernel peaks approximately
    double zref = 0.35;
    for (int ii=0; ii<3; ii++)
    {
        out[ii] = concp[3*ii+0] * pow(mref, concp[3*ii+1])
                  * pow(1.0+zref, concp[3*ii+2]);
    }
}//}}}

int varconc(hmpdf_obj *d, void *p, int idx/*0..5*/)
{//{{{
    printf("-----------varconc-------------\n");
                      //    A       B      C
    double def_conc[] = { 5.71, -0.087, -0.47,   // M200c
                          7.85, -0.081, -0.71,   // Mvir
                         10.14, -0.081, -1.01, };// M200m
    double ref_conc[3];
    conc_at_peak(def_conc, ref_conc);
    double new_conc[3];

    char names[] = {'A', 'B', 'C'};
    char pm[] = {'p', 'm'};
    char outname[256];
    double this_conc[9];

    int ii = idx/2;
    int jj = idx%2;

    printf("\t-----%c%c------\n", names[ii], pm[jj]);
    memcpy(this_conc, def_conc, 9*sizeof(double));
    for (int kk=0; kk<3; kk++)
    {
        this_conc[ii+3*kk] *= (jj==0) ? 1.1 : 0.9;
    }
    if (ii>0) // not changing the normalization
    {
        conc_at_peak(this_conc, new_conc);
        for (int kk=0; kk<3; kk++)
        {
            this_conc[3*kk+0] *= ref_conc[kk]/new_conc[kk];
        }
    }

    sprintf(outname, OUT_PATH"%c%c", names[ii], pm[jj]);
    if (COMMON_INIT(INI_PATH"fid.ini",
                    hmpdf_Duffy08_conc_params, this_conc))
        return -1;
    if (run(d, outname))
        return -1;

    return 0;
}//}}}

int main(int argc, char **argv)
{//{{{
    if (argc != 2)
        return -1;
    int IDX = atoi(argv[1]);

    // construct Wiener filter interpolator
    int Nlines;
    double **x = loadtxt(WIENER_FILE, &Nlines, 2);
    interp1d *Wiener_interp;
    new_interp1d(Nlines, x[0], x[1],
                 x[1][0]/* fill value when below lower bound */,
                 0.0/* fill value when above upper bound */,
                 interp_linear, NULL, &Wiener_interp);

    hmpdf_obj *d = hmpdf_new();
    if (!(d))
        return -1;

    if (IDX == 0)
    {
        if (fid(d, (void *)Wiener_interp))
            return -1;
    }
    else if (IDX <= 6)
    {
        if (varcosmo(d, (void *)Wiener_interp, IDX-1))
            return -1;
    }
    else if (IDX <= 12)
    {
        if (varconc(d, (void *)Wiener_interp, IDX-7))
            return -1;
    }
    else
    {
        return -1;
    }

    if (hmpdf_delete(d))
        return -1;

    return 0;
}//}}}

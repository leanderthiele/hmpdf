#ifndef UTILS_H
#define UTILS_H

#define SPEEDOFLIGHT 2997.92458  // 100 km/s
#define GNEWTON      4.30091e-13 // Mpc/Msun (100 km/s)^2
#define SIGMATHOMSON 6.98684e-74 // Mpc^2
#define MELECTRON    4.58110e-61 // Msun

#define HPLANCK      6.62607004e-34 // SI
#define KBOLTZMANN   1.38064852e-23 // SI

#include <complex.h>

#include <gsl/gsl_interp.h>

#include "hmpdf.h"

typedef enum//{{{
{
    interp_linear,
    interp_polynomial,
    interp_cspline,
    interp_cspline_periodic,
    interp_akima,
    interp_akima_periodic,
    interp_steffen,
}//}}}
interp_mode;

typedef struct gnuplot_s gnuplot;

typedef struct interp1d_s interp1d;

int ispwr2(int N, int *k);

void linspace(int N, double xmin, double xmax, double *x);
void logspace(int N, double xmin, double xmax, double *x);
void zero_real(int N, double *x);
void zero_comp(int N, complex *x);
void reverse(int N, double *in, double *out);

int wait(void);
gnuplot *plot(gnuplot *gp, int N, double *x, double *y);
gnuplot *plot_comp(gnuplot *gp, int N, double *x, complex *y, int mode);
void show(gnuplot *gp);

void savetxt(char *fname, int Nlines, int Nvec, ...);
double **loadtxt(char *fname, int *Nlines, int Nvec);
void tofile(char *fname, int Nlines, int Nvec, ...);
double **fromfile(char *fname, int *Nlines, int Nvec);
int isfile(char *fname);

int num_cores(void);
int this_core(void);

interp1d *new_interp1d(int N, double *x, double *y, double ylo, double yhi,
                       interp_mode m, gsl_interp_accel *a);
void delete_interp1d(interp1d *interp);
double interp1d_eval(interp1d *interp, double x);
double interp1d_eval_deriv(interp1d *interp, double x);
double interp1d_eval_integ(interp1d *interp, double a, double b);

#endif

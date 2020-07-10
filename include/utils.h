#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <complex.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

#include "hmpdf.h"

#define SPEEDOFLIGHT 2997.92458  // 100 km/s
#define GNEWTON      4.30091e-13 // Mpc/Msun (100 km/s)^2
#define SIGMATHOMSON 6.98684e-74 // Mpc^2
#define MELECTRON    4.58110e-61 // Msun

// branch prediction macros for some obvious cases
//     (check if those are available)
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

//ERRLOC{{{
#define ERRLOC                                   \
    fprintf(stderr, "Error in %s line %d "       \
                    "(see upward for cause) \n", \
                    __FILE__, __LINE__);         \
//}}}

// when calling external code (CLASS, GSL) we set errno=0 afterwards
//      in order not to catch internal errors the external code doesn't
//      consider important

//SAFECLASS{{{
#define SAFECLASS(expr, errmsg)                 \
    if (UNLIKELY(expr==_FAILURE_))              \
    {                                           \
        fprintf(stderr, "***CLASS error: %s\n", \
                        errmsg);                \
        fflush(stderr);                         \
        ERRLOC                                  \
        hmpdf_status = 1;                       \
        return hmpdf_status;                    \
    }                                           \
    errno = 0;                                  \
//}}}

int hmpdf_status_update(int *status, int result);

void new_gsl_error_handler(const char *reason, const char *file,
                           int line, int gsl_errno);

//SAFEGSL{{{
#define SAFEGSL(expr)                                \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status,  \
                                     expr)))         \
    {                                                \
        fprintf(stderr, "***GSL error: %s\n",        \
                        gsl_strerror(hmpdf_status)); \
        fflush(stderr);                              \
        ERRLOC                                       \
        return hmpdf_status;                         \
    }                                                \
    errno = 0;                                       \
//}}}

//SAFEGSL_NORETURN{{{
#define SAFEGSL_NORETURN(expr)                       \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status,  \
                            expr)))                  \
    {                                                \
        fprintf(stderr, "***GSL error: %s\n",        \
                        gsl_strerror(hmpdf_status)); \
        fflush(stderr);                              \
        ERRLOC                                       \
    }                                                \
    errno = 0;                                       \
//}}}

//SAFEHMPDF{{{
#define SAFEHMPDF(expr)                             \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status, \
                                     expr)))        \
    {                                               \
        ERRLOC                                      \
        return 1;                                   \
    }                                               \
//}}}

//SAFEHMPDF_NORETURN{{{
#define SAFEHMPDF_NORETURN(expr)                    \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status, \
                                     expr)))        \
    {                                               \
        ERRLOC                                      \
    }                                               \
//}}}

//SAFEALLOC{{{
#define SAFEALLOC(dt,var,expr)             \
    dt var = expr;                         \
    if (UNLIKELY(!(var)))                  \
    {                                      \
        fprintf(stderr, "***OOM error\n"); \
        fflush(stderr);                    \
    }                                      \
    SAFEHMPDF(!(var))                      \
//}}}

//SAFEALLOC_NORETURN{{{
#define SAFEALLOC_NORETURN(dt,var,expr)    \
    dt var = expr;                         \
    if (UNLIKELY(!(var)))                  \
    {                                      \
        fprintf(stderr, "***OOM error\n"); \
        fflush(stderr);                    \
    }                                      \
    SAFEHMPDF_NORETURN(!(var))             \
//}}}

//CHECKERR{{{
#define CHECKERR                                    \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status, \
                                     errno)))       \
    {                                               \
        perror("***C error (likely wrong "          \
               "description)");                     \
        fprintf(stderr, "%s\n", strerror(errno));   \
        fflush(stderr);                             \
        ERRLOC                                      \
    }                                               \
//}}}

//STARTFCT{{{
#define STARTFCT          \
    int hmpdf_status = 0; \
    errno = 0;            \
//}}}

//ENDFCT{{{
#define ENDFCT           \
    CHECKERR             \
    return hmpdf_status; \
//}}}

//HMPDFERR_NORETURN{{{
#define HMPDFERR_NORETURN(...)           \
    fprintf(stderr, "***hmpdf error: "); \
    fprintf(stderr, __VA_ARGS__);        \
    fprintf(stderr, "\n");               \
    fflush(stderr);                      \
    ERRLOC                               \
    hmpdf_status = 1;                    \
//}}}

//HMPDFERR{{{
#define HMPDFERR(...)              \
    HMPDFERR_NORETURN(__VA_ARGS__) \
    return hmpdf_status;           \
//}}}

//HMPDFPRINT{{{
#define HMPDFPRINT(level, ...)        \
    if (level <= d->verbosity)        \
    {                                 \
        fprintf(stdout, __VA_ARGS__); \
        fflush(stdout);               \
    }                                 \
//}}}

int ispwr2(int N, int *k);

void linspace(int N, double xmin, double xmax, double *x);
void logspace(int N, double xmin, double xmax, double *x);
void zero_real(int N, double *x);
void zero_comp(int N, complex *x);
void reverse(int N, double *in, double *out);
int roll1d(int N, int N1, int stride, double *in, double *out);
int roll2d(int N, int N1, double *in, double *out);
int not_monotonic(int N, double *x, int *problems);
int all_zero(int N, double *x, double threshold);

int wait(void);

typedef struct gnuplot_s gnuplot;
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
typedef struct interp1d_s interp1d;
int new_interp1d(int N, double *x, double *y,
                 double ylo, double yhi,
                 interp_mode m, gsl_interp_accel *a, interp1d **out);
void delete_interp1d(interp1d *interp);
int interp1d_eval(interp1d *interp, double x, double *out);
int interp1d_eval1(interp1d *interp, double x, int *inrange, double *out);
int interp1d_eval_deriv(interp1d *interp, double x, double *out);
int interp1d_eval_deriv1(interp1d *interp, double x, int *inrange, double *out);
int interp1d_eval_integ(interp1d *interp, double a, double b, double *out);

typedef enum//{{{
{
    interp2d_bilinear,
    interp2d_bicubic,
}//}}}
interp2d_mode;
typedef struct interp2d_s interp2d;
int new_interp2d(int N, double *x, double *z,
                 double zlo, double zhi,
                 interp2d_mode m, gsl_interp_accel *a, interp2d **out);
void delete_interp2d(interp2d *interp);
int interp2d_eval(interp2d *interp, double x, double y, double *out);

int bin_1d(int N, double *x, double *y,
           int Nbins, double *binedges, double *out, interp_mode m);

int bin_2d(int N, double *x, double *z, int Nsample,
           int Nbins, double *binedges, double *out, interp2d_mode m);

#endif

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

static inline int
hmpdf_status_update(int *status, int result)
{//{{{
    *status = result;
    return result;
}//}}}

void
new_gsl_error_handler(const char *reason, const char *file,
                      int line, int gsl_errno);

// branch prediction macros for some obvious cases
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

//ERRLOC{{{
#define ERRLOC                                     \
    do {                                           \
    fprintf(stderr, "Error in %s (%s; line %d) "   \
                    "(see upward for cause) \n",   \
                    __FILE__, __func__, __LINE__); \
    } while (0)                                    \
//}}}

// when calling external code (CLASS, GSL) we set errno=0 afterwards
//      in order not to catch internal errors the external code doesn't
//      consider important

//SAFECLASS{{{
#define SAFECLASS(expr, errmsg)                 \
    do {                                        \
    if (UNLIKELY(expr==_FAILURE_))              \
    {                                           \
        fprintf(stderr, "***CLASS error: %s\n", \
                        errmsg);                \
        fflush(stderr);                         \
        ERRLOC;                                 \
        hmpdf_status = 1;                       \
        return hmpdf_status;                    \
    }                                           \
    errno = 0;                                  \
    } while (0)                                 \
//}}}

//SAFEGSL{{{
#define SAFEGSL(expr)                                \
    do {                                             \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status,  \
                                     expr)))         \
    {                                                \
        fprintf(stderr, "***GSL error: %s\n",        \
                        gsl_strerror(hmpdf_status)); \
        fflush(stderr);                              \
        ERRLOC;                                      \
        return hmpdf_status;                         \
    }                                                \
    errno = 0;                                       \
    } while (0)                                      \
//}}}

//SAFEGSL_NORETURN{{{
#define SAFEGSL_NORETURN(expr)                       \
    do {                                             \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status,  \
                            expr)))                  \
    {                                                \
        fprintf(stderr, "***GSL error: %s\n",        \
                        gsl_strerror(hmpdf_status)); \
        fflush(stderr);                              \
        ERRLOC;                                      \
    }                                                \
    errno = 0;                                       \
    } while (0)                                      \
//}}}

//SAFEHMPDF{{{
#define SAFEHMPDF(expr)                             \
    do {                                            \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status, \
                                     expr)))        \
    {                                               \
        ERRLOC;                                     \
        return 1;                                   \
    }                                               \
    } while (0)                                     \
//}}}

//SAFEHMPDF_NORETURN{{{
#define SAFEHMPDF_NORETURN(expr)                    \
    do {                                            \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status, \
                                     expr)))        \
    {                                               \
        ERRLOC;                                     \
    }                                               \
    } while (0)                                     \
//}}}

//SAFEALLOC{{{
#define SAFEALLOC(var,expr)                \
    do {                                   \
    var = expr;                            \
    if (UNLIKELY(!(var)))                  \
    {                                      \
        fprintf(stderr, "***OOM error\n"); \
        fflush(stderr);                    \
    }                                      \
    SAFEHMPDF(!(var));                     \
    } while (0)                            \
//}}}

//SAFEALLOC_NORETURN{{{
#define SAFEALLOC_NORETURN(var,expr)       \
    do {                                   \
    var = expr;                            \
    if (UNLIKELY(!(var)))                  \
    {                                      \
        fprintf(stderr, "***OOM error\n"); \
        fflush(stderr);                    \
    }                                      \
    SAFEHMPDF_NORETURN(!(var));            \
    } while (0)                            \
//}}}

//CHECKERR{{{
#define CHECKERR                                    \
    do {                                            \
    if (UNLIKELY(hmpdf_status_update(&hmpdf_status, \
                                     errno)))       \
    {                                               \
        perror("***C error (likely wrong "          \
               "description)");                     \
        fprintf(stderr, "%s\n", strerror(errno));   \
        fflush(stderr);                             \
        ERRLOC;                                     \
    }                                               \
    } while (0)                                     \
//}}}

//STARTFCT{{{
#define STARTFCT          \
    int hmpdf_status = 0; \
    errno = 0;            \
//}}}

//ENDFCT{{{
#define ENDFCT           \
    CHECKERR;            \
    return hmpdf_status; \
//}}}

//HMPDFERR_NORETURN{{{
#define HMPDFERR_NORETURN(...)           \
    do {                                 \
    fprintf(stderr, "***hmpdf error: "); \
    fprintf(stderr, __VA_ARGS__);        \
    fprintf(stderr, "\n");               \
    fflush(stderr);                      \
    ERRLOC;                              \
    hmpdf_status = 1;                    \
    } while (0)                          \
//}}}

//HMPDFERR{{{
#define HMPDFERR(...)               \
    do {                            \
    HMPDFERR_NORETURN(__VA_ARGS__); \
    return hmpdf_status;            \
    } while (0)                     \
//}}}

//HMPDFCHECK_NORETURN{{{
#define HMPDFCHECK_NORETURN(expr, ...)  \
    do {                                \
    if (UNLIKELY(expr))                 \
    {                                   \
        HMPDFERR_NORETURN(__VA_ARGS__); \
    }                                   \
    } while (0)                         \
//}}}

//HMPDFCHECK{{{
#define HMPDFCHECK(expr, ...)  \
    do {                       \
    if (UNLIKELY(expr))        \
    {                          \
        HMPDFERR(__VA_ARGS__); \
    }                          \
    } while (0)                \
//}}}

//HMPDFPRINT{{{
#define HMPDFPRINT(level, ...)        \
    do {                              \
    if (level <= d->verbosity)        \
    {                                 \
        fprintf(stdout, __VA_ARGS__); \
        fflush(stdout);               \
    }                                 \
    } while (0)                       \
//}}}

//HMPDFWARN{{{
#define HMPDFWARN(...)                                      \
    do {                                                    \
    if (d->warn_is_err > 0)                                 \
    {                                                       \
        HMPDFERR_NORETURN(__VA_ARGS__);                     \
        fprintf(stderr, "\t\t>This error can be suppressed "\
                        "by changing the option "           \
                        "hmpdf_warn_is_err, "               \
                        "in case you know what you are "    \
                        "doing.\n");                        \
        fflush(stderr);                                     \
        return hmpdf_status;                                \
    }                                                       \
    else if (d->warn_is_err == 0)                           \
    {                                                       \
        fprintf(stderr, "***hmpdf warning: ");              \
        fprintf(stderr, __VA_ARGS__);                       \
        fprintf(stderr, "\t\t>This warning can be muted "   \
                        "or turned into an error by "       \
                        "changing the option "              \
                        "hmpdf_warn_is_err.\n");            \
        fflush(stderr);                                     \
    }                                                       \
    } while (0)                                             \
//}}}

//CONTINUE_IF_ERR{{{
#define CONTINUE_IF_ERR         \
    if (UNLIKELY(hmpdf_status)) \
    {                           \
        continue;               \
    }                           \
//}}}

int ispwr2(int N, int *k);

void linspace(int N, double xmin, double xmax, double *x);
void logspace(int N, double xmin, double xmax, double *x);
void zero_real(int N, double *x);
void zero_comp(int N, double complex *x);
void reverse(int N, double *in, double *out);
int not_monotonic(int N, double *x, int sgn);
int all_zero(int N, double *x, double threshold);

int wait(void);

#ifdef GNUPLOT
typedef struct gnuplot_s gnuplot;
gnuplot *plot(gnuplot *gp, int N, double *x, double *y);
gnuplot *plot_comp(gnuplot *gp, int N, double *x, double complex *y, int mode);
void show(gnuplot *gp);
#endif

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
const gsl_interp_type *interp1d_type(interp_mode m);
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

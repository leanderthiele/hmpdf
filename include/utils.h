#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <complex.h>

#ifdef _OPENMP
#   include <omp.h>
#endif

#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

#include "hmpdf.h"

#define SPEEDOFLIGHT 2997.92458  // 100 km/s
#define GNEWTON      4.30091e-13 // Mpc/Msun (100 km/s)^2
#define SIGMATHOMSON 6.98684e-74 // Mpc^2
#define MELECTRON    4.58110e-61 // Msun
#define RADPERARCMIN (M_PI/180.0/60.0) // rad/arcmin

static inline int
hmpdf_status_update(int *status, int result)
{//{{{
    *status |= result;
    return result;
}//}}}

void
new_gsl_error_handler(const char *reason, const char *file,
                      int line, int gsl_errno);

// branch prediction macros for some obvious cases
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

//ERRLOC{{{
#ifdef DEBUG
    #define ERRLOC                                         \
        do {                                               \
            fprintf(stderr, "Error in %s (%s; line %d) "   \
                            "(see upward for cause) \n",   \
                            __FILE__, __func__, __LINE__); \
        } while (0)
#endif
//}}}

// when calling external code (CLASS, GSL) we set errno=0 afterwards
//      in order not to catch internal errors the external code doesn't
//      consider important

//SAFECLASS{{{
#ifdef DEBUG
#   define SAFECLASS(expr, errmsg)                      \
        do {                                            \
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
        } while (0)
#else
#   define SAFECLASS(expr, errmsg)  \
        do {                        \
            expr;                   \
            errno = 0;              \
        } while (0)
#endif
//}}}

//SAFEGSL{{{
#ifdef DEBUG
#   define SAFEGSL(expr)                                     \
        do {                                                 \
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
        } while (0)
#else
#   define SAFEGSL(expr)  \
        do {              \
            expr;         \
            errno = 0;    \
        } while (0) 
#endif
//}}}

//SAFEGSL_NORETURN{{{
#ifdef DEBUG
#   define SAFEGSL_NORETURN(expr)                            \
        do {                                                 \
            if (UNLIKELY(hmpdf_status_update(&hmpdf_status,  \
                                    expr)))                  \
            {                                                \
                fprintf(stderr, "***GSL error: %s\n",        \
                                gsl_strerror(hmpdf_status)); \
                fflush(stderr);                              \
                ERRLOC;                                      \
            }                                                \
            errno = 0;                                       \
        } while (0)
#else
#   define SAFEGSL_NORETURN(expr)  \
        do {                       \
            expr;                  \
            errno = 0;             \
        } while (0)
#endif
//}}}

//SAFEHMPDF{{{
#ifdef DEBUG
#   define SAFEHMPDF(expr)                                  \
        do {                                                \
            if (UNLIKELY(hmpdf_status_update(&hmpdf_status, \
                                             expr)))        \
            {                                               \
                ERRLOC;                                     \
                return 1;                                   \
            }                                               \
        } while (0)
#else
#   define SAFEHMPDF(expr)  \
        do {                \
            expr;           \
        } while (0)
#endif
//}}}

//SAFEHMPDF_NORETURN{{{
#ifdef DEBUG
#   define SAFEHMPDF_NORETURN(expr)                         \
        do {                                                \
            if (UNLIKELY(hmpdf_status_update(&hmpdf_status, \
                                             expr)))        \
            {                                               \
                ERRLOC;                                     \
            }                                               \
        } while (0)
#else
#   define SAFEHMPDF_NORETURN(expr)  \
        do {                         \
            expr;                    \
        } while (0)
#endif
//}}}

#ifdef __GNUC__
#   define WCONVERSIONOFF _Pragma("GCC diagnostic ignored \"-Wsign-conversion\"")
#   define WCONVERSIONON  _Pragma("GCC diagnostic pop")
#else
#   define WCONVERSIONOFF 
#   define WCONVERSIONON
#endif

//SAFEALLOC{{{
#ifdef DEBUG
#   define SAFEALLOC(var,expr)                       \
        do {                                         \
            WCONVERSIONOFF                           \
            var = expr;                              \
            WCONVERSIONON                            \
            HMPDFCHECK(!(var),                       \
                       "memory allocation failed."); \
        } while (0)
#else
#   define SAFEALLOC(var,expr)  \
        do {                    \
            WCONVERSIONOFF      \
            var = expr;         \
            WCONVERSIONON       \
        } while (0)
#endif
//}}}

//SAFEALLOC_NORETURN{{{
#ifdef DEBUG
#   define SAFEALLOC_NORETURN(var,expr)              \
        do {                                         \
            WCONVERSIONOFF                           \
            var = expr;                              \
            WCONVERSIONON                            \
            HMPDFCHECK_NORETURN(!(var),              \
                       "memory allocation failed."); \
        } while (0)
#else
#   define SAFEALLOC_NORETURN(var,expr)  \
        do {                             \
            WCONVERSIONOFF               \
            var = expr;                  \
            WCONVERSIONON                \
        } while (0)
#endif
//}}}

//CHECKERR{{{
#ifdef DEBUG
#   define CHECKERR                                         \
        do {                                                \
            if (UNLIKELY(hmpdf_status_update(&hmpdf_status, \
                                             errno)))       \
            {                                               \
                perror("***C error (likely wrong "          \
                       "description)");                     \
                fprintf(stderr, "%s\n", strerror(errno));   \
                fflush(stderr);                             \
                ERRLOC;                                     \
            }                                               \
        } while (0)
#endif
//}}}

//CHECKINIT{{{
#ifdef DEBUG
#   define CHECKINIT                                         \
        do {                                                 \
            HMPDFCHECK(!(d->inited),                         \
                       "hmpdf_init has not been called yet " \
                       "or previous call failed.");          \
        } while (0)
#else
#   define CHECKINIT \
        do { } while (0)
#endif
//}}}

//STARTFCT -- no semicolon!{{{
#ifdef DEBUG
#   define STARTFCT           \
        int hmpdf_status = 0; \
        errno = 0;
#else
#   define STARTFCT
#endif
//}}}

//ENDFCT -- no semicolon!{{{
#ifdef DEBUG
#   define ENDFCT            \
        CHECKERR;            \
        return hmpdf_status;
#else
#   define ENDFCT \
        return 0;
#endif
//}}}

//HMPDFERR_NORETURN{{{
#ifdef DEBUG
#   define HMPDFERR_NORETURN(...)                \
        do {                                     \
            fprintf(stderr, "***hmpdf error: "); \
            fprintf(stderr, __VA_ARGS__);        \
            fprintf(stderr, "\n");               \
            fflush(stderr);                      \
            ERRLOC;                              \
            hmpdf_status = 1;                    \
        } while (0)
#else
#   define HMPDFERR_NORETURN(...) \
        do { } while (0)
#endif
//}}}

//HMPDFERR{{{
#ifdef DEBUG
#   define HMPDFERR(...)                    \
        do {                                \
            HMPDFERR_NORETURN(__VA_ARGS__); \
            return hmpdf_status;            \
        } while (0)
#else
#   define HMPDFERR(...) \
        do { } while (0)
#endif
//}}}

//HMPDFCHECK_NORETURN{{{
#ifdef DEBUG
#   define HMPDFCHECK_NORETURN(expr, ...)       \
        do {                                    \
            if (UNLIKELY(expr))                 \
            {                                   \
                HMPDFERR_NORETURN(__VA_ARGS__); \
            }                                   \
        } while (0)
#else
#   define HMPDFCHECK_NORETURN(expr, ...) \
        do { } while (0)
#endif
//}}}

//HMPDFCHECK{{{
#ifdef DEBUG
#   define HMPDFCHECK(expr, ...)       \
        do {                           \
            if (UNLIKELY(expr))        \
            {                          \
                HMPDFERR(__VA_ARGS__); \
            }                          \
        } while (0)
#else
#define HMPDFCHECK(expr, ...) \
        do { } while (0)
#endif
//}}}

//HMPDFPRINT{{{
#define HMPDFPRINT(level, ...)            \
    do {                                  \
        if (level <= d->verbosity)        \
        {                                 \
            fprintf(stdout, __VA_ARGS__); \
            fflush(stdout);               \
        }                                 \
    } while (0)
//}}}

//HMPDFWARN{{{
#ifdef DEBUG
#   define HMPDFWARN(...)                                            \
        do {                                                         \
            if (d->warn_is_err > 0)                                  \
            {                                                        \
                HMPDFERR_NORETURN(__VA_ARGS__);                      \
                fprintf(stderr, "\t\t>This error can be suppressed " \
                                "by changing the option "            \
                                "hmpdf_warn_is_err, "                \
                                "in case you know what you are "     \
                                "doing.\n");                         \
                fflush(stderr);                                      \
                return hmpdf_status;                                 \
            }                                                        \
            else if (d->warn_is_err == 0)                            \
            {                                                        \
                fprintf(stderr, "***hmpdf warning: ");               \
                fprintf(stderr, __VA_ARGS__);                        \
                fprintf(stderr, "\t\t>This warning can be muted "    \
                                "or turned into an error by "        \
                                "changing the option "               \
                                "hmpdf_warn_is_err.\n");             \
                fflush(stderr);                                      \
            }                                                        \
        } while (0)
#else
#   define HMPDFWARN(...) \
        do { } while (0)
#endif
//}}}

//CONTINUE_IF_ERR -- no semicolon!{{{
#ifdef DEBUG
#   define CONTINUE_IF_ERR          \
        if (UNLIKELY(hmpdf_status)) \
        {                           \
            continue;               \
        }
#else
#   define CONTINUE_IF_ERR
#endif
//}}}

// UNUSED{{{
#ifdef __GNUC__
#   define UNUSED(x) x##_UNUSED __attribute__((unused))
#else
#   define UNUSED(x) x
#endif
//}}}

//SETARRNULL{{{
#define SETARRNULL(arr, len)    \
    do {                        \
        for (int arridx=0;      \
             arridx < len;      \
             arridx++)          \
        {                       \
            arr[arridx] = NULL; \
        }                       \
    } while(0)
//}}}

// produces int hrs, min, sec
#define TIMESPLIT(t)                             \
    do {                                         \
        hrs = (int)floor(t/60.0/60.0);           \
        min = (int)floor(t/60.0                  \
                         - 60.0*(double)hrs);    \
        sec = (int)round(t                       \
                         - 60.0*60.0*(double)hrs \
                         - 60.0*(double)min);    \
    } while (0)

// assumes hmpdf_obj *d and time_t start_time are in scope
#define TIMEREMAIN(done, tot, fctname)                        \
    do {                                                      \
        time_t t1 = time(NULL);                               \
        double delta_time = difftime(t1, start_time);         \
        double remains = delta_time/(double)done              \
                         *(double)(tot - done);               \
        int done_perc = (int)round(100.0*(double)(done)       \
                                   / (double)(tot));          \
        int hrs, min, sec;                                    \
        TIMESPLIT(remains);                                   \
        HMPDFPRINT(1, "\t\t%3d %% done, "                     \
                      "%.2d hrs %.2d min %.2d sec remaining " \
                      "in "fctname".\n",                      \
                      done_perc, hrs, min, sec);              \
    } while (0)
#define TIMEELAPSED(fctname)                                  \
    do {                                                      \
        time_t t1 = time(NULL);                               \
        double elapsed = difftime(t1, start_time);            \
        int hrs, min, sec;                                    \
        TIMESPLIT(elapsed);                                   \
        HMPDFPRINT(1, "\t\tspent %.2d hrs %.2d min %.2d sec " \
                      "in "fctname".\n",                      \
                      hrs, min, sec);                         \
    } while (0)

int ispwr2(int N, int *k);

int linspace(int N, double xmin, double xmax, double *x);
int logspace(int N, double xmin, double xmax, double *x);
void zero_real(long N, double *x);
void zero_comp(long N, double complex *x);
void reverse(int N, double *in, double *out);
int not_monotonic(int N, double *x, int sgn);
int all_zero(int N, double *x, double threshold);

#define WAVENR(N, grid, idx) \
    (idx <= N/2) ? grid[idx] : -grid[N-idx]

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

#ifdef _OPENMP
#   define THIS_THREAD omp_get_thread_num()
#else
#   define THIS_THREAD 0
#endif

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

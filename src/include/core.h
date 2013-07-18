#ifndef __CORE_H__
#define __CORE_H__

#ifdef __cplusplus
  #include <cmath>
#else
  #include <math.h>
#endif

typedef float REAL;

// mathmatical constants
static const REAL pi = 2 * asin(1.0);  // 
static const REAL e  = exp(1.0);       // Napier's constant

#ifdef __cplusplus
extern "C" {
#endif
void electric_eqn_(REAL *Enx, REAL *Eny, REAL *Enz,
                   REAL *Ex, REAL *Ey, REAL *Ez,
                   REAL *Hx, REAL *Hy, REAL *Hz,
                   REAL *epx, REAL *epy, REAL *epz,
                   REAL *sx, REAL *sy, REAL *sz,
                   REAL *dx, REAL *dy, REAL *dz, REAL *dt,
                   int *imax, int *jmax, int *kmax, int *vc);


void magnetic_eqn_(REAL *Hnx, REAL *Hny, REAL *Hnz,
                   REAL *Hx, REAL *Hy, REAL *Hz,
                   REAL *Ex, REAL *Ey, REAL *Ez,
                   REAL *mx, REAL *my, REAL *mz,
                   REAL *sx, REAL *sy, REAL *sz,
                   REAL *dx, REAL *dy, REAL *dz, REAL *dt,
                   int *imax, int *jmax, int *kmax, int *vc);

void calc_dt_(REAL *dt, REAL *c, REAL *dx, REAL *dy, REAL *dz);

void calc_poynting_(REAL *S, REAL *E, REAL *H,
                    int *imax, int *jmax, int *kmax);

void array_add_(REAL *dst, REAL *src1, REAL *src2, REAL *k, int *sz);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
  #include <3darray.h>
  #include <boundary.h>
#endif

#endif /* __CORE_H__ */


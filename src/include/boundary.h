#ifndef __BOUNDARY_H__
#define __BOUNDARY_H__

#include <3darray.h>

int PEC_boundary_X(s3darray<REAL>& Ex, s3darray<REAL>& Ey, s3darray<REAL>& Ez);

int Mur_boundary(s3darray<REAL>& Enx, s3darray<REAL>& Eny, s3darray<REAL>& Enz,
                 const s3darray<REAL>& Ex, const s3darray<REAL>& Ey, const s3darray<REAL>& Ez,
                 const REAL dx, const REAL dy, const REAL dz,
                 const REAL dt, const REAL c);

#endif // __BOUNDARY_H__


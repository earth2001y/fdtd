#include <core.h>
#include <3darray.h>
#include <boundary.h>

int PEC_boundary(s3darray<REAL>& Enx, s3darray<REAL>& Eny, s3darray<REAL>& Enz)
{
  const size_t imax = Enx.isize();
  const size_t jmax = Enx.jsize();
  const size_t kmax = Enx.ksize();

  // PEC boundary in electric field
  // x-face
  for (size_t k = 0; k < kmax; k++) {
  for (size_t j = 0; j < jmax; j++) {
    Eny(     0,j,k) = 0.0;
    Eny(imax-1,j,k) = 0.0;
    Enz(     0,j,k) = 0.0;
    Enz(imax-1,j,k) = 0.0;
  }}

  // y-face
  for (size_t k = 0; k < kmax; k++) {
  for (size_t i = 0; i < imax; i++) {
    Enx(i,     0,k) = 0.0;
    Enx(i,jmax-1,k) = 0.0;
    Enz(i,     0,k) = 0.0;
    Enz(i,jmax-1,k) = 0.0;
  }}

  // z-face
  for (size_t j = 0; j < jmax; j++) {
  for (size_t i = 0; i < imax; i++) {
    Enx(i,j,     0) = 0.0;
    Enx(i,j,kmax-1) = 0.0;
    Eny(i,j,     0) = 0.0;
    Eny(i,j,kmax-1) = 0.0;
  }}

  return 1;
}


int Mur_boundary(s3darray<REAL>& Enx, s3darray<REAL>& Eny, s3darray<REAL>& Enz,
                 const s3darray<REAL>& Ex, const s3darray<REAL>& Ey, const s3darray<REAL>& Ez,
                 const REAL dx, const REAL dy, const REAL dz,
                 const REAL dt, const REAL c)
{
  //
  // Mur boundary: Enx(i,j,0) = Ex(i,j,1) + ((c*dt - dx) / (c*dt + dx)) * (Enx(i,j,1) - Ex(i,j,0))
  //

  const size_t imax = Ex.isize();
  const size_t jmax = Ex.jsize();
  const size_t kmax = Ex.ksize();

  const REAL mcx = (c * dt - dx) / (c * dt + dx);
  const REAL mcy = (c * dt - dy) / (c * dt + dy);
  const REAL mcz = (c * dt - dz) / (c * dt + dz);

  // x-face
  for (size_t k = 0; k < kmax; k++) {
  for (size_t j = 0; j < jmax; j++) {
    Eny(0,j,k)      = 0.0;
    Enz(0,j,k)      = 0.0;
    Eny(imax-1,j,k) = 0.0;
    Enz(imax-1,j,k) = 0.0;
  }}
  // y-face
  for (size_t k = 0; k < kmax; k++) {
  for (size_t i = 0; i < imax; i++) {
    Enx(i,0,k)      = 0.0;
    Enz(i,0,k)      = 0.0;
    Enx(i,jmax-1,k) = 0.0;
    Enz(i,jmax-1,k) = 0.0;
  }}
  // z-face
  for (size_t j = 0; j < jmax; j++) {
  for (size_t i = 0; i < imax; i++) {
    Enx(i,j,0)      = 0.0;
    Eny(i,j,0)      = 0.0;
    Enx(i,j,kmax-1) = 0.0;
    Eny(i,j,kmax-1) = 0.0;
  }}

  // x-face
  for (size_t k = 0; k < kmax; k++) {
  for (size_t j = 0; j < jmax; j++) {
    Eny(0,j,k)      = Ey(1,j,k) + mcx * (Eny(1,j,k) - Ey(0,j,k) );
    Enz(0,j,k)      = Ez(1,j,k) + mcx * (Enz(1,j,k) - Ez(0,j,k) );
    Eny(imax-1,j,k) = Ey(imax-2,j,k) + mcx * (Eny(imax-2,j,k) - Ey(imax-1,j,k) );
    Enz(imax-1,j,k) = Ez(imax-2,j,k) + mcx * (Enz(imax-2,j,k) - Ez(imax-1,j,k) );
  }}

  // y-face
  for (size_t k = 0; k < kmax; k++) {
  for (size_t i = 0; i < imax; i++) {
    Enx(i,0,k)      = Ex(i,1,k) + mcy * (Enx(i,1,k) - Ex(i,0,k));
    Enz(i,0,k)      = Ez(i,1,k) + mcy * (Enz(i,1,k) - Ez(i,0,k));
    Enx(i,jmax-1,k) = Ex(i,jmax-2,k) + mcy * (Enx(i,jmax-2,k)- Ex(i,jmax-1,k));
    Enz(i,jmax-1,k) = Ez(i,jmax-2,k) + mcy * (Enz(i,jmax-2,k)- Ez(i,jmax-1,k));
  }}

  // z-face
  for (size_t j = 0; j < jmax; j++) {
  for (size_t i = 0; i < imax; i++) {
    Enx(i,j,0)      = Ex(i,j,1) + mcz * (Enx(i,j,1) - Ex(i,j,0));
    Eny(i,j,0)      = Ey(i,j,1) + mcz * (Eny(i,j,1) - Ey(i,j,0));
    Enx(i,j,kmax-1) = Ex(i,j,kmax-2) + mcz * (Enx(i,j,kmax-2) - Ex(i,j,kmax-1));
    Eny(i,j,kmax-1) = Ey(i,j,kmax-2) + mcz * (Eny(i,j,kmax-2) - Ey(i,j,kmax-1));
  }}

  return 1;
}


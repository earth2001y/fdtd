
#include <core.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <iostream>

// physical constant values
REAL c;
REAL ep  = 1.0; // 8.85418782e-12;
REAL mu  = 1.0; // 1.25663706e-6;
REAL esg = 1.0;
REAL msg = 0.0;

REAL dx, dy, dz;
REAL dt;

v3darray<REAL> *E, *H, *Hn, *Ho, *S, *es;
s3darray<REAL> *Enx, *Eny, *Enz;
s3darray<REAL> *Ex,  *Ey,  *Ez;
s3darray<REAL> *Hnx, *Hny, *Hnz;
s3darray<REAL> *Hx,  *Hy,  *Hz;
s3darray<REAL> *epx, *epy, *epz;
s3darray<REAL> *mux, *muy, *muz;
s3darray<REAL> *esx, *esy, *esz;
s3darray<REAL> *msx, *msy, *msz;

int dealloc()
{
  return 1;
}

int alloc(const size_t imax, const size_t jmax, const size_t kmax)
{
  E  = new v3darray<REAL>(imax, jmax, kmax);
  H  = new v3darray<REAL>(imax, jmax, kmax);
  Ho = new v3darray<REAL>(imax, jmax, kmax);
  Hn = new v3darray<REAL>(imax, jmax, kmax);
  S  = new v3darray<REAL>(imax, jmax, kmax);
  es = new v3darray<REAL>(imax, jmax, kmax);

  Hnx = new s3darray<REAL>(imax, jmax, kmax);
  Hny = new s3darray<REAL>(imax, jmax, kmax);
  Hnz = new s3darray<REAL>(imax, jmax, kmax);
  Enx = new s3darray<REAL>(imax, jmax, kmax);
  Eny = new s3darray<REAL>(imax, jmax, kmax);
  Enz = new s3darray<REAL>(imax, jmax, kmax);
  Hx  = new s3darray<REAL>(imax, jmax, kmax);
  Hy  = new s3darray<REAL>(imax, jmax, kmax);
  Hz  = new s3darray<REAL>(imax, jmax, kmax);
  Ex  = new s3darray<REAL>(imax, jmax, kmax);
  Ey  = new s3darray<REAL>(imax, jmax, kmax);
  Ez  = new s3darray<REAL>(imax, jmax, kmax);
  epx = new s3darray<REAL>(imax, jmax, kmax);
  epy = new s3darray<REAL>(imax, jmax, kmax);
  epz = new s3darray<REAL>(imax, jmax, kmax);
  mux = new s3darray<REAL>(imax, jmax, kmax);
  muy = new s3darray<REAL>(imax, jmax, kmax);
  muz = new s3darray<REAL>(imax, jmax, kmax);
  esx = new s3darray<REAL>(imax, jmax, kmax);
  esy = new s3darray<REAL>(imax, jmax, kmax);
  esz = new s3darray<REAL>(imax, jmax, kmax);
  msx = new s3darray<REAL>(imax, jmax, kmax);
  msy = new s3darray<REAL>(imax, jmax, kmax);
  msz = new s3darray<REAL>(imax, jmax, kmax);

  return 1;
}

int test_calc_dt(int imax, int jmax, int kmax)
{
  dx = 1e-3;
  dy = 1e-3;
  dz = 1e-3;

  c = 1.0 / sqrt(ep * mu);
  calc_dt_(&dt, &c, &dx, &dy, &dz);
  printf("pi = %1.8e\n", pi);
  printf("e  = %1.8e\n", e);
  printf("c  = %1.8e\n", c);
  printf("dt = %1.8e\n", dt);
  return 1;
}

void init_fields()
{
  Ex->set(0.0);
  Ey->set(0.0);
  Ez->set(0.0);
  Hx->set(0.0);
  Hy->set(0.0);
  Hz->set(0.0);
  Enx->set(0.0);
  Eny->set(0.0);
  Enz->set(0.0);
  Hnx->set(0.0);
  Hny->set(0.0);
  Hnz->set(0.0);
  mux->set(mu);
  muy->set(mu);
  muz->set(mu);
  epx->set(ep);
  epy->set(ep);
  epz->set(ep);
  msx->set(msg);
  msy->set(msg);
  msz->set(msg);
  esx->set(0.0);
  esy->set(0.0);
  esz->set(0.0);
}

int e2m(int imax, int jmax, int kmax, int step, void (*loc_init_field)(), void (*bind_e)(int), const char* name)
{
  int vc   = 1;
  size_t size = imax * jmax * kmax;

  alloc(imax, jmax, kmax);
  init_fields();

  loc_init_field();

  for (int s = 0; s < step; s++) {

    {
      char fname[256];

      E->reduction(*Ex,*Ey,*Ez);
      H->reduction(*Hx,*Hy,*Hz);
      es->reduction(*epx,*epy,*epz);

      sprintf(fname, "./sph/%s_e_%010d.sph", name, s);
      E->outsph(fname,1.0,1.0,1.0,s,s*dt,1);

      sprintf(fname, "./sph/%s_h_%010d.sph", name, s);
      H->outsph(fname,1.0,1.0,1.0,s,s*dt,1);

//      sprintf(fname, "./sph/%s_s_%010d.sph", name, s);
//      S->outsph(fname,1.0,1.0,1.0,s,s*dt,1);

//      sprintf(fname,"./sph/%s_es_%010d.sph", name, s);
//      es->outsph(fname,1.0,1.0,1.0,s,s*dt,0);
    }


    // magnetic field
    magnetic_eqn_(*Hnx,*Hny,*Hnz,*Hx,*Hy,*Hz,*Ex,*Ey,*Ez,
            *mux,*muy,*muz,*msx,*msy,*msz,
            &dx,&dy,&dz,&dt,&imax,&jmax,&kmax,&vc);

//    REAL p = 0.5;
//    int sz = size;
//    array_add_(*Ho,*Hn,*H,&p,&sz);
//
//    /* Poynting */
//    calc_poynting_(*S,*E,*Ho,&imax,&jmax,&kmax);

    Hx->copy(Hnx);
    Hy->copy(Hny);
    Hz->copy(Hnz);

    // electric field
    electric_eqn_(*Enx,*Eny,*Enz,*Ex,*Ey,*Ez,*Hx,*Hy,*Hz,
            *epx,*epy,*epz,*esx,*esy,*esz,
            &dx,&dy,&dz,&dt,&imax,&jmax,&kmax,&vc);


//    PEC_boundary(*Enx,*Eny,*Enz);
    Mur_boundary(*Enx,*Eny,*Enz,*Ex,*Ey,*Ez,dx,dy,dz,dt,c);

    bind_e(s);

    Ex->copy(Enx);
    Ey->copy(Eny);
    Ez->copy(Enz);
  }

  dealloc();

  return 1;
}

//
void case1_init_field()  { Ez->set(1.0); }
void case1_bind_e(int s) { }

//
void case2_init_field()  {
  int imax = Enz->isize();
  int jmax = Enz->jsize();
  int kmax = Enz->ksize();
  for (int k = 0; k < kmax; k++) {
    int i = (imax / 2) - 1;
    int j = (jmax / 2) - 1;
    (*Ez)(i  ,j  ,k) = 1.0;
    (*Ez)(i-1,j  ,k) = 1.0;
    (*Ez)(i  ,j-1,k) = 1.0;
    (*Ez)(i-1,j-1,k) = 1.0;
  }
}
void case2_bind_e(int s) {
  int imax = Enz->isize();
  int jmax = Enz->jsize();
  int kmax = Enz->ksize();
  for (int k = 0; k < kmax; k++) {
    int i = (imax / 2) - 1;
    int j = (jmax / 2) - 1;
    (*Enz)(i  ,j  ,k) = 1.0;
    (*Enz)(i-1,j  ,k) = 1.0;
    (*Enz)(i  ,j-1,k) = 1.0;
    (*Enz)(i-1,j-1,k) = 1.0;
  }
}

//
void case3_init_field()  { }
void case3_bind_e(int s) {
  int imax = Enz->isize();
  int jmax = Enz->jsize();
  int kmax = Enz->ksize();
  for (int k = 0; k < kmax; k++) {
    int i = (imax / 2) - 1;
    int j = (jmax / 2) - 1;
    (*Enz)(i  ,j  ,k) = sin((1.0/16.0) * pi * s);
    (*Enz)(i-1,j  ,k) = sin((1.0/16.0) * pi * s);
    (*Enz)(i  ,j-1,k) = sin((1.0/16.0) * pi * s);
    (*Enz)(i-1,j-1,k) = sin((1.0/16.0) * pi * s);
  }
}

//
void case4_init_field()  {
  int imax = Enz->isize();
  int jmax = Enz->jsize();
  int kmax = Enz->ksize();
  for (int k = 0; k < kmax; k++) {
  for (int j = 0; j < jmax; j++) {
  for (int i = 0; i < imax; i++) {
    if (i + (j/4) - (imax/8) < imax / 2) {
      (*epz)(i,j,k) = 2.0 * ep;
      (*epy)(i,j,k) = 2.0 * ep;
      (*epx)(i,j,k) = 2.0 * ep;
      (*muz)(i,j,k) = 2.0 * mu;
      (*muy)(i,j,k) = 2.0 * mu;
      (*mux)(i,j,k) = 2.0 * mu;
    }
  }}}
}
void case4_bind_e(int s) {
  int imax = Enz->isize();
  int jmax = Enz->jsize();
  int kmax = Enz->ksize();
  for (int k = 0; k < kmax; k++) {
  for (int j = 0; j < jmax; j++) {
    (*Enz)(0,j,k) = sin((1.0/16.0) * pi * s);
  }}
}

int main(int argc, char *argv[])
{
  test_calc_dt(100,100,3);
  e2m(20, 20, 3, 100, case1_init_field, case1_bind_e, "c1");   // non-rot
  e2m(100,100,3, 200, case2_init_field, case2_bind_e, "c2");   // Ampare's raw
  e2m(100,100,3, 200, case3_init_field, case3_bind_e, "c3");   // Radiation
  e2m(100,100,3, 300, case4_init_field, case4_bind_e, "c4");   // Optical

  return 0;
}


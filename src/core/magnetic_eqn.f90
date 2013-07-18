
subroutine magnetic_eqn(Hnx,Hny,Hnz, &
                        Hx,Hy,Hz,Ex,Ey,Ez, &
                        mx,my,mz,sx,sy,sz, &
                        dx,dy,dz,dt, &
                        imax,jmax,kmax,vc)
  implicit none

  integer,intent(in):: imax,jmax,kmax,vc

  real,intent(out):: Hnx(imax,jmax,kmax)
  real,intent(out):: Hny(imax,jmax,kmax)
  real,intent(out):: Hnz(imax,jmax,kmax)

  real,intent(in):: Hx(imax,jmax,kmax)
  real,intent(in):: Hy(imax,jmax,kmax)
  real,intent(in):: Hz(imax,jmax,kmax)

  real,intent(in):: sx(imax,jmax,kmax)
  real,intent(in):: sy(imax,jmax,kmax)
  real,intent(in):: sz(imax,jmax,kmax)
  real,intent(in):: mx(imax,jmax,kmax)
  real,intent(in):: my(imax,jmax,kmax)
  real,intent(in):: mz(imax,jmax,kmax)

  real,intent(in):: Ex(imax,jmax,kmax)
  real,intent(in):: Ey(imax,jmax,kmax)
  real,intent(in):: Ez(imax,jmax,kmax)

  real,intent(in):: dx,dy,dz,dt

  real:: dxi,dyi,dzi,dti,i2
  integer:: i, j, k

  dxi = 1.d0 / dx
  dyi = 1.d0 / dy
  dzi = 1.d0 / dz
  dti = 1.d0 / dt
  i2  = 1.d0 / 2.d0

  do k=1,kmax-1
  do j=1,jmax-1
  do i=1,imax
    Hnx(i,j,k) = ( (mx(i,j,k) * dti - sx(i,j,k) * i2) * Hx(i,j,k)  &
                 - (Ez(i,j+1,k) - Ez(i,j,k)) * dyi                 &
                 + (Ey(i,j,k+1) - Ey(i,j,k)) * dzi                 &
                 ) / (mx(i,j,k) * dti + sx(i,j,k) * i2)
  enddo
  enddo
  enddo

  do k=1,kmax-1
  do j=1,jmax
  do i=1,imax-1
    Hny(i,j,k) = ( (my(i,j,k) * dti - sy(i,j,k) * i2) * Hy(i,j,k)  &
                 - (Ex(i,j,k+1) - Ex(i,j,k)) * dzi                 &
                 + (Ez(i+1,j,k) - Ez(i,j,k)) * dxi                 &
                 ) / (my(i,j,k) * dti + sy(i,j,k) * i2)
  enddo
  enddo
  enddo

  do k=1,kmax
  do j=1,jmax-1
  do i=1,imax-1
    Hnz(i,j,k) = ( (mz(i,j,k) * dti - sz(i,j,k) * i2) * Hz(i,j,k)  &
                 - (Ey(i+1,j,k) - Ey(i,j,k)) * dxi                 &
                 + (Ex(i,j+1,k) - Ex(i,j,k)) * dyi                 &
                 ) / (mz(i,j,k) * dti + sz(i,j,k) * i2)
  enddo
  enddo
  enddo

end subroutine magnetic_eqn


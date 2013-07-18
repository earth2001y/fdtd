
subroutine electric_eqn(Enx,Eny,Enz, &
                        Ex,Ey,Ez,Hx,Hy,Hz, &
                        epx,epy,epz,sx,sy,sz, &
                        dx,dy,dz,dt, &
                        imax,jmax,kmax,vc)
  implicit none

  integer,intent(in):: imax,jmax,kmax,vc

  real,intent(out):: Enx(imax,jmax,kmax)
  real,intent(out):: Eny(imax,jmax,kmax)
  real,intent(out):: Enz(imax,jmax,kmax)

  real,intent(in):: Ex(imax,jmax,kmax)
  real,intent(in):: Ey(imax,jmax,kmax)
  real,intent(in):: Ez(imax,jmax,kmax)

  real,intent(in):: sx(imax,jmax,kmax)
  real,intent(in):: sy(imax,jmax,kmax)
  real,intent(in):: sz(imax,jmax,kmax)
  real,intent(in):: epx(imax,jmax,kmax)
  real,intent(in):: epy(imax,jmax,kmax)
  real,intent(in):: epz(imax,jmax,kmax)

  real,intent(in):: Hx(imax,jmax,kmax)
  real,intent(in):: Hy(imax,jmax,kmax)
  real,intent(in):: Hz(imax,jmax,kmax)

  real,intent(in):: dx,dy,dz,dt

  real:: dxi,dyi,dzi,dti,i2
  integer:: i,j,k

  dxi = 1.d0 / dx
  dyi = 1.d0 / dy
  dzi = 1.d0 / dz
  dti = 1.d0 / dt
  i2  = 1.d0 / 2.d0

  Enx = 0.0
  Eny = 0.0
  Enz = 0.0

  do k=2,kmax
  do j=2,jmax
  do i=1,imax
    Enx(i,j,k) = ( (epx(i,j,k) * dti - sx(i,j,k) * i2) * Ex(i,j,k) &
                 + (Hz(i,j,k) - Hz(i,j-1,k)) * dyi                 &
                 - (Hy(i,j,k) - Hy(i,j,k-1)) * dzi                 &
                 ) / (epx(i,j,k) * dti + sx(i,j,k) * i2)
  enddo
  enddo
  enddo

  do k=2,kmax
  do j=1,jmax
  do i=2,imax
    Eny(i,j,k) = ( (epy(i,j,k) * dti - sy(i,j,k) * i2) * Ey(i,j,k) &
                 + (Hx(i,j,k) - Hx(i,j,k-1)) * dzi                 &
                 - (Hz(i,j,k) - Hz(i-1,j,k)) * dxi                 &
                 ) / (epy(i,j,k) * dti + sy(i,j,k) * i2)
  enddo
  enddo
  enddo

  do k=1,kmax
  do j=2,jmax
  do i=2,imax
    Enz(i,j,k) = ( (epz(i,j,k) * dti - sz(i,j,k) * i2) * Ez(i,j,k) &
                 + (Hy(i,j,k) - Hy(i-1,j,k)) * dxi                 &
                 - (Hx(i,j,k) - Hx(i,j-1,k)) * dyi                 &
                 ) / (epz(i,j,k) * dti + sz(i,j,k) * i2)
  enddo
  enddo
  enddo

end subroutine electric_eqn


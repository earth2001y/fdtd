
subroutine H2Hxyz(H, Hx, Hy, Hz, &
                  imax, jmax, kmax, vc)

  implicit none

  integer,intent(in):: imax, jmax, kmax, vc

  real,intent(in)::  H(3,imax,jmax,kmax)
  real,intent(out):: Hx(imax,jmax,kmax)
  real,intent(out):: Hy(imax,jmax,kmax)
  real,intent(out):: Hz(imax,jmax,kmax)

  integer:: i,j,k

  do k=1,kmax
  do j=1,jmax
  do i=1,imax-1
    Hx(i,j,k) = (H(1,i,j,k) + H(1,i+1,j,k)) * 0.d5
  enddo
  enddo
  enddo

  do k=1,kmax
  do j=1,jmax-1
  do i=1,imax
    Hy(i,j,k) = (H(2,i,j,k) + H(2,i,j+1,k)) * 0.d5
  enddo
  enddo
  enddo

  do k=1,kmax
  do j=1,jmax-1
  do i=1,imax
    Hz(i,j,k) = (H(3,i,j,k) + H(3,i,j,k+1)) * 0.d5
  enddo
  enddo
  enddo

end subroutine H2Hxyz

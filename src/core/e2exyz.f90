
subroutine E2Exyz(E, Ex, Ey, Ez, &
                  imax, jmax, kmax, vc)

  implicit none

  integer,intent(in):: imax, jmax, kmax, vc

  real,intent(in)::  E(3,imax,jmax,kmax)
  real,intent(out):: Ex(imax,jmax,kmax)
  real,intent(out):: Ey(imax,jmax,kmax)
  real,intent(out):: Ez(imax,jmax,kmax)

  integer:: i,j,k

  do k=1,kmax
  do j=1,jmax
  do i=1,imax-1
    Ex(i,j,k) = (E(1,i,j,k) + E(1,i+1,j,k)) * 0.d5
  enddo
  enddo
  enddo

  do k=1,kmax
  do j=1,jmax-1
  do i=1,imax
    Ey(i,j,k) = (E(2,i,j,k) + E(2,i,j+1,k)) * 0.d5
  enddo
  enddo
  enddo

  do k=1,kmax-1
  do j=1,jmax
  do i=1,imax
    Ez(i,j,k) = (E(3,i,j,k) + E(3,i,j,k+1)) * 0.d5
  enddo
  enddo
  enddo

end subroutine E2Exyz

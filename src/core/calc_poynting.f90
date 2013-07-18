subroutine calc_poynting(S, E, H, imax, jmax, kmax)
    implicit none
    integer, intent(in):: imax, jmax, kmax
    real, intent(in):: E(3,imax,jmax,kmax)
    real, intent(in):: H(3,imax,jmax,kmax)
    real, intent(out):: S(3,imax,jmax,kmax)

    integer i, j, k

    do k = 1,kmax
    do j = 1,jmax
    do i = 1,imax
        S(1,i,j,k) = E(2,i,j,k) * H(3,i,j,k) - E(3,i,j,k) * H(2,i,j,k)
        S(2,i,j,k) = E(3,i,j,k) * H(1,i,j,k) - E(1,i,j,k) * H(3,i,j,k)
        S(3,i,j,k) = E(1,i,j,k) * H(2,i,j,k) - E(2,i,j,k) * H(1,i,j,k)
    enddo
    enddo
    enddo

end subroutine calc_poynting


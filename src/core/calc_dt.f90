
subroutine calc_dt(dt,c,dx,dy,dz)

  implicit none

  real,intent(out):: dt

  real,intent(in):: c
  real,intent(in):: dx, dy, dz

  dt = 1.d0 / (sqrt(1.d0/(dx**2) + 1.d0/(dy**2) + 1.d0/(dz**2)) * c)

end subroutine calc_dt


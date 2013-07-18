subroutine array_add(dst, src1, src2, k, sz)
  implicit none

  integer,intent(in):: sz
  real,intent(in)::  src1(sz), src2(sz), k
  real,intent(out):: dst(sz)

  dst = (src1 + src2) * k

end subroutine array_add


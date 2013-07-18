program e2m_xy
  implicit none

  integer:: imax = 128
  integer:: jmax = 128
  integer:: kmax = 1
  integer:: vc = 1

  real,allocatable:: ep(:,:,:), epx(:,:,:), epy(:,:,:), epz(:,:,:)
  real,allocatable:: sg(:,:,:), sgx(:,:,:), sgy(:,:,:), sgz(:,:,:)
  real,allocatable:: mu(:,:,:), mux(:,:,:), muy(:,:,:), muz(:,:,:)
  real,allocatable:: ss(:,:,:), ssx(:,:,:), ssy(:,:,:), ssz(:,:,:)

  real,allocatable:: E(:,:,:,:), En(:,:,:,:)
  real,allocatable:: H(:,:,:,:), Hn(:,:,:,:)

  allocate(E(3,nx,ny,nz), En(3,nx,ny,nz))
  allocate(H(3,nx,ny,nz), Hn(3,nx,ny,nz))
  allocate(ep(nx,ny,nz), epx(nx,ny,nz), epy(nx,ny,nz), epz(nx,ny,nz))
  allocate(sg(nx,ny,nz), sgx(nx,ny,nz), sgy(nx,ny,nz), sgz(nx,ny,nz))
  allocate(mu(nx,ny,nz), mux(nx,ny,nz), muy(nx,ny,nz), muz(nx,ny,nz))
  allocate(ss(nx,ny,nz), ssx(nx,ny,nz), ssy(nx,ny,nz), ssz(nx,ny,nz))

  deallocate(E, En)
  deallocate(H, Hn)
  deallocate(ep, epx, epy, epz)
  deallocate(sg, sgx, sgy, sgz)
  deallocate(mu, mux, muy, muz)
  deallocate(ss, ssx, ssy, ssz)

end program e2m



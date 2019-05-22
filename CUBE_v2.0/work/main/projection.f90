#define proj3d
! single node, fine cell, NGP projection
subroutine projection
  use variables
  implicit none
  save
  integer(8) idxf(3)
  real proj_yz(nf,nf), proj_xz(nf,nf), proj_xy(nf,nf)
#ifdef proj3d
    real den(nf,nf,nf)
#endif

  if (head) print*, 'projection'

  proj_yz=0
  proj_xz=0
  proj_xy=0
#ifdef proj3d
    den=0
#endif

  ip=0

  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
  do k=1,nt
  do j=1,nt
  do i=1,nt
    np=rhoc(i,j,k,itx,ity,itz)
    do l=1,np
      ip=ip+1
      !idxf=ceiling( nft*((/itx,ity,itz/)-1) &
      !             +4.0*((/i,j,k/)-1) &
      !             +0.015625*((x(:,ip)+int(-128,1))+128.5) )
      idxf=ceiling(nft*((/itx,ity,itz/)-1)+ncell*((/i,j,k/)-1) &
                   +(int(xp(:,ip)+ishift,izipx)+rshift)*ncell*x_resolution)
      proj_yz(idxf(2),idxf(3))=proj_yz(idxf(2),idxf(3))+1
      proj_xz(idxf(1),idxf(3))=proj_xz(idxf(1),idxf(3))+1
      proj_xy(idxf(1),idxf(2))=proj_xy(idxf(1),idxf(2))+1
#ifdef proj3d
      den(idxf(1),idxf(2),idxf(3))=den(idxf(1),idxf(2),idxf(3))+1
#endif
    enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  sync all

  open(11,file=output_name('den'),status='replace',access='stream')
  write(11) den
  close(11)


  !open(11,file='./output/proj_yz.dat',status='replace',access='stream')
  !write(11) proj_yz
  !close(11)

  !open(12,file='./output/proj_xz.dat',status='replace',access='stream')
  !write(12) proj_xz
  !close(12)

  !open(13,file='./output/proj_xy.dat',status='replace',access='stream')
  !write(13) proj_xy
  !close(13)

  !#ifdef proj3d
  !open(14,file='./output/proj_3d.dat',status='replace',access='stream')
  !write(14) proj_3d
  !close(14)
  !#endif
  sync all

endsubroutine

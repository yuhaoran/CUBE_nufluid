program zipconvert
  use parameters
  implicit none
  save

  integer :: i,j,k,l,nplocal,itx,ity,itz
  integer(8) :: nlast,ip,np

  real(4) :: pos(3),vel(3)

  integer(izipx),allocatable :: xp(:,:)
  integer(izipv),allocatable :: vp(:,:)

  integer(4) :: rhoc(nt,nt,nt,nnt,nnt,nnt)

  real, dimension(6,256**3) :: xv

  call geometry

  open(11,file=output_name('info'),status='old',action='read',access='stream')
  read(11) sim
  close(11)
  !call print_header(sim); stop
  if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
     print*, 'zip format incompatable'
     close(11)
     stop
  endif
  nplocal=sim%nplocal
  print*, 'nplocal =',nplocal
  allocate(xp(3,nplocal),vp(3,nplocal))

  cur_checkpoint=1
  write(*,*) 'output_name',output_name('xp')
  open(11,file='../output/image1/100.000_xp_1.bin',status='old',action='read',access='stream')
  read(11) xp
  close(11)
  open(11,file='../output/image1/100.000_np_1.bin',status='old',action='read',access='stream')
  read(11) rhoc
  close(11)
  open(11,file='../output/image1/100.000_vp_1.bin',status='old',action='read',access='stream')
  read(11) vp
  close(11)

  nlast=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
     do k=1,nt
     do j=1,nt
     do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
           ip=nlast+l
           
           !positions
           pos=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
           pos=pos*ncell !fine cells
           write(*,*) 'pos',pos
           xv(1:3,ip) = pos

           !velocities
           vel=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sim%sigma_vi*vrel_boost))
           write(*,*) 'vel',vel
           xv(4:6,ip) = vel

        end do
        nlast=nlast+np
        if (nlast.gt.10) goto 10
     end do
     end do
     end do

  end do
  end do
  end do

10 write(*,*) 'fin'

end program zipconvert

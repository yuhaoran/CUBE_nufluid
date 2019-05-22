program initial_conditions_nu
  use omp_lib
  use pencil_fft
  implicit none
  save

  !Fermi-Dirac CDF
  integer, parameter :: ncdf = 10000
  real, dimension(2,ncdf) :: cdf

  !Units
  !real, parameter :: vp2s = 1.0/(300.*sqrt(omega_m)*box/a_i_nu/2./nc)

  !real, parameter :: vp2s = 1.0/(150.*h0*sqrt(omega_m)*box/a_i_nu/nf_global)
  real, parameter :: vsim2phys=150.*sqrt(omega_m)*box/a_i_nu/2./nc ! previous version. missed h0 and nf_global
  !real, parameter :: vsim2phys=(150./a_i_nu)*box*h0*sqrt(omega_m)/nf_global
  ! because: sim%vsim2phys=(150./a)*box*h0*sqrt(omega_m)/nf_global in initial_conditions.f90
  real, parameter :: fdf = 25.8341 !kT/m for T=1K, m=1eV
  real, parameter :: fd = (1./vsim2phys)*fdf*maxval(Tnu/Mnu)/a_i_nu !kBcTnu/mass with temp in K and mass in eV
  real, parameter :: sigma_vi_nu = 3.59714*fd !fd velocity dispersion (45/3 * Zeta(5)/Zeta(3))**0.5

  !Seed
  integer(4) seedsize
  integer(4), allocatable :: iseed(:)
  real, allocatable :: rseed_all(:,:)

  !Useful variables and small arrays
  integer :: itx,ity,itz,i,j,k,ii,jj,kk,l,n,pii,pjj,pkk,b,b1,b2
  integer(8) ip

  !Particle information
  integer(8), parameter :: npt = np_nc_nu*nt
  integer(8), parameter :: npmax=npt**3
  integer(izipx_nu) xp(3,npmax)
  integer(izipv_nu) vp(3,npmax)
  integer(izipi), dimension(npt**3) :: pid
  real, dimension(3) :: xq,vq,rng,vmax

  integer(4), dimension(nt,nt,nt) :: rhoc
  real(4), dimension(3,nt,nt,nt) :: vfield

  real(8) vreal(3)
  real(8) std_vsim_c,std_vsim_res,std_vsim

  !Setup
  call geometry
  if (head) then
     write(*,*) ''
     write(*,*) 'Homogeneous Initial Conditions for Massive Neutrinos'
     print*, 'on',nn**3,' images'
     print*, 'Resolution', ng*nn
     print*, 'To genterate npglobal_nu =',int(np_nc_nu*nc*nn,2),'^3'
     print*, 'Box size', box
     print*, 'output: ', opath
     write(*,*) ''
     write(*,*) 'vp2s/(km/s)=',(1./vsim2phys)
     write(*,*) 'fd/(km/s)=',fd*vsim2phys
  end if
  sync all
  call system('mkdir -p '//opath//'image'//image2str(image))

  !Read seed
  if (head) write(*,*) 'Reading seeds'
  call random_seed(size=seedsize)
  seedsize=max(seedsize,36)
  allocate(iseed(seedsize))
  allocate(rseed_all(seedsize,nn**3))
  open(11,file=output_dir()//'seed'//output_suffix(),status='old',access='stream')
  read(11) iseed
  close(11)
  call random_seed(put=iseed)

!!$  !Read cdf table
!!$  if (head) write(*,*) 'Reading cdf'
!!$  open(11,file='./CDFTable.txt')
!!$  read(11,*) cdf
!!$  close(11)
!!$  cdf(1,:) = cdf(1,:)*fd

  if (head) write(*,*) 'Generating cdf'
  call compute_cdf
  if (head) then
     write(*,*) 'Writing CDF to file'
     write(*,*) 'FD factor used: ',fd
     open(11,file=opath//'cdf.txt')
     do i=1,ncdf
        write(11,*) cdf(1,i)/fd,cdf(2,i)
     end do
     close(11)
  end if

  !Create particles
  if (head) write(*,*) 'Computing particle positions and velocities'
  open(unit=11,file=ic_name_nu('xp_nu'),status='replace',access='stream')
  open(unit=12,file=ic_name_nu('vp_nu'),status='replace',access='stream')
  open(unit=13,file=ic_name_nu('np_nu'),status='replace',access='stream')
  open(unit=14,file=ic_name_nu('vc_nu'),status='replace',access='stream')
  open(unit=15,file=ic_name_nu('id_nu'),status='replace',access='stream')
  vfield=0
  rhoc=np_nc_nu**3
  vmax=0
  std_vsim_c=0; std_vsim_res=0; std_vsim=0;
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do kk=1,nt
    do jj=1,nt
    do ii=1,nt
      do pkk=1,np_nc_nu
      do pjj=1,np_nc_nu
      do pii=1,np_nc_nu
        ip=pii+np_nc_nu*(pjj-1)+np_nc_nu**2*(pkk-1) &
        +np_nc_nu**3*(ii-1)+np_nc_nu**3*nt*(jj-1)+np_nc_nu**3*nt**2*(kk-1)
        xq=([itx,ity,itz]-1)*nt+([ii,jj,kk]-1) + (([pii,pjj,pkk]-1d0)/np_nc_nu)+0.5d0/ncell
        xp(:,ip)=floor( xq/x_resolution_nu,kind=8 )
        !print*,ip,xq(1),xp(1,ip)
        !Random velocities
        call random_number(rng)
        !Bisection interpolate CDF
        !Type out fns here to save time in case compiler does not inline
        b1=1
        b2=ncdf
        do while(b2-b1>1)
          b=(b1+b2)/2
          if ( rng(1)>cdf(2,b) ) then
            b1=b
          else
            b2=b
          endif
        enddo
        n=merge(b1,b2,b1<b2)
        rng(1)=(cdf(1,n)*(cdf(2,n+1)-rng(1))+cdf(1,n+1)*(rng(1)-cdf(2,n)))/(cdf(2,n+1)-cdf(2,n))

        !!Store fraction of max velocity
        !!Fermi-Dirac CDF Approximated by Gaussian
        pid(ip)=nint(approxCDF(rng(1))*int(2,8)**(8*izipi)-int(2,8)**(8*izipi-1),kind=izipi)

        !!Amplitude and Angle
        rng(1)=rng(1)!*fd !!Holds velocity amplitude
        rng(2)=2.*rng(2)-1. !cosTheta in (-1,1)
        rng(3)=rng(3)*2.*pi !Phi in 0 to 2*pi
        !!Direction
        vq(1)=rng(1)*sqrt(1.-rng(2)**2.)*cos(rng(3))
        vq(2)=rng(1)*sqrt(1.-rng(2)**2.)*sin(rng(3))
        vq(3)=rng(1)*rng(2)

        vmax=max(vmax,abs(vq))

        vp(:,ip)=nint(real(nvbin_nu-1)*atan(sqrt(pi/2)/(sigma_vi_nu*vrel_boost)*vq)/pi,kind=izipv_nu)
      enddo
      enddo
      enddo
    !stop
    enddo
    enddo
    enddo

    write(11) xp
    write(12) vp
    write(13) rhoc
    write(14) vfield
    write(15) pid

    ! velocity analysis
    ip=0
    do k=1,nt
    do j=1,nt
    do i=1,nt
      std_vsim_c=std_vsim_c+sum(vfield(:,i,j,k)**2)
      do l=1,np_nc_nu**3
        ip=ip+1
        vreal=tan(pi*real(vp(:,ip))/real(nvbin_nu-1))/(sqrt(pi/2)/(sigma_vi_nu*vrel_boost))
        std_vsim_res=std_vsim_res+sum(vreal**2)
        vreal=vreal+vfield(:,i,j,k)
        std_vsim=std_vsim+sum(vreal**2)
      enddo
    enddo
    enddo
    enddo

  enddo
  enddo
  enddo ! end of tile loop

  close(11)
  close(12)
  close(13)
  close(14)
  close(15)

  ! update header information
  open(unit=10,file=ic_name('info'),access='stream')
  read(10) sim
  sim%sigma_vi_nu=sigma_vi_nu
  sim%dt_vmax_nu=vbuf*20./maxval(abs(vmax))
  rewind(10)
  write(10) sim
  close(10)

  if (head) then
    print*,''
    print*,'Velocity analysis on head node'
    std_vsim_res=sqrt(std_vsim_res/sim%nplocal_nu)
    std_vsim_c=sqrt(std_vsim_c/nc/nc/nc)
    std_vsim=sqrt(std_vsim/sim%nplocal_nu)
    print*,'  std_vsim         ',real(std_vsim*sim%vsim2phys,4),'km/s'
    print*,'  std_vsim_c       ',real(std_vsim_c*sim%vsim2phys,4),'km/s'
    print*,'  std_vsim_res     ',real(std_vsim_res*sim%vsim2phys,4),'km/s'
    print*,'  std_vi (sim unit)',real(std_vsim_res/sqrt(3.),4),'(simulation unit)'
    print*,''
  endif
  sync all

  call print_header(sim)
  if (head) write(*,*) 'Finished neutrino ic'


contains

  function approxCDF(v) result(c)
    implicit none
    real, intent(in) :: v
    real :: c
    real, parameter :: s=3.5
    c=1.-exp(-(v/s)**2.)
  end

  function invertCDF(c) result(v)
    implicit none
    real, intent(in) :: c
    real :: v
    real, parameter :: s=3.5
    v=s*sqrt(log(1./(1.-c)))
  end

  subroutine compute_cdf
    implicit none
    real(8), dimension(2,ncdf) :: cdf0
    real, dimension(2,ncdf) :: cdfn
    integer, parameter :: ni = 1000 !How many points to integrate per cdf
    real, dimension(ni) :: x,y
    real, parameter :: maxu = 15.0
    real, parameter :: cdfinf = 1.80309
    integer :: i,j,n
    real :: l,u,fnu,fdnu

    cdf0 = 0.
    do i=2,ncdf
       !Limits of integration
       l=maxu*(1.0*i-1.)/ncdf
       u=maxu*(1.0*i)/ncdf
       cdf0(1,i)=u
       !Integral
       do j=1,ni
          x(j)=l+(j-1)*(u-l)/(ni-1)
          y(j)=f0(x(j))
       enddo
       cdf0(2,i)=cdf0(2,i-1)+integrate(x,y)

    enddo

    write(*,*) 'cdf: u->inf = ',cdf0(2,ncdf),cdfinf
    cdf0(2,:) = cdf0(2,:)/cdf0(2,ncdf)

    !Now sum up for each neutrino
    cdf=0
    cdf(1,:)=fd*cdf0(1,:)
    cdfn(2,:)=cdf0(2,:)
    do n=1,Nnu
       fnu=Mnu(n)*(Tnu(n)/Tcnb)**3/Meff ! fraction of energy in this neutrino
       fdnu=(1./vsim2phys)*fdf*Tnu(n)/Mnu(n)/a_i_nu
       cdfn(1,:)=cdf0(1,:)*fdnu
       do i=1,ncdf
          j=nearest_loc(cdf(1,i),cdfn(1,:))
          cdf(2,i) = cdf(2,i)+fnu*interp(cdf(1,i),cdfn(1,j),cdfn(1,j+1),cdfn(2,j),cdfn(2,j+1))
       enddo
    enddo
  endsubroutine compute_cdf

  function f0(u) result(f)
    implicit none
    real, intent(in) :: u
    real :: f
    f= u**2./(exp(u)+1.)
  end

  function integrate(x,y) result(s)
    implicit none
    real, dimension(:), intent(in) :: x,y
    real :: s
    integer :: i
    s=0
    do i=2,size(x)
       s=s+0.5*(x(i)-x(i-1))*(y(i)+y(i-1))
    enddo
  end

  function nearest_loc(u,c) result(nl)
    implicit none
    real, intent(in) :: u
    real, dimension(:) :: c
    integer :: b1,b2,b,nl
    b1=1
    b2=size(c)
    do while(b2-b1>1)
       b=(b1+b2)/2
       if ( u.gt.c(b) ) then
          b1=b
       else
          b2=b
       endif
    enddo
    nl=merge(b1,b2,b1<b2)
  end

  function interp(x,x1,x2,y1,y2) result(y)
    implicit none
    real, intent(in) :: x,x1,x2,y1,y2
    real :: y
    y = (y1*(x2-x)+y2*(x-x1))/(x2-x1)
  end

end

#define q_EQN_OF_STATE
module HydroGR
  use parameters
  use hydro3d
  implicit none

  !Units
  !!External Grid parameters
  integer, parameter :: hg_nf = nf
  integer, parameter :: hg_nf_tot = hg_nf*nn !nf_global
  integer, parameter :: hg_nfb = ncb*ncell!merge(6,2,h_hrt)
  !!Velocity -- always matches cubep3m
  real(kind=h_fpp), parameter :: hg_dv = (2*hg_nf_tot)/(3.*100.*sqrt(omega_m)*box)

  !Verbosity (-1=none, 0=little, 1=more, 2=debug)
  integer :: hg_verb = 0

  type hydro
     !real(kind=h_fpp), allocatable :: fld(:,:,:,:)[*]
     ! coarray conflicts with pointer in hydrodf.f90
     real(kind=h_fpp), allocatable :: fld(:,:,:,:)
     integer :: nc,ncb
     real(kind=h_fpp) :: g,cs2,dx,dv
     logical :: isothermal
   contains
     procedure :: setup => hg_setup
     procedure :: timestep => hg_timestep
     procedure :: evolve => hg_evolve
     procedure :: gravity => hg_gravity
     procedure :: buffer => hg_buffer
     procedure :: density => hg_density
     procedure :: set_fld => hg_set_fld
     procedure :: checkpoint => hg_checkpoint
     procedure :: restart => hg_restart
  end type hydro

contains

  subroutine hg_setup(h,nc,g,cs,iso)
    implicit none
    class(hydro) :: h
    integer, intent(in) :: nc
    real(kind=h_fpp), intent(in), optional :: g,cs
    logical, intent(in), optional :: iso

    call hg_write('Setting up hydro with nc=',0,real(nc,kind=h_fpp))
    h%nc=nc
    h%dx=1.*hg_nf/h%nc
    h%dv=hg_dv
    h%ncb=max( merge(6,2,h_hrt), 1+floor(hg_nfb/h%dx) )

    if (allocated(h%fld)) deallocate(h%fld)
    allocate(h%fld(5,1-h%ncb:h%nc+h%ncb,1-h%ncb:h%nc+h%ncb,1-h%ncb:h%nc+h%ncb))
    print*,'  h%nc =',h%nc
    print*,'  h%ncb =',h%ncb
    print*,'  res =',(h%nc+h%ncb)-(1-h%ncb)+1
    print*,'  size of [h%fld] =',size(h%fld)
    if (present(g)) then
       h%g=g
    else
       h%g=5./3.
    end if
    call hg_write('  >gamma=',1,h%g)

    if (present(cs)) then
       h%cs2=max(1e-6,(cs*h%dv)**2)
    else
       h%cs2=1e-6
    end if
    call hg_write('  >cs (km/s)=',1,sqrt(h%cs2)/h%dv)
    call hg_write('  >cs (sim)=',1,sqrt(h%cs2))

    if (present(iso)) then
       h%isothermal=iso
    else
       h%isothermal=.false.
    end if

    h%fld(1,:,:,:)=1.
    h%fld(2:4,:,:,:)=0.
    h%fld(5,:,:,:)=h%cs2*h%fld(1,:,:,:)/h%g/(h%g-1.)

#ifdef q_EQN_OF_STATE
    h%fld(5,:,:,:)=h%cs2*h%fld(1,:,:,:)/(h%g-1.)
#endif
    
    call hg_write('  >complete',1)
    print*,''

  end subroutine hg_setup

  subroutine hg_timestep(h,dt)
    implicit none
    class(hydro) :: h
    real, intent(inout) :: dt
    real(kind=h_fpp) :: dtn
    
    if (h%isothermal) then
       call h3_timestep(h%fld,h%g,dtn,h%dx,h%cs2)
    else
       call h3_timestep(h%fld,h%g,dtn,h%dx)
    end if
    dt=min(dt,dtn)
    call hg_global_min(dt)

    call hg_write('Timestep=',0,real(dt,kind=h_fpp))

  end subroutine hg_timestep

  subroutine hg_evolve(h,dt,dtg,dir)
    implicit none
    class(hydro) :: h
    real, intent(in) :: dt
    real, intent(inout) :: dtg
    integer, intent(in) :: dir

    call hg_write('Evolving hydro w/ cs=',0,sqrt(h%cs2)/h%dv)

    call h%buffer
    if (h%isothermal) then
       call h3_evolve(h%fld,h%g,real(dt,kind=h_fpp),h%dx,dir,h%cs2)
    else
       call h3_evolve(h%fld,h%g,real(dt,kind=h_fpp),h%dx,dir)
    end if
    call h%buffer
    call h%timestep(dtg)

    call hg_hydro_properties(h)
    
    call hg_write('>Complete',1)

  end subroutine hg_evolve

  !gravity subroutine -- adds gravitational acceleration to momentum and energy
  !gc controls NGP (gc=0) or CIC (gc!=0) interpolation, the latter requires buffered forces
  subroutine hg_gravity(h,x,acc,gc)
    implicit none
    class(hydro) :: h
    real, dimension(3), intent(in) :: x
    real(8), dimension(3), intent(in) :: acc
    integer, intent(in) :: gc
    real(kind=h_fpp), dimension(3) :: accr,xx,dx1,dx2,dp
    real(kind=h_fpp) :: dV,dE
    integer, dimension(3) :: idx1,idx2
    integer :: i,j,k

    !assume we always receive hg_nf accelerations, so have to downweight them
    accr=acc/(h%dx)**3

    if (gc.eq.0) then

       i=1+floor(x(1)/h%dx)
       j=1+floor(x(2)/h%dx)
       k=1+floor(x(3)/h%dx)

       dV=1
       dE=sum(h%fld(2:4,i,j,k)*accr)+0.5*sum(accr**2*h%fld(1,i,j,k))
       dp=accr*h%fld(1,i,j,k)

       h%fld(5,i,j,k)=h%fld(5,i,j,k)+dV*dE
       h%fld(2:4,i,j,k)=h%fld(2:4,i,j,k)+dV*dp

    else
       
       xx=x/h%dx-0.5

       idx1=1+floor(xx)
       idx2=1+idx1
       dx1=idx1-xx
       dx2=1.-dx1

       !1,1,1
       i=idx1(1);j=idx1(2);k=idx1(3)
       dV=dx1(1)*dx1(2)*dx1(3)
       dE=sum(h%fld(2:4,i,j,k)*accr)+0.5*sum(accr**2*h%fld(1,i,j,k))
       dp=accr*h%fld(1,i,j,k)
       h%fld(5,i,j,k)=h%fld(5,i,j,k)+dV*dE
       h%fld(2:4,i,j,k)=h%fld(2:4,i,j,k)+dV*dp

       !2,1,1
       i=idx2(1);j=idx1(2);k=idx1(3)
       dV=dx2(1)*dx1(2)*dx1(3)
       dE=sum(h%fld(2:4,i,j,k)*accr)+0.5*sum(accr**2*h%fld(1,i,j,k))
       dp=accr*h%fld(1,i,j,k)
       h%fld(5,i,j,k)=h%fld(5,i,j,k)+dV*dE
       h%fld(2:4,i,j,k)=h%fld(2:4,i,j,k)+dV*dp

       !1,2,1
       i=idx1(1);j=idx2(2);k=idx1(3)
       dV=dx1(1)*dx2(2)*dx1(3)
       dE=sum(h%fld(2:4,i,j,k)*accr)+0.5*sum(accr**2*h%fld(1,i,j,k))
       dp=accr*h%fld(1,i,j,k)
       h%fld(5,i,j,k)=h%fld(5,i,j,k)+dV*dE
       h%fld(2:4,i,j,k)=h%fld(2:4,i,j,k)+dV*dp

       !2,2,1
       i=idx2(1);j=idx2(2);k=idx1(3)
       dV=dx2(1)*dx2(2)*dx1(3)
       dE=sum(h%fld(2:4,i,j,k)*accr)+0.5*sum(accr**2*h%fld(1,i,j,k))
       dp=accr*h%fld(1,i,j,k)
       h%fld(5,i,j,k)=h%fld(5,i,j,k)+dV*dE
       h%fld(2:4,i,j,k)=h%fld(2:4,i,j,k)+dV*dp

       !1,1,2
       i=idx1(1);j=idx1(2);k=idx2(3)
       dV=dx1(1)*dx1(2)*dx2(3)
       dE=sum(h%fld(2:4,i,j,k)*accr)+0.5*sum(accr**2*h%fld(1,i,j,k))
       dp=accr*h%fld(1,i,j,k)
       h%fld(5,i,j,k)=h%fld(5,i,j,k)+dV*dE
       h%fld(2:4,i,j,k)=h%fld(2:4,i,j,k)+dV*dp

       !2,1,2
       i=idx2(1);j=idx1(2);k=idx2(3)
       dV=dx2(1)*dx1(2)*dx2(3)
       dE=sum(h%fld(2:4,i,j,k)*accr)+0.5*sum(accr**2*h%fld(1,i,j,k))
       dp=accr*h%fld(1,i,j,k)
       h%fld(5,i,j,k)=h%fld(5,i,j,k)+dV*dE
       h%fld(2:4,i,j,k)=h%fld(2:4,i,j,k)+dV*dp

       !1,2,2
       i=idx1(1);j=idx2(2);k=idx2(3)
       dV=dx1(1)*dx2(2)*dx2(3)
       dE=sum(h%fld(2:4,i,j,k)*accr)+0.5*sum(accr**2*h%fld(1,i,j,k))
       dp=accr*h%fld(1,i,j,k)
       h%fld(5,i,j,k)=h%fld(5,i,j,k)+dV*dE
       h%fld(2:4,i,j,k)=h%fld(2:4,i,j,k)+dV*dp
       
       !2,2,2
       i=idx2(1);j=idx2(2);k=idx2(3)
       dV=dx2(1)*dx2(2)*dx2(3)
       dE=sum(h%fld(2:4,i,j,k)*accr)+0.5*sum(accr**2*h%fld(1,i,j,k))
       dp=accr*h%fld(1,i,j,k)
       h%fld(5,i,j,k)=h%fld(5,i,j,k)+dV*dE
       h%fld(2:4,i,j,k)=h%fld(2:4,i,j,k)+dV*dp

    end if

  end subroutine hg_gravity

  subroutine hg_buffer(h)
    implicit none
    class(hydro) :: h
    integer hg_ex
    real(h_fpp),allocatable :: slicex(:,:,:,:)[:],slicey(:,:,:,:)[:],slicez(:,:,:,:)[:]
   
    hg_ex=h%nc+2*h%ncb
    allocate(slicex(5,hg_nfb,hg_ex,hg_ex)[*])
    allocate(slicey(5,hg_ex,hg_nfb,hg_ex)[*])
    allocate(slicez(5,hg_ex,hg_ex,hg_nfb)[*])

    !real(kind=h_fpp), codimension[*], save :: slicex(5,1:hg_nfb,1-hg_nfb:hg_nf+hg_nfb,1-hg_nfb:hg_nf+hg_nfb)
    !real(kind=h_fpp), codimension[*], save :: slicey(5,1-hg_nfb:hg_nf+hg_nfb,1:hg_nfb,1-hg_nfb:hg_nf+hg_nfb)
    !real(kind=h_fpp), codimension[*], save :: slicez(5,1-hg_nfb:hg_nf+hg_nfb,1-hg_nfb:hg_nf+hg_nfb,1:hg_nfb)
    
    !integer :: i

    !+x
    slicex=h%fld(:,h%nc-h%ncb+1:h%nc,:,:)
    h%fld(:,:0,:,:)=slicex(:,:,:,:)[image1d(inx,icy,icz)]
    !h%fld(1:5,1-h%ncb:0,1:h%nc,1:h%nc)=h%fld(1:5,h%nc-h%ncb+1:h%nc,1:h%nc,1:h%nc)
    !-x
    slicex=h%fld(:,1:h%ncb,:,:)
    h%fld(:,h%nc+1:,:,:)=slicex(:,:,:,:)[image1d(ipx,icy,icz)]
    !h%fld(1:5,h%nc+1:h%nc+h%ncb,1:h%nc,1:h%nc)=h%fld(1:5,1:h%ncb,1:h%nc,1:h%nc)
    sync all

    !+y
    slicey=h%fld(:,:,h%nc-h%ncb+1:h%nc,:)
    h%fld(:,:,:0,:)=slicey(:,:,:,:)[image1d(icx,iny,icz)]
    !h%fld(1:5,1-h%ncb:h%nc+h%ncb,1-h%ncb:0,1:h%nc)=h%fld(1:5,1-h%ncb:h%nc+h%ncb,h%nc-h%ncb+1:h%nc,1:h%nc)
    !-y
    slicey=h%fld(:,:,1:h%ncb,:)
    h%fld(:,:,h%nc+1:,:)=slicey(:,:,:,:)[image1d(icx,ipy,icz)]
    !h%fld(1:5,1-h%ncb:h%nc+h%ncb,h%nc+1:h%nc+h%ncb,1:h%nc)=h%fld(1:5,1-h%ncb:h%nc+h%ncb,1:h%ncb,1:h%nc)
    sync all

    !+z
    slicez=h%fld(:,:,:,h%nc-h%ncb+1:h%nc)
    h%fld(:,:,:,:0)=slicez(:,:,:,:)[image1d(icx,icy,inz)]
    !h%fld(1:5,1-h%ncb:h%nc+h%ncb,1-h%ncb:h%nc+h%ncb,1-h%ncb:0)=h%fld(1:5,1-h%ncb:h%nc+h%ncb,1-h%ncb:h%nc+h%ncb,h%nc-h%ncb+1:h%nc)
    !-z
    slicez=h%fld(:,:,:,1:h%ncb)
    h%fld(:,:,:,h%nc+1:)=slicez(:,:,:,:)[image1d(icx,icy,ipz)]
    !h%fld(1:5,1-h%ncb:h%nc+h%ncb,1-h%ncb:h%nc+h%ncb,h%nc+1:h%nc+h%ncb)=h%fld(1:5,1-h%ncb:h%nc+h%ncb,1-h%ncb:h%nc+h%ncb,1:h%ncb)
    sync all
    deallocate(slicex,slicey,slicez)
  end subroutine hg_buffer

  !Subroutine to set a fluid element given an input array that may be of different size
  !assumes size of r3d is divisible by h%nc
  subroutine hg_set_fld(h,r3d,d)
    implicit none
    class(hydro) :: h
    real, dimension(:,:,:), intent(in) :: r3d
    integer, intent(in) :: d
    integer :: n,i,j,k
    real(8) :: r

    n=size(r3d,dim=1)/h%nc
    if (n.ne.size(r3d,dim=2)/h%nc .or. n.ne.size(r3d,dim=3)/h%nc) error stop
    r=(1.d0*n)**3.

    do k=1,h%nc
    do j=1,h%nc
    do i=1,h%nc
       h%fld(d,i,j,k)=sum( r3d((i-1)*n+1:i*n,(j-1)*n+1:j*n,(k-1)*n+1:k*n) )/r
    end do
    end do
    end do

  end subroutine hg_set_fld

  !Subroutine to compute interpolated density at coordinate x
  !gc controls NGP vs CIC interpolation: gc=0 yields NGP, gc!=0 yields CIC
  function hg_density(h,x,gc) result(d)
    implicit none
    class(hydro) :: h
    real, dimension(3), intent(in) :: x
    integer, intent(in) :: gc
    integer, dimension(3) :: idx1,idx2
    real, dimension(3) :: xx,dx1,dx2
    real(kind=h_fpp) :: d

    if (gc.eq.0) then

       idx1=1+floor(x/h%dx)
       d=h%fld(1,idx1(1),idx1(2),idx1(3))

    else

       xx=x/h%dx-0.5

       idx1=1+floor(xx)
       idx2=1+idx1
       dx1=idx1-xx
       dx2=1.-dx1

       d     = h%fld(1,idx1(1),idx1(2),idx1(3))*dx1(1)*dx1(2)*dx1(3) &
            &+ h%fld(1,idx2(1),idx1(2),idx1(3))*dx2(1)*dx1(2)*dx1(3) &
            &+ h%fld(1,idx1(1),idx2(2),idx1(3))*dx1(1)*dx2(2)*dx1(3) &
            &+ h%fld(1,idx2(1),idx2(2),idx1(3))*dx2(1)*dx2(2)*dx1(3) &
            &+ h%fld(1,idx1(1),idx1(2),idx2(3))*dx1(1)*dx1(2)*dx2(3) &
            &+ h%fld(1,idx2(1),idx1(2),idx2(3))*dx2(1)*dx1(2)*dx2(3) &
            &+ h%fld(1,idx1(1),idx2(2),idx2(3))*dx1(1)*dx2(2)*dx2(3) &
            &+ h%fld(1,idx2(1),idx2(2),idx2(3))*dx2(1)*dx2(2)*dx2(3) 

    end if

  end function hg_density

  !Generally useful subroutines
  subroutine hg_write(str,vl,val)
    implicit none
    character(len=*), intent(in) :: str
    integer, intent(in) :: vl
    real(kind=h_fpp), intent(in), optional :: val
    if (present(val)) then
       if (hg_verb.ge.vl) write(*,*) str,val
    else
       if (hg_verb.ge.vl) write(*,*) str
    end if
  end subroutine hg_write

  subroutine hg_global_min(t)
    implicit none
    real, intent(inout) :: t
    real, codimension[*], save :: tl
    integer :: n
    tl=t
    do n=1,nn**3
       t=min(t,tl[n])
    end do
    sync all
  end subroutine hg_global_min

  subroutine hg_hydro_properties(h)
    implicit none
    class(hydro) :: h
    call hg_write('      >Hydro properties after evolution',1)
    call hg_array_properties(h%fld(1,1:h%nc,1:h%nc,1:h%nc),'      density',hg_verb.ge.1)
    call hg_array_properties(h%fld(2,1:h%nc,1:h%nc,1:h%nc)/h%fld(1,1:h%nc,1:h%nc,1:h%nc),'      x-velocity',hg_verb.ge.2)
    call hg_array_properties(h%fld(3,1:h%nc,1:h%nc,1:h%nc)/h%fld(1,1:h%nc,1:h%nc,1:h%nc),'      y-velocity',hg_verb.ge.2)
    call hg_array_properties(h%fld(4,1:h%nc,1:h%nc,1:h%nc)/h%fld(1,1:h%nc,1:h%nc,1:h%nc),'      z-velocity',hg_verb.ge.2)
    call hg_array_properties(h%fld(5,1:h%nc,1:h%nc,1:h%nc),'      energy',hg_verb.ge.2)
  end subroutine hg_hydro_properties

  subroutine hg_array_properties(arr,str,ver)
    implicit none
    real(kind=h_fpp), dimension(:,:,:), intent(in) :: arr
    character(len=*), intent(in) :: str
    logical, intent(in) :: ver
    
    real(8), codimension[*], save :: s8,s8g,max8,min8
    integer :: n

    s8g=0.
    s8=sum(arr*1.d0)
    max8=maxval(arr*1.d0)
    min8=minval(arr*1.d0)
    do n=1,nn**3
       s8g=s8g+s8[n]
       max8=max(max8,max8[n])
       min8=min(min8,min8[n])
    end do
    sync all

    if (ver) write(*,*) str//' avg, max, min',real(s8g/size(arr),kind=h_fpp),real(max8,kind=h_fpp),real(min8,kind=h_fpp)

  end subroutine hg_array_properties

  subroutine hg_checkpoint(h,fn)
    implicit none
    class(hydro) :: h
    character(len=*), intent(in) :: fn
    integer :: stat

    call hg_write('Checkpointing to file: '//fn,0)
    
    !Open file and write out
    open(unit=11,file=trim(adjustl(fn)),access='stream',iostat=stat,status='replace')
    if (stat.ne.0) then
       write(*,*) 'ERROR in checkpoint: could not open outfile file: '//trim(adjustl(fn))
       error stop
    end if
    write(11) h%nc,h%ncb,h%g,h%cs2,h%dx,h%dv,h%isothermal
    write(11) h%fld
    close(11)

  end subroutine hg_checkpoint

  subroutine hg_restart(h,fn)
    implicit none
    class(hydro) :: h
    character(len=*), intent(in) :: fn
    integer :: stat

    call hg_write('Checkpointing to file: '//fn,0)

    !Open file and write out
    open(unit=11,file=trim(adjustl(fn)),access='stream',iostat=stat,status='old')
    if (stat.ne.0) then
       write(*,*) 'ERROR in restart: could not open input file: '//trim(adjustl(fn))
       error stop
    end if

    if (allocated(h%fld)) deallocate(h%fld)
    read(11) h%nc,h%ncb,h%g,h%cs2,h%dx,h%dv,h%isothermal
    allocate(h%fld(5,1-h%ncb:h%nc+h%ncb,1-h%ncb:h%nc+h%ncb,1-h%ncb:h%nc+h%ncb))
    read(11) h%fld
    close(11)
    
  end subroutine hg_restart

end module HydroGR

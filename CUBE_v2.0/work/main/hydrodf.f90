module HydroDF
  use hydrogr
  implicit none

  !Constants
  !!Beta = m/T with m=1eV
  real(kind=h_fpp), parameter :: dfld_beta_eV = 5966.34
  !!Speed of light
  real(kind=h_fpp), parameter :: dfld_sol = 2.998e5
  !!When to start evolution of components
  real(kind=h_fpp), parameter :: dfld_ai_kms = (8.*pi**2./3.)/(box**2.*(h0*100.)**2*omega_m)/hg_dv**2/2.
  real(kind=h_fpp), parameter :: dfld_ai_max = 0.05
  real(kind=h_fpp), parameter :: dfld_ai_soft_g = 1.
  real(kind=h_fpp), parameter :: dfld_ai_soft_h = 0.9*dfld_ai_soft_g !Start hydro slightly before gravity

  !Number of isothermal fluids in dispersive sum
  integer, parameter :: dfld_n_hydro = 3
  !Gauss-Laguerre coordinates and weights
  real(kind=h_fpp), dimension(dfld_n_hydro), parameter :: dfld_GL_v = (/0.41577456,2.29428036,6.28994508/)
  real(kind=h_fpp), dimension(dfld_n_hydro), parameter :: dfld_GL_w = (/0.71109301,0.27851773,0.01038926/)
  !Cell sizes
  integer, dimension(dfld_n_hydro), parameter :: dfld_nc = (/ hg_nf/2,hg_nf/2,hg_nf/4 /)

  !Verbosity level (-1=nothing, 0=a little, 1=a lot)
  integer :: dfld_verbosity=1

  !Dispersive fluid = sum of fluids
  type dfld
     !Neutrino mass parameter: beta=m/T
     real(kind=h_fpp) :: beta
     !Array of fluids
     type(hydro), dimension(dfld_n_hydro) :: n_hydro
     !Weights
     real(kind=h_fpp), dimension(dfld_n_hydro) :: dwgts
   contains
     procedure :: setup => dfld_setup
     procedure :: evolve => dfld_evolve
     procedure :: gravity => dfld_gravity
     procedure :: density => dfld_density
     procedure :: density_field => dfld_density_field
     procedure :: checkpoint => dfld_checkpoint
     procedure :: restart => dfld_restart
  end type dfld

contains

  subroutine dfld_setup(d,g,m,iso)
    implicit none
    class(dfld) :: d
    real(kind=h_fpp), intent(in) :: m,g
    logical, intent(in), optional :: iso
    integer :: n
    real(kind=h_fpp) :: v

    if (dfld_verbosity.ge.0) write(*,*) 'Setting up dispersive sum with m/eV = ',m
    !Compute beta=m/T
    d%beta = dfld_beta_eV*m

    if (present(iso)) then
       write(*,*) 'Isothermal: ',iso
    end if

    !Compute density weights, setup isothermal fluids
    do n=1,dfld_n_hydro
       v = dfld_GL_v(n)/d%beta*dfld_sol !km/s
       if (present(iso)) then
          call d%n_hydro(n)%setup(dfld_nc(n),g,v,iso)
       else
          call d%n_hydro(n)%setup(dfld_nc(n),g,v)
       end if
       d%dwgts(n) = dfld_GL_w(n)*f0(dfld_GL_v(n))*exp(dfld_GL_v(n))
    end do 

    if (dfld_verbosity.ge.0) write(*,*) '> sum of density weights',sum(d%dwgts),1.80305
    !Normalize weights
    d%dwgts = d%dwgts/sum(d%dwgts)

  end subroutine dfld_setup

  subroutine dfld_evolve(d,dt,dtg,dir,a)
    implicit none
    class(dfld) :: d
    real, intent(in) :: dt,a
    real, intent(inout) :: dtg
    integer, intent(in) :: dir
    integer :: n
    real :: ai

    do n=1,dfld_n_hydro
       ai=min(dfld_ai_kms*d%n_hydro(n)%cs2,dfld_ai_max)
       if (a.ge.ai) then
          call d%n_hydro(n)%evolve(dt,dtg,dir)
       else if (a.ge.dfld_ai_soft_h*ai) then
          !Start hydro slightly earlier than gravity to include it in the timestep
          if (dfld_verbosity.ge.0) write(*,*) 'Evolving dfld component with softened force'
          call d%n_hydro(n)%evolve(dt,dtg,dir)
       else
          if (dfld_verbosity.ge.0) write(*,*) 'Not evolving dfld component with cs until zi',sqrt(d%n_hydro(n)%cs2)/hg_dv,1./ai-1.
       end if
    end do
    
  end subroutine dfld_evolve

  subroutine dfld_gravity(d,x,acc,a)
    implicit none
    class(dfld) :: d
    real, intent(in) :: a
    real, dimension(3) :: x
    real(8), dimension(3) :: acc
    integer :: n
    real :: ai    

    do n=1,dfld_n_hydro
       ai=min(dfld_ai_kms*d%n_hydro(n)%cs2,dfld_ai_max)
       if (a.ge.ai) then
          call d%n_hydro(n)%gravity(x,acc)
       else if (a.ge.dfld_ai_soft_g*ai) then
          !call d%n_hydro(n)%gravity(x,acc*(a/ai)**4.)
          call d%n_hydro(n)%gravity(x,acc)
       end if
    end do

  end subroutine dfld_gravity

  subroutine dfld_checkpoint(d,fn)
    implicit none
    class(dfld) :: d
    character(len=*), intent(in) :: fn
    integer :: n,stat

    call hg_write('Checkpointing to file: '//fn,0)

    !Open file and write out
    open(unit=11,file=trim(adjustl(fn)),access='stream',iostat=stat,status='replace')
    if (stat.ne.0) then
       write(*,*) 'ERROR in checkpoint: could not open outfile file: '//trim(adjustl(fn))
       error stop
    end if
    write(11) d%beta/dfld_beta_eV,d%n_hydro(1)%g,d%n_hydro(1)%isothermal,d%dwgts
    
    do n=1,dfld_n_hydro
       write(11) d%n_hydro(n)%fld
    end do

    close(11)

  end subroutine dfld_checkpoint

  subroutine dfld_restart(d,fn)
    implicit none
    class(dfld) :: d
    character(len=*), intent(in) :: fn
    integer :: n,stat
    real(kind=h_fpp) :: m,g
    logical :: iso

    call hg_write('Restarting from file: '//fn,0)

    !Open file and write out
    open(unit=11,file=trim(adjustl(fn)),access='stream',iostat=stat,status='old')
    if (stat.ne.0) then
       write(*,*) 'ERROR in checkpoint: could not open inout file: '//trim(adjustl(fn))
       error stop
    end if
   
    read(11) m,g,iso,d%dwgts
    call d%setup(g,m,iso) 
    do n=1,dfld_n_hydro
       read(11) d%n_hydro(n)%fld
    end do

    close(11)

  end subroutine dfld_restart

  function dfld_density_field(d) result(density)
    implicit none
    class(dfld) :: d
    real(kind=h_fpp), dimension(hg_nf/2,hg_nf/2,hg_nf/2) :: density
    integer :: n,i,j,k
    real, dimension(3) :: x

    density=0.
    do n=1,dfld_n_hydro
       do k=1,hg_nf/2
          do j=1,hg_nf/2
             do i=1,hg_nf/2
                x=2.*(/i,j,k/)-1.
                density(i,j,k)=density(i,j,k)+d%dwgts(n)*d%n_hydro(n)%density(x)
             end do
          end do
       end do
    end do

  end function dfld_density_field

  function dfld_density(d,x) result(density)
    implicit none
    class(dfld) :: d
    real, dimension(3), intent(in) :: x
    real(kind=h_fpp) :: density
    integer :: n
    
    density=0.
    do n=1,dfld_n_hydro
       density=density+d%dwgts(n)*d%n_hydro(n)%density(x)
    end do

  end function dfld_density

  !Background distribution function
  function f0(v)
    implicit none
    real(kind=h_fpp), intent(in) :: v
    real(kind=h_fpp) :: f0
    f0=v**2/(exp(v)+1.)
  end function f0

end module HydroDF

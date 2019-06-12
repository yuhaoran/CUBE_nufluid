!1D hydrodynamics code
!Based on homework set from Andrew MacFadyen

!#define HIGH_RES_T
#define HIGH_RES_X
!#define P_SEMILINEAR
!#define q_EQN_OF_STATE
!#define q_SEMILINEAR
module hydro1d
  use, intrinsic :: ISO_FORTRAN_ENV
  implicit none

  !Numerical parameters
  !!Floating point precision
  integer, parameter :: h_fpp = REAL32
  !!CFL condition
  real(kind=h_fpp), parameter :: h_cfl = 0.71!85
  !!Flux limiter
  real(kind=h_fpp) :: h_lim = 1.5

#ifdef HIGH_RES_T
  logical, parameter :: h_hrt = .true.
#else
  logical, parameter :: h_hrt = .false.
#endif

contains

  !Fluid variables
  !DIR$ ATTRIBUTES INLINE :: h_density
  pure function h_density(h) result(d)
    implicit none
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    real(kind=h_fpp), dimension(size(h,2)) :: d
    d=h(1,:)
  end function h_density

  !DIR$ ATTRIBUTES INLINE :: h_momentum
  pure function h_momentum(h,d) result(m)
    implicit none
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    integer, intent(in) :: d
    real(kind=h_fpp), dimension(size(h,2)) :: m
    m=h(d,:)
  end function h_momentum

  !DIR$ ATTRIBUTES INLINE :: h_energy
  pure function h_energy(h) result(e)
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    real(kind=h_fpp), dimension(size(h,2)) :: e
    e=h(5,:)
  end function h_energy

  !Derived fluid variables
  !DIR$ ATTRIBUTES INLINE :: h_velocity
  pure function h_velocity(h,d) result(v)
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    integer, intent(in) :: d
    real(kind=h_fpp), dimension(size(h,2)) :: v
    v=h_momentum(h,d)/h_density(h)
  end function h_velocity

  !DIR$ ATTRIBUTES INLINE :: h_kinetic
  pure function h_kinetic(h,cs2) result(k)
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    real(kind=h_fpp), intent(in), optional :: cs2
    real(kind=h_fpp), dimension(size(h,2)) :: k
    k=0.5*sum(h(2:4,:)**2,1)/h(1,:)
  end function h_kinetic

  !DIR$ ATTRIBUTES INLINE :: h_pressure
  pure function h_pressure(h,g,cs2) result(p)
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    real(kind=h_fpp), intent(in) :: g
    real(kind=h_fpp), intent(in), optional :: cs2
    real(kind=h_fpp), dimension(size(h,2)) :: p 
    p=max(0.,(h_energy(h)-h_kinetic(h))*(g-1.))
#ifdef P_SEMILINEAR
    if (present(cs2)) p=cs2*h_density(h)
#endif
  end function h_pressure

  !DIR$ ATTRIBUTES INLINE :: h_soundspd
  pure function h_soundspd(h,g,cs2) result(s)
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    real(kind=h_fpp), intent(in) :: g
    real(kind=h_fpp), intent(in), optional :: cs2
    real(kind=h_fpp), dimension(size(h,2)) :: s
    s=sqrt(g*h_pressure(h,g)/h_density(h))
#ifdef P_SEMILINEAR
    s=sqrt(cs2)
#endif
#ifdef q_EQN_OF_STATE
    s=sqrt(h_pressure(h,g,cs2)/h_density(h))
#endif
  end function h_soundspd

  !Subroutine to compute timestep
  pure subroutine h_timestep(h,g,dt,dx,cs2)
    implicit none
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    real(kind=h_fpp), intent(in) :: g,dx
    real(kind=h_fpp), intent(out) :: dt
    real(kind=h_fpp), intent(in), optional :: cs2
    real(kind=h_fpp) :: v,s

    v=max(0.,maxval(abs(h_velocity(h,2))),maxval(abs(h_velocity(h,3))),maxval(abs(h_velocity(h,4))))
    s=maxval(h_soundspd(h,g,cs2))
    dt=h_cfl*dx/(v+s)

  end subroutine h_timestep

  !Subroutine to update fluid array
  pure subroutine h_evolve(h,g,dt,dx,cs2)
    real(kind=h_fpp), dimension(:,:), intent(inout) :: h
    real(kind=h_fpp), intent(in) :: g,dt,dx
    real(kind=h_fpp), intent(in), optional :: cs2
    real(kind=h_fpp), dimension(size(h,1),size(h,2)) :: h12

#ifdef HIGH_RES_T
    h12=h+dt*h_riemann(h,g,cs2)/dx
    h12=(3./4.)*h+(1./4.)*h12+(1./4.)*dt*h_riemann(h,g,cs2)/dx
    h=(1./3.)*h+(2./3.)*h12+(2./3.)*dt*h_riemann(h12,g,cs2)/dx
#else
    h=h+dt*h_riemann(h,g,cs2)/dx
#endif

  end subroutine h_evolve

  !Function to compute fluxes
  pure function h_riemann(h,g,cs2) result(f)
    implicit none
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    real(kind=h_fpp), intent(in) :: g
    real(kind=h_fpp), intent(in), optional :: cs2
    real(kind=h_fpp), dimension(size(h,1),size(h,2)) :: hl,hr,fl,fr,f
    real(kind=h_fpp), dimension(size(h,2)) :: ap,am,al,ar
    integer :: d

    !Left biased density, velocity, pressure + state and flux
    hl(1,:)=h_left_biased(h(1,:))
    hl(2,:)=h_density(hl)*h_left_biased(h_velocity(h,2))
    hl(3,:)=h_density(hl)*h_left_biased(h_velocity(h,3))
    hl(4,:)=h_density(hl)*h_left_biased(h_velocity(h,4))
    hl(5,:)=h_left_biased(h_pressure(h,g,cs2))/(g-1.)+h_kinetic(hl,cs2)
    fl=h_flux(hl,g,cs2)

    !Right biased density, velocity, pressure + state and flux
    hr(1,:)=h_right_biased(h(1,:))
    hr(2,:)=h_density(hr)*h_right_biased(h_velocity(h,2))
    hr(3,:)=h_density(hr)*h_right_biased(h_velocity(h,3))
    hr(4,:)=h_density(hr)*h_right_biased(h_velocity(h,4))
    hr(5,:)=h_right_biased(h_pressure(h,g,cs2))/(g-1.)+h_kinetic(hr,cs2)
    fr=h_flux(hr,g,cs2)

    !Eigenvalues
    ap=max(0.,h_velocity(hl,2)+h_soundspd(hl,g,cs2),h_velocity(hr,2)+h_soundspd(hr,g,cs2))
    am=max(0.,-h_velocity(hl,2)+h_soundspd(hl,g,cs2),-h_velocity(hr,2)+h_soundspd(hr,g,cs2))

    do d=1,5
       f(d,:)=(ap*fl(d,:)+am*fr(d,:)-ap*am*(hr(d,:)-hl(d,:)))/(ap+am)
    end do

    f=-f+cshift(f,-1,2)

  end function h_riemann

  !DIR$ ATTRIBUTES INLINE :: h_flux
  pure function h_flux(h,g,cs2) result(f)
    real(kind=h_fpp), dimension(:,:), intent(in) :: h
    real(kind=h_fpp), intent(in) :: g
    real(kind=h_fpp), intent(in), optional :: cs2
    real(kind=h_fpp), dimension(size(h,1),size(h,2)) :: f

    f(1,:)=h_momentum(h,2)
    f(2,:)=h_momentum(h,2)*h_velocity(h,2)+h_pressure(h,g,cs2)
    f(3,:)=h_momentum(h,2)*h_velocity(h,3)
    f(4,:)=h_momentum(h,2)*h_velocity(h,4)
    f(5,:)=(h_energy(h)+h_pressure(h,g,cs2))*h_velocity(h,2)
#ifdef q_EQN_OF_STATE
    f(5,:)=h_energy(h)*h_velocity(h,2)
#ifdef q_SEMILINEAR
    if (present(cs2)) f(5,:)=(h_energy(h)+h_pressure(h,g,cs2)-h_density(h)*cs2)*h_velocity(h,2)
#endif
#endif
  end function h_flux

  !DIR$ ATTRIBUTES INLINE :: h_left_biased
  pure function h_left_biased(u) result(l)
    implicit none
    real(kind=h_fpp), dimension(:), intent(in) :: u
    real(kind=h_fpp), dimension(size(u)) :: l

#ifdef HIGH_RES_X
    l=u+0.5*h_minmod( h_lim*(u-cshift(u,-1)), &
         &          0.5*(cshift(u,1)-cshift(u,-1)), &
         &          h_lim*(cshift(u,1)-u) )
#else
    l=u !First order    
#endif

  end function h_left_biased
  
  !DIR$ ATTRIBUTES INLINE :: h_right_biased
  pure function h_right_biased(u) result(r)
    implicit none
    real(kind=h_fpp), dimension(:), intent(in) :: u
    real(kind=h_fpp), dimension(size(u)) :: r

#ifdef HIGH_RES_X
    r=cshift(u,1)-0.5*h_minmod( h_lim*(cshift(u,1)-u), &
         &                    0.5*(cshift(u,2)-u), &
         &                    h_lim*(cshift(u,2)-cshift(u,1)) )
#else
    r=cshift(u,1) !First order
#endif
  end function h_right_biased

  !DIR$ ATTRIBUTES INLINE :: h_minmod
  elemental function h_minmod(a,b,c) result(m)
    implicit none
    real(kind=h_fpp), intent(in) :: a,b,c
    real(kind=h_fpp) :: m
    real(kind=h_fpp), parameter :: one = real(1.,kind=h_fpp)

    m=0.25*abs(sign(one,a)+sign(one,b))*(sign(one,a)+sign(one,c))&
         &*min( abs(a),abs(b),abs(c) )

  end function h_minmod

end module hydro1d

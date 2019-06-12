!3D hydro code
!Based on Trac and Pen 2003 with the relaxing algorithm 
!replaced by the hydro1d approximate Riemann solver 
module hydro3d
  use hydro1d
  implicit none

contains

  subroutine h3_timestep(h,g,dt,dx,cs2)
    implicit none
    real(kind=h_fpp), dimension(:,:,:,:), intent(in) :: h
    real(kind=h_fpp), intent(in) :: g,dx
    real(kind=h_fpp), intent(out) :: dt
    real(kind=h_fpp), intent(in), optional :: cs2
    real(kind=h_fpp) :: dtn
    integer :: j,k

    dt=1e9
    !$omp parallel do default(none) shared(h,g,dx,cs2) private(j,k,dtn) reduction(min:dt)
    do k=1,size(h,4)
       do j=1,size(h,3)
          call h_timestep(h(:,:,j,k),g,dtn,dx,cs2)
          dt=min(dt,dtn)
       end do
    end do
    !$omp end parallel do
    !dt=2.*dt
    if (dt.eq.1e9) dt=1.

  end subroutine h3_timestep

  subroutine h3_evolve(h,g,dt,dx,dir,cs2)
    implicit none
    real(kind=h_fpp), dimension(:,:,:,:), intent(inout) :: h
    real(kind=h_fpp), intent(in) :: g,dt,dx
    integer, intent(in) :: dir
    real(kind=h_fpp), intent(in), optional :: cs2

    if (dir.gt.0) then
       call h3_evolve_xyz(h,g,dt,dx,1,cs2)
       call h3_evolve_xyz(h,g,dt,dx,2,cs2)
       call h3_evolve_xyz(h,g,dt,dx,3,cs2)
    else
       call h3_evolve_xyz(h,g,dt,dx,3,cs2)
       call h3_evolve_xyz(h,g,dt,dx,2,cs2)
       call h3_evolve_xyz(h,g,dt,dx,1,cs2)
    end if

  end subroutine h3_evolve

  subroutine h3_evolve_xyz(h,g,dt,dx,xyz,cs2)
    real(kind=h_fpp), dimension(:,:,:,:), intent(inout) :: h
    real(kind=h_fpp), intent(in) :: g,dt,dx
    integer, intent(in) :: xyz
    real(kind=h_fpp), intent(in), optional :: cs2

    real(kind=h_fpp), dimension(size(h,1),size(h,2)) :: h1dx
    real(kind=h_fpp), dimension(size(h,1),size(h,3)) :: h1dy
    real(kind=h_fpp), dimension(size(h,1),size(h,4)) :: h1dz
    integer :: j,k

    if (xyz.eq.1) then
       if (size(h,2).le.7) return
       !$omp parallel do default(none) shared(h,g,dt,dx,cs2) private(j,k,h1dx)
       do k=1,size(h,4)
          do j=1,size(h,3)
             h1dx(:,:)=h(:,:,j,k)
             call h_evolve(h1dx,g,dt,dx,cs2)
             h(:,:,j,k)=h1dx(:,:)
          end do
       end do
       !$omp end parallel do
    else if (xyz.eq.2) then
       if (size(h,3).le.7) return
       !$omp parallel do default(none) shared(h,g,dt,dx,cs2) private(j,k,h1dy)
       do k=1,size(h,4)
          do j=1,size(h,2)
             h1dy(:,:)=h(:,j,:,k)
             h1dy(2,:)=h(3,j,:,k)
             h1dy(3,:)=h(2,j,:,k)
             call h_evolve(h1dy,g,dt,dx,cs2)
             h(:,j,:,k)=h1dy(:,:)
             h(3,j,:,k)=h1dy(2,:)
             h(2,j,:,k)=h1dy(3,:)
          end do
       end do
       !$omp end parallel do
    else if (xyz.eq.3) then
       if (size(h,4).le.7) return
       !$omp parallel do default(none) shared(h,g,dt,dx,cs2) private(j,k,h1dz)
       do k=1,size(h,3)
          do j=1,size(h,2)
             h1dz(:,:)=h(:,j,k,:)
             h1dz(2,:)=h(4,j,k,:)
             h1dz(4,:)=h(2,j,k,:)
             call h_evolve(h1dz,g,dt,dx,cs2)
             h(:,j,k,:)=h1dz(:,:)
             h(4,j,k,:)=h1dz(2,:)
             h(2,j,k,:)=h1dz(4,:)
          end do
       end do
       !$omp end parallel do
    else
       error stop 4
    end if
  
  end subroutine h3_evolve_xyz

end module hydro3d

!#define record_Fc
subroutine particle_mesh
  use omp_lib
  use variables
  use cubefft
  use pencil_fft
  use neutrinos
  use hydrogr
  implicit none
  save

#ifdef record_Fc
  character(10) :: str_a
#endif

  ! force settings
  !logical,parameter :: fine_force=.true.
  !logical,parameter :: coarse_force=.true.
  !logical,parameter :: pp_force=.false.
  !logical,parameter :: ext_pp_force=.false.

integer(4) t01,t02,t0rate

  integer(4),parameter :: nlayer=3 ! thread save for CIC integpolation
  integer(4) ilayer
  integer(8) idx1(3),idx2(3)
  real tempx(3),dx1(3),dx2(3)
  real r3t(-1:nt+2,-1:nt+2,-1:nt+2) ! coarse density on tile, with buffer=2
  real(8) testrhof, testrhoc
  if (head) then
    print*, ''
    print*, 'particle mesh'
    print*, '  mass_p cdm/nu =',sim%mass_p_cdm,sim%mass_p_nu
    print*, '  a_mid =',a_mid
    print*, '  dt =',dt
    call system_clock(tt1,t_rate)
  endif

  vmax=0; vmax_nu=0
  f2_max_fine(1:nnt,1:nnt,1:nnt)=0
  f2_max_coarse=0

  if (head) print*, '  pm fine over',nnt**3,'tiles'
  call system_clock(t1,t_rate)
  testrhof=0; testrhoc=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    !if (head) print*,'    tile',int(itx,1),int(ity,1),int(itz,1)
    ! fine_cic_mass ------------------------------------------------------------
    rho_f=0
    crho_f=0

    if (neutrino_flag) then
       if (1./a_mid-1. .gt. z_i_nu) then
          !add homogeneous
          rho_f=rho_f+sum(f_neu)
       else
          !add perturbation
          do k=1,nfe
          do j=1,nfe
          do i=1,nfe
             tempx=((/itx,ity,itz/)-1)*nt*ncell+(/i,j,k/)-nfb
             do nu=1,Nneu
                rho_f(i,j,k)=rho_f(i,j,k)+neu(nu)%density(tempx)*f_neu(nu)
             end do
          end do
          end do
          end do
       end if
    end if

    do ilayer=0,nlayer-1
      !$omp paralleldo default(shared) schedule(static,2)&
      !$omp& private(k,j,i,np,nzero,l,ip,tempx,idx1,idx2,dx1,dx2)
      do k=2-ncb+ilayer,nt+ncb-1,nlayer
      !do k=2-ncb,nt+ncb-1
      do j=2-ncb,nt+ncb-1
      do i=2-ncb,nt+ncb-1
        np=rhoc(i,j,k,itx,ity,itz)
        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
        do l=1,np ! loop over cdm particles
          ip=nzero+l
          tempx=ncell*((/i,j,k/)-1)+ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution !-0.5
          idx1 = floor(tempx) + 1
          idx2 = idx1 + 1
          dx1 = idx1 - tempx
          dx2 = 1 - dx1
          idx1=idx1+nfb
          idx2=idx2+nfb
          rho_f(idx1(1),idx1(2),idx1(3))=rho_f(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*sim%mass_p_cdm
          rho_f(idx2(1),idx1(2),idx1(3))=rho_f(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*sim%mass_p_cdm
          rho_f(idx1(1),idx2(2),idx1(3))=rho_f(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*sim%mass_p_cdm
          rho_f(idx1(1),idx1(2),idx2(3))=rho_f(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*sim%mass_p_cdm
          rho_f(idx1(1),idx2(2),idx2(3))=rho_f(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*sim%mass_p_cdm
          rho_f(idx2(1),idx1(2),idx2(3))=rho_f(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*sim%mass_p_cdm
          rho_f(idx2(1),idx2(2),idx1(3))=rho_f(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*sim%mass_p_cdm
          rho_f(idx2(1),idx2(2),idx2(3))=rho_f(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*sim%mass_p_cdm
        enddo
      enddo
      enddo
      enddo
      !$omp endparalleldo
    enddo ! ilayer
    testrhof=testrhof+sum(rho_f(nfb+1:nft+nfb,nfb+1:nft+nfb,nfb+1:nft+nfb)*1d0)
    ! fine force ---------------------------------------------------------------
    !if (head) print*,'      fine_fft'
    !call system_clock(t01,t0rate)
    call sfftw_execute(plan_fft_fine)
    !call system_clock(t02,t0rate)
    !print*, '  elapsed time =',real(t02-t01)/t0rate,'secs'; stop
    crho_f(:,:,:)=rho_f(:,:,:) ! back up
    do i_dim=1,3
      !if (head) print*,'      fine_ifft dim',int(i_dim,1)
      !$omp workshare
      rho_f(::2,:,:)=-crho_f(2::2,:,:)*kern_f(:,:,:,i_dim)
      rho_f(2::2,:,:)=crho_f(::2,:,:)*kern_f(:,:,:,i_dim)
      !$omp endworkshare
      call sfftw_execute(plan_ifft_fine)
      rho_f=rho_f/real(nfe)/real(nfe)/real(nfe)
      force_f(i_dim,:,:,:)=rho_f(nfb:nfe-nfb+1,nfb:nfe-nfb+1,nfb:nfe-nfb+1)
    enddo
    f2_max_fine(itx,ity,itz)=maxval(sum(force_f(:,:,:,:)**2,1))
    ! fine velocity ------------------------------------------------------------
    !if (head) print*,'      fine velocity'
    !$omp paralleldo default(shared) schedule(static,2)&
    !$omp& private(k,j,i,nzero,np,l,ip,tempx,idx1,idx2,dx1,dx2,vreal)
    do k=1,nt
    do j=1,nt
    do i=1,nt ! loop over coarse cell
      np=rhoc(i,j,k,itx,ity,itz)
      nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
      do l=1,np ! loop over cdm particles
        ip=nzero+l
        tempx=ncell*((/i,j,k/)-1)+ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution !-0.5
        idx1 = floor(tempx) + 1
        idx2 = idx1 + 1
        dx1 = idx1 - tempx
        dx2 = 1 - dx1
        idx1=idx1+nfb
        idx2=idx2+nfb
        vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
        vreal=vreal+force_f(:,idx1(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx1(3)
        vreal=vreal+force_f(:,idx2(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx1(3)
        vreal=vreal+force_f(:,idx1(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx1(3)
        vreal=vreal+force_f(:,idx1(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx2(3)
        vreal=vreal+force_f(:,idx1(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx2(3)
        vreal=vreal+force_f(:,idx2(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx2(3)
        vreal=vreal+force_f(:,idx2(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx1(3)
        vreal=vreal+force_f(:,idx2(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx2(3)
        vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi_new*vrel_boost)*vreal)/pi,kind=izipv)
      enddo

      if (neutrino_flag) then
         do kk=1,ncell
         do jj=1,ncell
         do ii=1,ncell
            tempx=ncell*((/i,j,k/)-1) + ((/ii,jj,kk/))
            idx1(:)=floor(tempx(:))+1
            vreal=force_f(:,idx1(1),idx1(2),idx1(3))*a_mid*dt/6/pi
            tempx=((/itx,ity,itz/)-1)*nt*ncell+((/i,j,k/)-1)*ncell+((/ii,jj,kk/))
            do nu=1,Nneu
               call neu(nu)%gravity(tempx,vreal,a_mid)
            end do
         end do
         end do
         end do
      end if

    enddo
    enddo
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo ! itz
  sigma_vi=sigma_vi_new
  if (neutrino_flag) sigma_vi_nu=sigma_vi_new_nu
  call system_clock(t2,t_rate)
  print*, '    elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all
  !-----------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (head) print*, '  pm coarse'
  ! coarse_cic_mass ------------------------------------------------------------
  if (head) print*, '    coarse cic mass'
  call system_clock(t1,t_rate)
  r3=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt ! loop over tile
    r3t=0
    
    if (neutrino_flag) then
       if (1./a_mid-1..gt.z_i_nu) then
          r3t=r3t+sum(f_neu)
       else
          do k=1,nt
          do j=1,nt
          do i=1,nt
             do kk=1,ncell
             do jj=1,ncell
             do ii=1,ncell
                tempx=((/itx,ity,itz/)-1)*nt*ncell+((/i,j,k/)-1)*ncell+((/ii,jj,kk/))
                do nu=1,Nneu
                   r3t(i,j,k)=r3t(i,j,k)+neu(nu)%density(tempx)*f_neu(nu)/ncell**3
                end do
             end do
             end do
             end do
          end do
          end do
          end do
       end if
    end if

    do ilayer=0,nlayer-1
      !$omp paralleldo default(shared) schedule(static,2)&
      !$omp& private(k,j,i,np,nzero,l,ip,tempx,idx1,idx2,dx1,dx2)
      do k=0+ilayer,nt+1,nlayer
      !do k=0,nt+1
      do j=0,nt+1
      do i=0,nt+1
        np=rhoc(i,j,k,itx,ity,itz)
        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
        do l=1,np ! loop over particle
          ip=nzero+l
          tempx=((/i,j,k/)-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution-0.5
          idx1(:)=floor(tempx(:))+1
          idx2(:)=idx1(:)+1
          dx1(:)=idx1(:)-tempx(:) ! CIC contribution to idx1
          dx2(:)=1-dx1(:) ! CIC contribution to idx2
          r3t(idx1(1),idx1(2),idx1(3))=r3t(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*sim%mass_p_cdm
          r3t(idx2(1),idx1(2),idx1(3))=r3t(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*sim%mass_p_cdm
          r3t(idx1(1),idx2(2),idx1(3))=r3t(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*sim%mass_p_cdm
          r3t(idx1(1),idx1(2),idx2(3))=r3t(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*sim%mass_p_cdm
          r3t(idx1(1),idx2(2),idx2(3))=r3t(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*sim%mass_p_cdm
          r3t(idx2(1),idx1(2),idx2(3))=r3t(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*sim%mass_p_cdm
          r3t(idx2(1),idx2(2),idx1(3))=r3t(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*sim%mass_p_cdm
          r3t(idx2(1),idx2(2),idx2(3))=r3t(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*sim%mass_p_cdm
        enddo
      enddo
      enddo
      enddo
      !$omp endparalleldo
    enddo
    ! put center part of r3t into subset of r3
    r3((itx-1)*nt+1:itx*nt,(ity-1)*nt+1:ity*nt,(itz-1)*nt+1:itz*nt)=r3t(1:nt,1:nt,1:nt)
  enddo
  enddo
  enddo
  testrhoc=sum(r3*1d0)
  call system_clock(t2,t_rate)
  print*, '    elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all
  ! coarse force ---------------------------------------------------------------
  if (head) print*,'    coarse cic force'
  if (head) print*,'      coarse_fft'
  call system_clock(t1,t_rate)
  call pencil_fft_forward
  ! save complex rho_c into crho_c
  crho_c(::2,:,:)=real(cxyz)
  crho_c(2::2,:,:)=imag(cxyz)
  do i_dim=1,3
    if (head) print*,'      coarse_ifft dim',int(i_dim,1)
    !$omp workshare
    rxyz(::2,:,:)=-crho_c(2::2,:,:)*kern_c(:,:,:,i_dim)
    rxyz(2::2,:,:)=crho_c(::2,:,:)*kern_c(:,:,:,i_dim)
    !$omp endworkshare
    call pencil_fft_backward
    force_c(i_dim,1:nc,1:nc,1:nc)=r3
  enddo

#ifdef record_Fc
  write(str_a,'(f7.3)') a_mid
  write(65,*) str_a
  open(66,file=opath//'image1/'//trim(adjustl(str_a))//'_Fc'//output_suffix(),status='replace',access='stream')
  write(66) force_c(1,1:nc,1:nc,1:nc)
  write(66) force_c(2,1:nc,1:nc,1:nc)
  write(66) force_c(3,1:nc,1:nc,1:nc)
  close(66)
#endif

  call system_clock(t2,t_rate)
  print*, '      elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all
  ! sync force_c buffer for CIC force
  if (head) print*, '      sync force_c buffer'
  call system_clock(t1,t_rate)
  force_c(:,0,:,:)=force_c(:,nc,:,:)[image1d(inx,icy,icz)]
  force_c(:,nc+1,:,:)=force_c(:,1,:,:)[image1d(ipx,icy,icz)]
  sync all
  force_c(:,:,0,:)=force_c(:,:,nc,:)[image1d(icx,iny,icz)]
  force_c(:,:,nc+1,:)=force_c(:,:,1,:)[image1d(icx,ipy,icz)]
  sync all
  force_c(:,:,:,0)=force_c(:,:,:,nc)[image1d(icx,icy,inz)]
  force_c(:,:,:,nc+1)=force_c(:,:,:,1)[image1d(icx,icy,ipz)]
  sync all
  ! coarse_max_dt
  f2_max_coarse=maxval(sum(force_c**2,1))
  sync all
  call system_clock(t2,t_rate)
  print*, '      elapsed time =',real(t2-t1)/t_rate,'secs';
  ! coarse velocity ------------------------------------------------------------
  if (head) print*, '    coarse cic velocity'
  call system_clock(t1,t_rate)
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    !$omp paralleldo default(shared) schedule(static,2)&
    !$omp& private(k,j,i,nlast,np,l,ip,tempx,idx1,idx2,dx1,dx2,vreal) &
    !$omp& reduction(max:vmax,vmax_nu)
    do k=1,nt
    do j=1,nt
    do i=1,nt
      np=rhoc(i,j,k,itx,ity,itz)
      nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
      do l=1,np ! loop over cdm particles
        ip=nzero+l
        tempx=((/itx,ity,itz/)-1)*nt+((/i,j,k/)-1)+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution-0.5
        idx1(:)=floor(tempx(:))+1
        idx2(:)=idx1(:)+1
        dx1(:)=idx1(:)-tempx(:)
        dx2(:)=1-dx1(:)
        vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
        vreal=vreal+force_c(:,idx1(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx1(3)
        vreal=vreal+force_c(:,idx2(1),idx1(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx1(3)
        vreal=vreal+force_c(:,idx1(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx1(3)
        vreal=vreal+force_c(:,idx1(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx1(2)*dx2(3)
        vreal=vreal+force_c(:,idx1(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx1(1)*dx2(2)*dx2(3)
        vreal=vreal+force_c(:,idx2(1),idx1(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx1(2)*dx2(3)
        vreal=vreal+force_c(:,idx2(1),idx2(2),idx1(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx1(3)
        vreal=vreal+force_c(:,idx2(1),idx2(2),idx2(3))*a_mid*dt/6/pi*dx2(1)*dx2(2)*dx2(3)
        vmax=max(vmax,abs(vreal+vfield(:,i,j,k,itx,ity,itz)))
        vp(:,ip)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
      enddo

      if (neutrino_flag) then
         tempx=((/itx,ity,itz/)-1)*nt+((/i,j,k/))
         idx1(:)=floor(tempx(:))+1
         vreal=force_c(:,idx1(1),idx1(2),idx1(3))*a_mid*dt/6/pi
         do kk=1,ncell
         do jj=1,ncell
         do ii=1,ncell
            tempx=((/itx,ity,itz/)-1)*nt*ncell+((/i,j,k/)-1)*ncell+((/ii,jj,kk/))
            do nu=1,Nneu
               call neu(nu)%gravity(tempx,vreal,a_mid)
            end do
         end do
         end do
         end do
      end if

    enddo
    enddo
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  call system_clock(t2,t_rate)
  print*, '      elapsed time =',real(t2-t1)/t_rate,'secs';

  sim%vsim2phys=(1.5/a)*box*h0*100.*sqrt(omega_m)/nf_global
  sync all

  if (head) print*, '  constrain dt'
  sim%dt_fine=sqrt( 1.0 / (sqrt(maxval(f2_max_fine))*a_mid*GG) )
  sim%dt_coarse=sqrt( real(ncell) / (sqrt(f2_max_coarse)*a_mid*GG) )
  sim%dt_vmax=vbuf*20/maxval(vmax)
  sim%vz_max=vmax(3)
  if (neutrino_flag) sim%dt_vmax_nu=vbuf*20/maxval(vmax_nu)
  sync all

  do i=1,nn**3
    sim%dt_fine=min(sim%dt_fine,sim[i]%dt_fine)
    sim%dt_coarse=min(sim%dt_coarse,sim[i]%dt_coarse)
    sim%dt_vmax=min(sim%dt_vmax,sim[i]%dt_vmax)
    if (neutrino_flag) sim%dt_vmax_nu=min(sim%dt_vmax_nu,sim[i]%dt_vmax_nu)
  enddo
  if (head) then
    call system_clock(tt2,t_rate)
    print*, '  mass_f,c =',testrhof,testrhoc
    print*, '  elapsed time =',real(tt2-tt1)/t_rate,'secs'
    print*, ''
  endif
  sync all
endsubroutine

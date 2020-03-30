!#define onehalo
#define sigma_8
!#define READ_SEED
!#define READ_NOISE
!#define READ_RECO
!#define filter_phi
!#define force_power

!#define ELUCID
!#define IniC_R0
!#define IniC_E

program initial_conditions
  use omp_lib
  use variables, only: spine_tile
  use pencil_fft
  use iso_fortran_env, only : int64
# ifdef force_power
    use powerspectrum
# endif
  implicit none
  save

  ! nc: coarse grid per image per dim
  ! nf: fine grid per image per dim, ng=nf
  ! nyquest: Nyquest frequency
  logical,parameter :: correct_kernel=.true.
  logical,parameter :: write_potential=.true.
  logical,parameter :: write_noise=.false.
  logical,parameter :: write_velocity=.true.

#ifdef sigma_8
  real, parameter :: s8 = 0.090
  character(len=*), parameter :: fntf = '../../tf/neutrino_tf/mnu3x100meV/mnu3x100meV_transfer_out_z10.dat'
  integer(8),parameter :: nk=1675
  real :: tf(13,nk)
#else
  integer(8),parameter :: nk=132
  real tf(14,nk)
#endif

#ifdef onehalo
  integer pid_halo(100000),np_halo,ipos(3),nw,iw
  real xv(6,100000),xv_mean(6),ang_mom(3)
#endif

  complex, allocatable :: cxyz_r(:,:,:)
#ifdef ELUCID
# define delta_small
  integer,parameter :: nr=500
  character(*),parameter :: file_deltak='../../S500_5001/cxyz_251_500_500.bin'
#endif

#ifdef IniC_R0
# define delta_big
  integer,parameter :: nr=400
  character(*),parameter :: file_deltak='../../IniC/cxyz_201_400_400.bin'
#endif
#ifdef IniC_R1
# define delta_big
  integer,parameter :: nr=200
  character(*),parameter :: file_deltak='../../IniC/cxyz_101_200_200.bin'
#endif
#ifdef IniC_R2
# define delta_big
  integer,parameter :: nr=300
  character(*),parameter :: file_deltak='../../IniC/cxyz_151_300_300.bin'
#endif

  integer(8) i,j,k,ip,l,nzero
  integer(8) ind,dx,dxy,kg,mg,jg,ig,ii,jj,kk,itx,ity,itz,idx,imove,g(3),iq(3)
  integer(4) seedsize,t1,t2,tt1,tt2,ttt1,ttt2,t_rate,ilayer,nlayer
  real a_i,kr,kx,ky,kz,kmax,temp_r,temp_theta,pow,phi8,temp8[*]
  real(8) v8,norm,xq(3),gradphi(3),vreal(3),dvar[*],dvarg
  integer(int64) time64

  integer(4),allocatable :: iseed(:)
  real,allocatable :: rseed_all(:,:)

  complex delta_k(nyquest+1,nf,npen)
  real phi(-nfb:nf+nfb+1,-nfb:nf+nfb+1,-nfb:nf+nfb+1)[*]

  ! zip arrays
  integer(8),parameter :: npt=nt*np_nc ! np / tile / dim !64
  integer(8),parameter :: npb=ncb*np_nc !24
  integer(8),parameter :: npmax=2*(npt+2*npb)**3
  integer(4) rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(4) rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  real(4) vfield(3,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(8) idx_ex_r(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(8),dimension(nt,nt) :: pp_l,pp_r,ppe_l,ppe_r

  integer(izipx) xp(3,npmax)
  integer(izipv) vp(3,npmax)
#ifdef PID
    integer(4) pid(npmax)
#endif
#ifdef force_power
  real xi(10,nbin)[*]
#endif
  real grad_max(3)[*],vmax(3),vf
  real(4) svz(500,2),svr(100,2)
  real(8) sigma_vc,sigma_vf
  real(8) std_vsim_c,std_vsim_res,std_vsim
  character (10) :: img_s, z_s

  call system_clock(ttt1,t_rate)
  call omp_set_num_threads(ncore)
  call geometry
  call system('mkdir -p '//opath//'image'//image2str(image))

#ifdef onehalo
  nw=0
#endif

  cur_checkpoint=1
  open(16,file='../main/z_checkpoint.txt',status='old')
  read(16,fmt='(f8.4)') z_checkpoint(cur_checkpoint)
# ifdef force_power
    z_checkpoint(cur_checkpoint)=5.0
# endif
  close(16)

  if (head) then
    print*, ''
    print*, 'CUBE Initial Conditions'
    print*, 'on',nn**3,' images'
    print*, 'Resolution', ng*nn
    print*, 'at redshift=',z_checkpoint(cur_checkpoint)
    if (body_centered_cubic) then
      print*, 'To genterate npglobal = 2 x', int(np_nc*nc*nn,2),'^3'
    else
      print*, 'To genterate npglobal =', int(np_nc*nc*nn,2),'^3'
    endif
    print*, 'Box size', box
    print*, 'body_centered_cubic =',body_centered_cubic
    print*, 'output: ', opath
    print*, 'head image number',icx,icy,icz
    print*, '-----------------------------------------'
    call system('mkdir -p '//opath//'code')
    call system('cp initial_conditions*.f90 '//opath//'code/')
    call system('cp ../main/*.f90 '//opath//'code/')
    call system('cp ../main/z_*.txt '//opath//'code/')
  endif

  sync all

  sim%nplocal=0
  sim%nplocal_nu=0
  sim%a=1./(1+z_checkpoint(cur_checkpoint))
  sim%t=0
  sim%tau=0

  sim%timestep=0

  sim%dt_pp=1000
  sim%dt_fine=1000
  sim%dt_coarse=1000
  sim%dt_vmax=1000
  sim%dt_vmax_nu=1000

  sim%cur_checkpoint=0
  sim%box=box
  sim%image=image
  sim%nn=nn
  sim%nnt=nnt
  sim%nt=nt
  sim%ncell=ncell
  sim%ncb=ncb
  sim%izipx=izipx
  sim%izipv=izipv
  sim%izipx_nu=izipx_nu
  sim%izipv_nu=izipv_nu

  sim%h0=h0
  sim%omega_m=omega_m
  sim%omega_l=omega_l
  sim%s8=s8
  sim%vsim2phys=(150./sim%a)*box*sqrt(omega_m)/nf_global
  sim%z_i=z_checkpoint(cur_checkpoint)
  sim%z_i_nu=z_i_nu
  sync all
  phi=0
  tf=0
  sync all

  if (head) print*,''
  if (head) print*,'Creating FFT plans'
  call system_clock(t1,t_rate)
  call create_penfft_plan
  call system_clock(t2,t_rate)
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all

  ! transferfnc --------------------------------------
  if (head) print*,''
  if (head) print*,'Transfer function'
  call system_clock(t1,t_rate)
#ifdef sigma_8
  open(11,file=fntf,form='formatted')
  read(11,*) tf
  close(11)
  ! normalization
  norm=1
  if (head) print*, 'Normalization factor: norm =', norm
  tf(2,:)=tf(8,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
  tf(3,:)=tf(3,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
  tf(6,:)=tf(6,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
  tf(7,:)=tf(7,:)**2.0 * tf(1,:)**(3+n_s) * norm / (2.0*pi**2)
  !dk
  tf(4,1)=tf(1,2)/2
  do k=2,nk-1
    tf(4,k)=(tf(1,k+1)-tf(1,k-1))/2
  enddo
  tf(4,nk)=tf(1,nk)-tf(1,nk-1)
  v8=0
  kmax=2*pi*sqrt(3.)*nyquest/box
  do k=1,nk
    if (tf(1,k)>kmax) exit
    v8=v8+tf(7,k)*tophat(tf(1,k)*8)**2*tf(4,k)/tf(1,k)
  enddo
  if (head) print*, 's8**2/v8:', v8, s8**2/v8,nyquest
  tf(2:3,:)=tf(2:3,:)*(s8**2/v8)*DgrowRatio(sim%z_i,sim%z_i_nu)**2

  if (head) then
     print*,'omega_neu,f_neu',omega_neu,f_neu
     print*,'sum f_neu',sum(f_neu)
     print*,'z_i,z_i_nu,DgrowRatio',sim%z_i,sim%z_i_nu,DgrowRatio(sim%z_i,sim%z_i_nu)
  end if

  sync all

#else
  ! remark: requires "CLASS" format for tf ("CAMB"="CLASS"/(-k^2) with k in 1/Mpc)
  open(11,file='../../tf/caf_z10_tk.dat',form='formatted')
  read(11,*) !header
  read(11,*) tf
  close(11)
  ! replace T_g with T_cb = f_c T_c + f_b T_b
  tf(2,:) = (omega_bar*tf(3,:)+omega_cdm*tf(4,:))/(omega_bar+omega_cdm)
  ! compute power spectrum @ z_tf
  tf(2,:) = A_s*(tf(1,:)/k_o)**(n_s-1.)*tf(2,:)**2
  ! propagate to starting redshift
  tf(2,:) = tf(2,:)*DgrowRatio(z_i,z_tf)**2

  sync all
#endif

  ! noisemap -------------------------------------
  if (head) print*,''
  if (head) print*,'Generating random noise'
  call random_seed(size=seedsize)
  if (head) print*,'  min seedsize =', seedsize
  seedsize=max(seedsize,12)
  allocate(iseed(seedsize))
  allocate(rseed_all(seedsize,nn**3))
#ifdef READ_SEED
    if (head) print*, '  Copy and read seeds from ../confings/'
    !call system('cp ../../configs/seed_'//image2str(image)//'.bin '//opath//'image'//image2str(image))
    open(11,file=output_dir()//'seed'//output_suffix(),status='old',access='stream')
    read(11) iseed
    close(11)
    ! Input iseed
    call random_seed(put=iseed)
    if (head) print*, 'iseed', iseed
#else
    ! Generate at least 12 seeds according to system clock
    call system_clock(time64)
    do i = 1, seedsize
      iseed(i) = lcg(time64) + image*137
      !print*,'time64,iseed(',int(i,1),')=',time64,iseed(i)
    enddo
    ! Input iseed to system
    call random_seed(put=iseed)
    !print*, 'iseed', iseed
    ! Write iseed into file
    open(11,file=output_dir()//'seed'//output_suffix(),status='replace',access='stream')
    write(11) iseed
    close(11)
    ! execute the following line if you want to save seeds
    !call system('cp ../output/universe1/image*/seed* ../configs/')
#endif

  call random_number(r3)
  deallocate(iseed)
  deallocate(rseed_all)
  sync all
# ifdef READ_NOISE
    open(11,file=output_dir()//'noise'//output_suffix(),access='stream')
    read(11) r3
    close(11)
    print*, '  READ IN NOISE MAP:', r3(1,1,1), r3(ng,ng,ng)
# else
    open(11,file=output_dir()//'noise'//output_suffix(),status='replace',access='stream')
    write(11) r3
    close(11)
    print*, '  noise',int(image,1),r3(1:2,1,1)
# endif
  sync all

  ! Box-Muller transform ----------------------------------------------
  if (head) print*,'  Box-Muller transform'
  !$omp paralleldo&
  !$omp& default(shared) &
  !$omp& private(k,j,i,temp_theta,temp_r)
  do k=1,nf
  do j=1,nf
  do i=1,nf,2
    temp_theta=2*pi*r3(i,j,k)
    temp_r=sqrt(-2*log(1-r3(i+1,j,k)))
    r3(i,j,k)=temp_r*cos(temp_theta)
    r3(i+1,j,k)=temp_r*sin(temp_theta)
  enddo
  enddo
  enddo
  !$omp endparalleldo
  sync all

  if (write_noise) then
     open(11,file=output_dir()//'gaussian_noise'//output_suffix(),status='replace',access='stream')
     write(11) r3
     close(11)
     print*, '  noise',int(image,1),r3(1:2,1,1)
     sync all
  end if

  call system_clock(t2,t_rate)
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  ! delta_field ----------------------------------------------------
  if (head) print*, ''
  if (head) print*, 'Delta field'
  call system_clock(t1,t_rate)
  if (head) print*, '  ftran'
  call pencil_fft_forward
  if (head) print*, '  Wiener filter'
  !$omp paralleldo&
  !$omp& default(shared) &
  !$omp& private(k,j,i,kg,jg,ig,kz,ky,kx,kr,pow)
  do k=1,npen
  do j=1,nf
  do i=1,nyquest+1
    ! global grid in Fourier space for i,j,k
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nyquest-1,nf_global)-nyquest
    ky=mod(jg+nyquest-1,nf_global)-nyquest
    kx=ig-1
    kr=sqrt(kx**2+ky**2+kz**2)
    kr=max(kr,1.0)
    !pow=interp_tf(2*pi*kr/box,1,2)/(4*pi*kr**3)
    !cxyz(i,j,k)=cxyz(i,j,k)*sqrt(pow*nf_global*nf_global*nf_global)
    cxyz(i,j,k)=cxyz(i,j,k)*sqrt(interp_tf(2*pi*kr/box,1,2)/(4*pi)*(nf_global/kr))*(nf_global/kr)
    if (kr==1) then
      print*,cxyz(i,j,k)
    endif
  enddo
  enddo
  enddo
  !$omp endparalleldo

#ifdef delta_small
  print*, 'replace cxyz with the center part of cxyz_r'
  allocate(cxyz_r(nr/2+1,nr,nr))
  open(11,file=file_deltak,status='old',access='stream')
  read(11) cxyz_r
  cxyz(:ng/2+1,:ng/2+1,:ng/2+1)=cxyz_r(:ng/2+1,:ng/2+1,:ng/2+1)
  cxyz(:ng/2+1,ng/2+2:,:ng/2+1)=cxyz_r(:ng/2+1,nr-ng/2+2:,:ng/2+1)
  cxyz(:ng/2+1,:ng/2+1,ng/2+2:)=cxyz_r(:ng/2+1,:ng/2+1,nr-ng/2+2:)
  cxyz(:ng/2+1,ng/2+2:,ng/2+2:)=cxyz_r(:ng/2+1,nr-ng/2+2:,nr-ng/2+2:)
  close(11)
  deallocate(cxyz_r)
  cxyz=cxyz*ng*ng*ng*Dgrow(sim%a)
#endif

#ifdef delta_big
  print*, 'replace the center part of cxyz with cxyz_r'
  allocate(cxyz_r(nr/2+1,nr,nr))
  open(11,file=file_deltak,status='old',access='stream')
  read(11) cxyz_r
  cxyz_r=cxyz_r*ng*ng*ng*Dgrow(sim%a)
  cxyz(:nr/2,:nr/2,:nr/2)=cxyz_r(:nr/2,:nr/2,:nr/2)
  cxyz(:nr/2,ng-nr/2+2:,:nr/2)=cxyz_r(:nr/2,nr/2+2:,:nr/2+1)
  cxyz(:nr/2,:nr/2,ng-nr/2+2:)=cxyz_r(:nr/2,:nr/2+1,nr/2+2:)
  cxyz(:nr/2,ng-nr/2+2:,ng-nr/2+2:)=cxyz_r(:nr/2,nr/2+2:,nr/2+2:)
  close(11)
  deallocate(cxyz_r)
#endif

#ifdef IniC_E
  print*,'use delta_E'
  open(11,file='/mnt/raid-cita/haoran/CUBEnu/output/universe16/image1/0.000_delta_E_1.bin',access='stream')
  read(11) r3
  close(11); sync all
  r3=r3*Dgrow(sim%a)
  call pencil_fft_forward
#endif

#ifdef READ_RECO
  open(11,file='/home/yuyu22/0.000_recon_1.bin',access='stream')
  read(11) r3
  close(11)
  r3=r3*Dgrow(sim%a)
  call pencil_fft_forward
#endif


  if (head) cxyz(1,1,1)=0 ! DC frequency
  sync all
  delta_k=cxyz ! backup k-space delta_L

  if (head) print*,'  btran'
  call pencil_fft_backward
  print*,'  delta_L',r3(1:4,1,1)
  print*,'  rms of delta',sqrt(sum(r3**2*1.d0)/nf_global/nf_global/nf_global)

  if (head) print*,'  write delta_L into file'
  if (head) print*,'  growth factor Dgrow(',sim%a,') =',Dgrow(sim%a)
  open(11,file=output_dir()//'delta_L'//output_suffix(),status='replace',access='stream')
  write(11) r3/DgrowRatio(sim%z_i,sim%z_i_nu)
  close(11); sync all

  open(11,file=output_dir()//'delta_L_proj'//output_suffix(),status='replace',access='stream')
  write(11) sum(r3(:,:,:13),dim=3)/Dgrow(sim%a)
  close(11); sync all
  call system_clock(t2,t_rate)
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  ! Potential field ----------------------------------------------------
  if (head) print*, ''
  if (head) print*, 'Potential field'
  call system_clock(t1,t_rate)
  !$omp paralleldo&
  !$omp& default(shared) &
  !$omp& private(k,j,i,kg,jg,ig,kz,ky,kx,kr)
  do k=1,npen
  do j=1,nf
  do i=1,nyquest+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nyquest-1,nf_global)-nyquest
    ky=mod(jg+nyquest-1,nf_global)-nyquest
    kx=ig-1
    kz=2*sin(pi*kz/nf_global)
    ky=2*sin(pi*ky/nf_global)
    kx=2*sin(pi*kx/nf_global)
    kr=kx**2+ky**2+kz**2
    kr=max(kr,1.0/nf_global**2) ! avoid kr being 0
    cxyz(i,j,k)=-4*pi/kr
  enddo
  enddo
  enddo
  !$omp endparalleldo
  if (head) cxyz(1,1,1)=0 ! DC frequency
  sync all

  if (correct_kernel) then
    if (head) print*, '  correct kernel'
    call pencil_fft_backward
    temp8=0
    if (image==1) temp8=temp8+r3(9,1,1)+r3(1,9,1)+r3(1,1,9)
    sync all
    if (icx==nn .and. icy==1 .and. icz==1) temp8=temp8+r3(nf-7,1,1)
    sync all
    if (icx==1 .and. icy==nn .and. icz==1) temp8=temp8+r3(1,nf-7,1)
    sync all
    if (icx==1 .and. icy==1 .and. icz==nn) temp8=temp8+r3(1,1,nf-7)
    sync all
    phi8=0
    do i=1,nn**3
      phi8=phi8+temp8[i]
      sync all
    enddo
    sync all
    phi8=phi8/6
    if (head) print*,'  phi8 =',phi8
    sync all
    if (head) print*, '  construct Ewald potential kernel in real space'
    !$omp paralleldo&
    !$omp& default(shared) &
    !$omp& private(k,j,i,kg,jg,ig,kx,ky,kz,kr)
    do k=1,nf
    do j=1,nf
    do i=1,nf
      kg=k+nf*(icz-1)
      jg=j+nf*(icy-1)
      ig=i+nf*(icx-1)
      kx=mod(kg+nyquest-1,nf_global)-nyquest
      ky=mod(jg+nyquest-1,nf_global)-nyquest
      kz=mod(ig+nyquest-1,nf_global)-nyquest
      kr=sqrt(kx**2+ky**2+kz**2)
      if (kr>8) then
        r3(i,j,k)=r3(i,j,k)-(phi8+1/8.)
      elseif (kr>0) then
        r3(i,j,k)=-1/kr
      elseif (kr==0) then
        r3(i,j,k)=-2.5
      endif
    enddo
    enddo
    enddo
    !$omp endparalleldo
    sync all
    call pencil_fft_forward
  endif

  ! Complex multiply delta_L with potential kernel
  cxyz=real(cxyz)*delta_k
  delta_k=cxyz  ! backup phi(k)
  call pencil_fft_backward

  phi=0
  phi(1:nf,1:nf,1:nf)=r3 ! phi1
  print*,'  phi',phi(1:4,1,1)
#ifdef READ_RECO
  open(11,file=output_name('phireco'),status='replace',access='stream')
  write(11) r3
  close(11)
  stop
#endif
  if (write_potential) then
    if (head) print*, '  write phi1 into file'
    open(11,file=output_name('phi1'),status='replace',access='stream')
    write(11) r3
    close(11)
  endif

  sync all

  ! buffer phi ---------------------------------------------------
  if (head) print*, '  buffer phi'
  phi(:0,:,:)=phi(nf-nfb:nf,:,:)[image1d(inx,icy,icz)]
  sync all
  phi(nf+1:,:,:)=phi(1:nfb+1,:,:)[image1d(ipx,icy,icz)]
  sync all
  phi(:,:0,:)=phi(:,nf-nfb:nf,:)[image1d(icx,iny,icz)]
  sync all
  phi(:,nf+1:,:)=phi(:,1:nfb+1,:)[image1d(icx,ipy,icz)]
  sync all
  phi(:,:,:0)=phi(:,:,nf-nfb:nf)[image1d(icx,icy,inz)]
  sync all
  phi(:,:,nf+1:)=phi(:,:,1:nfb+1)[image1d(icx,icy,ipz)]
  sync all

#ifdef filter_phi
  cxyz=0
  do k=1,npen
  do j=1,nf
  do i=1,nyquest+1
    ! global grid in Fourier space for i,j,k
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nyquest-1,nf_global)-nyquest
    ky=mod(jg+nyquest-1,nf_global)-nyquest
    kx=ig-1
    kr=sqrt(kx**2+ky**2+kz**2)
    kr=max(kr,1.0)
    kr=2*pi*kr/box
    pow=exp(-kr**2*1.8**2/2)**0.25 ! apply E-mode window function
    cxyz(i,j,k)=delta_k(i,j,k)*pow
  enddo
  enddo
  enddo
  if (head) cxyz(1,1,1)=0 ! DC frequency
  call pencil_fft_backward
  if (head) print*, '  write filtered phi1 into file'
  if (head) print*, '  ',output_name('phiE1')
  open(11,file=output_name('phiE1'),status='replace',access='stream')
  write(11) r3
  close(11)
!#   ifdef force_power
      ! now diff and compute F_x
!      do k=1,nf
!      do j=1,nf
!      do i=1,nf
!        r3(i,j,k)=-(phi(i+1,j,k)-phi(i-1,j,k))/2
!      enddo
!      enddo
!      enddo
!      ! compute power spectrum of F_x
!      call cross_power(xi,r3,r3)
!      print*, 'called cross_power'
!      sync all
!      if (head) then
!        open(15,file=output_name('power_Fx'),status='replace',access='stream')
!        write(15) xi
!        close(15)
!      endif
!      print*,'wrote',output_name('power_Fx')
!      stop
!#   endif

  cxyz=0
  do k=1,npen
  do j=1,nf
  do i=1,nyquest+1
    ! global grid in Fourier space for i,j,k
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*nf+j
    ig=i
    kz=mod(kg+nyquest-1,nf_global)-nyquest
    ky=mod(jg+nyquest-1,nf_global)-nyquest
    kx=ig-1
    kr=sqrt(kx**2+ky**2+kz**2)
    kr=max(kr,1.0)
    kr=2*pi*kr/box
    pow=exp(-kr**2*3.0**2/2)**0.25 ! apply E-mode window function
    cxyz(i,j,k)=delta_k(i,j,k)*pow
  enddo
  enddo
  enddo
  if (head) cxyz(1,1,1)=0 ! DC frequency
  call pencil_fft_backward
  if (head) print*, '  write filtered phi1 into file'
  if (head) print*, '  ',output_name('phiE2')
  open(11,file=output_name('phiE2'),status='replace',access='stream')
  write(11) r3
  close(11)
#endif

  if (head) print*, '  destroying FFT plans'
  call destroy_penfft_plan
  sync all
  call system_clock(t2,t_rate)
  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all

  !write velocity field
  if (write_velocity) then

     vf=vfactor(1./(1.+sim%z_i_nu))/( -(1./4.)+(5./4.)*sqrt(1-24.*sum(f_neu)/25.) ) !vfactor at z_i_nu

     !!x
     r3=0
     do k=1,nf
        do j=1,nf
           do i=1,nf
              r3(i,j,k)=-(phi(i+1,j,k)-phi(i-1,j,k))/(8*pi)*vf/DgrowRatio(sim%z_i,sim%z_i_nu)
           end do
        end do
     end do

     if (head) print*, '  write v_x into file'
     open(11,file=output_dir()//'v_x_L'//output_suffix(),status='replace',access='stream')
     write(11) r3
     close(11)

     !!y
     r3=0
     do k=1,nf
        do j=1,nf
           do i=1,nf
              r3(i,j,k)=-(phi(i,j+1,k)-phi(i,j-1,k))/(8*pi)*vf/DgrowRatio(sim%z_i,sim%z_i_nu)
           end do
        end do
     end do

     if (head) print*, '  write v_y into file'
     open(11,file=output_dir()//'v_y_L'//output_suffix(),status='replace',access='stream')
     write(11) r3
     close(11)

     !!z
     r3=0
     do k=1,nf
        do j=1,nf
           do i=1,nf
              r3(i,j,k)=-(phi(i,j,k+1)-phi(i,j,k-1))/(8*pi)*vf/DgrowRatio(sim%z_i,sim%z_i_nu)
           end do
        end do
     end do

     if (head) print*, '  write v_z into file'
     open(11,file=output_dir()//'v_z_L'//output_suffix(),status='replace',access='stream')
     write(11) r3
     close(11)

  end if

  ! zip checkpoints ------------------------------------------------
  if (head) print*, ''
  if (head) print*, 'zip checkpoints'
  vf=vfactor(sim%a)
  if (head) print*, '  vf =',vf
  grad_max(1)=maxval(abs(phi(-nfb:nf+nfb-1,:,:)-phi(-nfb+2:nf+nfb+1,:,:)))
  grad_max(2)=maxval(abs(phi(:,-nfb:nf+nfb-1,:)-phi(:,-nfb+2:nf+nfb+1,:)))
  grad_max(3)=maxval(abs(phi(:,:,-nfb:nf+nfb-1)-phi(:,:,-nfb+2:nf+nfb+1)))
  sync all
  do i=1,nn**3 ! co_max
    grad_max=max(grad_max,grad_max(:)[i])
    sync all
  enddo
  sync all
  vmax=grad_max/2/(4*pi)*vf

  sim%dt_vmax=vbuf*20./maxval(abs(vmax))
  sim%vz_max=vmax(3)
  nlayer=2*ceiling(grad_max(3)*np_nc/8/pi/ncell)+1
  if (head) then
    print*, '  grad_max',grad_max
    print*, '  max dsp',grad_max/2/(4*pi)
    print*, '  vmax',vmax
    print*, '  vz_max',sim%vz_max
    if (maxval(grad_max)/2/(4*pi)>=nfb) then
      print*, '  particle dsp > buffer'
      print*, maxval(grad_max)/2/(4*pi),nfb
      stop
    endif
    print*,'  Thread save nlayer =',nlayer
    print*,''
  endif
  sync all

  open(11,file='../../velocity_conversion/sigmav_z.bin',access='stream')
  read(11) svz
  close(11)
  open(11,file='../../velocity_conversion/sigmav_r.bin',access='stream')
  read(11) svr
  close(11)

  sigma_vf=interp_sigmav(sim%a,box/nf_global) ! sigma(v) on scale of fine grid, in km/s
  sigma_vc=interp_sigmav(sim%a,box/nc_global) ! sigma(v) on scale of coarse grid, in km/s
  sim%sigma_vres=sqrt(sigma_vf**2-sigma_vc**2) ! sigma(v) residual, in km/s
  sim%sigma_vi=sim%sigma_vres/sim%vsim2phys/sqrt(3.) ! sigma(v_i) residual, in sim unit
  sync all
  if (head) then
    print*, ''
    print*, 'Read velocity dispersion prediction'
    print*,'sigma_vf(a=',sim%a,', r=',box/nf_global,'Mpc/h)=',real(sigma_vf,4),'km/s'
    print*,'sigma_vc(a=',sim%a,', r=',box/nc_global,'Mpc/h)=',real(sigma_vc,4),'km/s'
    print*,'sigma_vres=',real(sim%sigma_vres,4),'km/s'
    print*,'sigma_vi =',real(sim%sigma_vi,4),'(simulation unit)'
  endif

  sync all

  ! create particles (no communication) ----------------------------
  if (head) print*,''
  if (head) print*, 'Create particles'
  call system_clock(t1,t_rate)
  open(11,file=output_name('xp'),status='replace',access='stream')
  open(12,file=output_name('vp'),status='replace',access='stream')
  open(13,file=output_name('np'),status='replace',access='stream')
  open(14,file=output_name('vc'),status='replace',access='stream')
#ifdef PID
  if (head) print*, '  also create PID'
  open(15,file=output_name('id'),status='replace',access='stream')
#endif

  if (body_centered_cubic .and. ncell/np_nc/2==0) stop 'ncell/np_nc/2 = 0, unsuitable for body centered cubic'

  vfield=0
  std_vsim_c=0; std_vsim_res=0; std_vsim=0; !np_prev=0

  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    rhoce=0
    rholocal=0

    do ilayer=0,nlayer-1
      !$omp paralleldo default(shared) schedule(static,8)&
      !$omp& private(k,j,i,imove,kk,jj,ii,xq,gradphi,g,vreal)
      do k=1-npb+ilayer,npt+npb,nlayer
      do j=1-npb,npt+npb
      do i=1-npb,npt+npb
      do imove=0,merge(1,0,body_centered_cubic)
        kk=nft*(itz-1)+(ncell/np_nc)*(k-1)+1+imove*(ncell/np_nc/2)
        jj=nft*(ity-1)+(ncell/np_nc)*(j-1)+1+imove*(ncell/np_nc/2)
        ii=nft*(itx-1)+(ncell/np_nc)*(i-1)+1+imove*(ncell/np_nc/2)
        xq=((/i,j,k/)-1d0)/np_nc + (0.5d0+imove*(ncell/np_nc/2))/ncell ! Lagrangian position q, in coarse grid
        gradphi(1)=phi(ii+1,jj,kk)-phi(ii-1,jj,kk)
        gradphi(2)=phi(ii,jj+1,kk)-phi(ii,jj-1,kk)
        gradphi(3)=phi(ii,jj,kk+1)-phi(ii,jj,kk-1)
        g=ceiling(xq-gradphi/(8*pi*ncell))
        rhoce(g(1),g(2),g(3))=rhoce(g(1),g(2),g(3))+1
        vreal=-gradphi/(8*pi)*vf
        vfield(:,g(1),g(2),g(3))=vfield(:,g(1),g(2),g(3))+vreal ! record vfield according to real particles
      enddo ! imove
      enddo
      enddo
      enddo !k
      !$omp endparalleldo
    enddo ! ilayer
    vfield(1,:,:,:)=vfield(1,:,:,:)/merge(1,rhoce,rhoce==0)
    vfield(2,:,:,:)=vfield(2,:,:,:)/merge(1,rhoce,rhoce==0)
    vfield(3,:,:,:)=vfield(3,:,:,:)/merge(1,rhoce,rhoce==0)
    call spine_tile(rhoce,idx_ex_r,pp_l,pp_r,ppe_l,ppe_r)

    do ilayer=0,nlayer-1
      !$omp paralleldo default(shared) schedule(static,8)&
      !$omp& private(k,j,i,imove,kk,jj,ii,xq,gradphi,g,idx,vreal,iq)
      do k=1-npb+ilayer,npt+npb,nlayer
      do j=1-npb,npt+npb
      do i=1-npb,npt+npb
      do imove=0,merge(1,0,body_centered_cubic)
        kk=nft*(itz-1)+(ncell/np_nc)*(k-1)+1+imove*(ncell/np_nc/2)
        jj=nft*(ity-1)+(ncell/np_nc)*(j-1)+1+imove*(ncell/np_nc/2)
        ii=nft*(itx-1)+(ncell/np_nc)*(i-1)+1+imove*(ncell/np_nc/2)
        xq=((/i,j,k/)-1d0)/np_nc + (0.5d0+imove*(ncell/np_nc/2))/ncell
        gradphi(1)=phi(ii+1,jj,kk)-phi(ii-1,jj,kk)
        gradphi(2)=phi(ii,jj+1,kk)-phi(ii,jj-1,kk)
        gradphi(3)=phi(ii,jj,kk+1)-phi(ii,jj,kk-1)
        g=ceiling(xq-gradphi/(8*pi*ncell))
        rholocal(g(1),g(2),g(3))=rholocal(g(1),g(2),g(3))+1
        idx=idx_ex_r(g(2),g(3))-sum(rhoce(g(1):,g(2),g(3)))+rholocal(g(1),g(2),g(3))
        xp(:,idx)=floor((xq-gradphi/(8*pi*ncell))/x_resolution,kind=8)
        vreal=-gradphi/(8*pi)*vf
        vreal=vreal-vfield(:,g(1),g(2),g(3)) ! save relative velocity
        vp(:,idx)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sim%sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
#       ifdef PID
          iq = ((/icx,icy,icz/)-1)*nf + ((/itx,ity,itz/)-1)*nft + (ncell/np_nc)*((/i,j,k/)-1)+imove
          iq = modulo(iq,nf_global)
          pid(idx)=iq(1)+nf_global*iq(2)+nf_global**2*iq(3)+1
          if (pid(idx)<1) stop
#       endif
      enddo ! imove
      enddo
      enddo
      enddo ! k
      !$omp endparalleldo
    enddo ! ilayer

    do k=1,nt ! delete buffer particles
    do j=1,nt
      xp(:,pp_l(j,k):pp_r(j,k))=xp(:,ppe_l(j,k):ppe_r(j,k))
      vp(:,pp_l(j,k):pp_r(j,k))=vp(:,ppe_l(j,k):ppe_r(j,k))
#     ifdef PID
        pid(pp_l(j,k):pp_r(j,k))=pid(ppe_l(j,k):ppe_r(j,k))
        if(minval(pid(pp_l(j,k):pp_r(j,k)))<1) then
          print*,'there is zero',k,j
          !print*,pp_l(j,k),pp_r(j,k)
          !print*,ppe_l(j,k),ppe_r(j,k)
          !print*,pid(ppe_l(j,k):ppe_r(j,k))
          !stop
        endif
#     endif
    enddo
    enddo

    ! velocity analysis
    !$omp paralleldo&
    !$omp& default(shared) &
    !$omp& private(k,j,i,nzero,l,ip,vreal)&
    !$omp& reduction(+:std_vsim_c,std_vsim_res,std_vsim)
    do k=1,nt
    do j=1,nt
    do i=1,nt
      nzero=pp_r(j,k)-sum(rhoce(i:nt,j,k))
      std_vsim_c=std_vsim_c+sum(vfield(:,i,j,k)**2)
      do l=1,rhoce(i,j,k)
        ip=nzero+l
        vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sim%sigma_vi*vrel_boost))
        std_vsim_res=std_vsim_res+sum(vreal**2)
        vreal=vreal+vfield(:,i,j,k)
        std_vsim=std_vsim+sum(vreal**2)
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo

    write(11) xp(:,1:pp_r(nt,nt))
    write(12) vp(:,1:pp_r(nt,nt))
    write(13) rhoce(1:nt,1:nt,1:nt)
    write(14) vfield(:,1:nt,1:nt,1:nt)
#   ifdef PID
      if (minval(pid(1:pp_r(nt,nt))) <1 ) then
         print*, 'nonpositive PID',itx,ity,itz
         print*, count(pid(1:pp_r(nt,nt))==0),minval(pid(1:pp_r(nt,nt)))
         stop
      endif
      write(15) pid(1:pp_r(nt,nt))
#   endif
    sim%nplocal=sim%nplocal+pp_r(nt,nt)
  enddo
  enddo
  enddo ! end of tile loop

  close(11)
  close(12)
  close(13)
  close(14)
# ifdef PID
    close(15)
# endif
  call system_clock(t2,t_rate)
  sync all

  if (head) print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  if (head) then
    print*,''
    print*,'Velocity analysis on head node'
    std_vsim_res=sqrt(std_vsim_res/sim%nplocal)
    std_vsim_c=sqrt(std_vsim_c/nc/nc/nc)
    std_vsim=sqrt(std_vsim/sim%nplocal)
    print*,'  std_vsim         ',real(std_vsim*sim%vsim2phys,4),'km/s'
    print*,'  std_vsim_c       ',real(std_vsim_c*sim%vsim2phys,4),'km/s'
    print*,'  std_vsim_res     ',real(std_vsim_res*sim%vsim2phys,4),'km/s'
    print*,'  std_vi (sim unit)',real(std_vsim_res/sqrt(3.),4),'(simulation unit)'
    print*,''
  endif
  sync all

  print*,'image',image,', nplocal',sim%nplocal
  sync all
  sim%npglobal=0
  do i=1,nn**3
    sim%npglobal=sim%npglobal+sim[i]%nplocal
    sync all
  enddo
  !sim%nplocal_nu=(np_nc_nu*nc)**3
  !sim%npglobal_nu=(np_nc_nu*nc*nn)**3
  if (head) then
    print*, 'npglobal =',sim%npglobal
    if (sim%npglobal/=merge(2,1,body_centered_cubic)*(np_nc*nc*nn)**3) then
      print*, 'warning: incorrect npglobal'
    endif
  endif
  sim%mass_p_cdm=real(f_cdm*nf_global**3,kind=8)/sim%npglobal
  sim%mass_p_nu=real(sum(f_neu)*nf_global**3,kind=8)/sim%npglobal_nu

#ifdef onehalo
  sim%mass_p_cdm=4
#endif
  call print_header(sim)

  sync all
  open(10,file=output_name('info'),status='replace',access='stream')
  write(10) sim
  close(10)
  sync all

  call system_clock(ttt2,t_rate)
  if (head) print*, 'total elapsed time =',real(ttt2-ttt1)/t_rate,'secs';
  if (head) print*, 'initial condition done'

  contains

  real function interp_sigmav(aa,rr)
    implicit none
    integer(8) ii,i1,i2
    real aa,rr,term_z,term_r
    i1=1
    i2=500
    do while (i2-i1>1)
      ii=(i1+i2)/2
      if (aa>svz(ii,1)) then
        i1=ii
      else
        i2=ii
      endif
    enddo
    term_z=svz(i1,2)+(svz(i2,2)-svz(i1,2))*(aa-svz(i1,1))/(svz(i2,1)-svz(i1,1))
    i1=1
    i2=100
    do while (i2-i1>1)
      ii=(i1+i2)/2
      if (rr>svz(ii,1)) then
        i1=ii
      else
        i2=ii
      endif
    enddo
    term_r=svr(i1,2)+(svr(i2,2)-svr(i1,2))*(rr-svr(i1,1))/(svr(i2,1)-svr(i1,1))
    interp_sigmav=term_z*term_r
  endfunction

  real function interp_tf(kr,ix,iy)
    implicit none
    integer(4) ix,iy
    integer(8) ii,i1,i2
    real kr,xx,yy,x1,x2,y1,y2
    i1=1
    i2=nk
    do while (i2-i1>1)
      ii=(i1+i2)/2
      if (kr>tf(ix,ii)) then
        i1=ii
      else
        i2=ii
      endif
    enddo
    x1=log(tf(ix,i1))
    y1=log(tf(iy,i1))
    x2=log(tf(ix,i2))
    y2=log(tf(iy,i2))
    xx=log(kr)
    yy=y1+(y2-y1)*(xx-x1)/(x2-x1)
    interp_tf=exp(yy)
  endfunction


  function tophat(x)
    implicit none
    real :: x,tophat
    !if (x/=0) then
    !  tophat=3*(sin(x)-cos(x)*x)/x**3
    !else
    !  tophat=1
    !endif
    tophat=merge(1.,3*(sin(x)-cos(x)*x)/x**3,x==0)
  endfunction tophat

  function Dgrow(a)
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real a
    real Dgrow
    real g,ga,hsq,oma,ola
    hsq=om/a**3+(1-om-ol)/a**2+ol
    oma=om/(a**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
    Dgrow=a*ga/g
  end function Dgrow

  function DgrowRatio(z1,z2) result(Dgrow)
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real, parameter :: np = -(1./4.)+(5./4.)*sqrt(1-24.*sum(f_neu)/25.) !~1-3f/5
    real z1,z2
    real Dgrow
    real hsq,oma,ola,a1,a2,ga1,ga2

    a1=1./(1.+z1)
    hsq=om/a1**3+(1-om-ol)/a1**2+ol
    oma=om/(a1**3*hsq)
    ola=ol/hsq
    ga1=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))

    a2=1./(1.+z2)
    hsq=om/a2**3+(1-om-ol)/a2**2+ol
    oma=om/(a2**3*hsq)
    ola=ol/hsq
    ga2=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))

    Dgrow=(a1*ga1)/(a2*ga2)
    Dgrow=Dgrow**np!(1.-3.*(omega_nu/omega_m)/5.)
  end function DgrowRatio

  function vfactor(a)
    implicit none
    real, parameter :: np = -(1./4.)+(5./4.)*sqrt(1-24.*sum(f_neu)/25.) !~1-3f/5
    real :: a
    real :: H,km,lm
    real :: vfactor
    lm=omega_l/omega_m
    km=(1-omega_m-omega_l)/omega_m
    H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
    vfactor=a**2*H
    vfactor=vfactor*np!(1.-3.*(omega_nu/omega_m)/5.)
  endfunction vfactor

  function lcg(s) !// Linear congruential generator
    implicit none
    integer(4) :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  endfunction
end

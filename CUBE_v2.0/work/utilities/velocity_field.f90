program cic_velpower
  use parameters
  use pencil_fft
  use powerspectrum
  implicit none
  save
  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim
  integer(8),parameter :: npnode=nf**3
  real,parameter :: density_buffer=1.2
  integer(8),parameter :: npmax=npnode*density_buffer
  integer(8) i,j,k,l,i_dim,iq(3),nplocal,nplocal_nu,itx,ity,itz
  integer(8) nlast,ip,np,idx1(3),idx2(3)

  real(4) rho0(0:ng+1,0:ng+1,0:ng+1)[*]
  real(4) vel(3,0:ng+1,0:ng+1,0:ng+1)[*]
  real(4) rho_c(ng,ng,ng),rho_nu(ng,ng,ng)
  real(4) vel_c(3,ng,ng,ng),vel_nu(3,ng,ng,ng)
  real(4) mass_p,pos1(3),dx1(3),dx2(3),vreal(3),sigma_vi,sigma_vi_nu
  real(8) rho8[*]

  integer(izipx),allocatable :: xp(:,:)
  integer(izipv),allocatable :: vp(:,:)
  integer(izipx_nu),allocatable :: xp_nu(:,:)
  integer(izipv_nu),allocatable :: vp_nu(:,:)
  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt),vfield(3,nt,nt,nt,nnt,nnt,nnt)
  real xi(10,nbin)[*]
  character(20) str_z,str_i

  call geometry
  if (head) then
    print*, 'cic_velpower on resolution:'
    print*, 'ng=',ng
    print*, 'ng*nn=',ng*nn
    print*, 'checkpoint at:'
    open(16,file='../main/redshifts.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
      print*, z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
    print*,''
  endif

  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]
  sync all

  call create_penfft_plan

  do cur_checkpoint= 1,n_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))

    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim; close(11)
    !call print_header(sim); stop
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) stop 'zip format incompatable'
    !mass_p=sim%mass_p
    mass_p=1. ! compute delta so mass_p will cancle out
    nplocal=sim%nplocal
    nplocal_nu=sim%nplocal_nu
    sigma_vi=sim%sigma_vi
    sigma_vi_nu=sim%sigma_vi_nu
    if (head) then
      print*, 'mass_p =',mass_p
      print*, 'nplocal =',nplocal
      print*, 'nplocal_nu =',nplocal_nu
    endif

    !cdm
    allocate(xp(3,nplocal),vp(3,nplocal))
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp; close(11)
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) vp; close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc; close(11)
    open(11,file=output_name('vc'),status='old',action='read',access='stream')
    read(11) vfield; close(11)
    rho0=0; vel=0; nlast=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l
          pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          pos1=pos1*real(ng)/real(nc) - 0.5
          idx1=floor(pos1)+1
          idx2=idx1+1
          dx1=idx1-pos1
          dx2=1-dx1
          ! density field
          rho0(idx1(1),idx1(2),idx1(3))=rho0(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p
          rho0(idx2(1),idx1(2),idx1(3))=rho0(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p
          rho0(idx1(1),idx2(2),idx1(3))=rho0(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p
          rho0(idx1(1),idx1(2),idx2(3))=rho0(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p
          rho0(idx1(1),idx2(2),idx2(3))=rho0(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p
          rho0(idx2(1),idx1(2),idx2(3))=rho0(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p
          rho0(idx2(1),idx2(2),idx1(3))=rho0(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p
          rho0(idx2(1),idx2(2),idx2(3))=rho0(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p
          ! velocity field
          vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
          vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
          vel(:,idx1(1),idx1(2),idx1(3))=vel(:,idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*vreal
          vel(:,idx2(1),idx1(2),idx1(3))=vel(:,idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*vreal
          vel(:,idx1(1),idx2(2),idx1(3))=vel(:,idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*vreal
          vel(:,idx1(1),idx1(2),idx2(3))=vel(:,idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*vreal
          vel(:,idx1(1),idx2(2),idx2(3))=vel(:,idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*vreal
          vel(:,idx2(1),idx1(2),idx2(3))=vel(:,idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*vreal
          vel(:,idx2(1),idx2(2),idx1(3))=vel(:,idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*vreal
          vel(:,idx2(1),idx2(2),idx2(3))=vel(:,idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*vreal
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all
    deallocate(xp,vp)

    if (head) print*, 'Start sync from buffer regions'
    sync all
    rho0(1,:,:)=rho0(1,:,:)+rho0(ng+1,:,:)[image1d(inx,icy,icz)]
    rho0(ng,:,:)=rho0(ng,:,:)+rho0(0,:,:)[image1d(ipx,icy,icz)]; sync all
    rho0(:,1,:)=rho0(:,1,:)+rho0(:,ng+1,:)[image1d(icx,iny,icz)]
    rho0(:,ng,:)=rho0(:,ng,:)+rho0(:,0,:)[image1d(icx,ipy,icz)]; sync all
    rho0(:,:,1)=rho0(:,:,1)+rho0(:,:,ng+1)[image1d(icx,icy,inz)]
    rho0(:,:,ng)=rho0(:,:,ng)+rho0(:,:,0)[image1d(icx,icy,ipz)]; sync all
    rho_c=rho0(1:ng,1:ng,1:ng)
    sync all
    vel(:,1,:,:)=vel(:,1,:,:)+vel(:,ng+1,:,:)[image1d(inx,icy,icz)]
    vel(:,ng,:,:)=vel(:,ng,:,:)+vel(:,0,:,:)[image1d(ipx,icy,icz)]; sync all
    vel(:,:,1,:)=vel(:,:,1,:)+vel(:,:,ng+1,:)[image1d(icx,iny,icz)]
    vel(:,:,ng,:)=vel(:,:,ng,:)+vel(:,:,0,:)[image1d(icx,ipy,icz)]; sync all
    vel(:,:,:,1)=vel(:,:,:,1)+vel(:,:,:,ng+1)[image1d(icx,icy,inz)]
    vel(:,:,:,ng)=vel(:,:,:,ng)+vel(:,:,:,0)[image1d(icx,icy,ipz)]; sync all
    vel_c=vel(:,1:ng,1:ng,1:ng)

    vel_c(1,:,:,:)=vel_c(1,:,:,:)/rho_c
    vel_c(2,:,:,:)=vel_c(2,:,:,:)/rho_c
    vel_c(3,:,:,:)=vel_c(3,:,:,:)/rho_c
    !print*, 'check: min,max,sum of rho0 = '
    !print*, minval(rho_c),maxval(rho_c),sum(rho_c*1d0)

    rho8=sum(rho_c*1d0); sync all
    ! co_sum
    if (head) then
      do i=2,nn**3
        rho8=rho8+rho8[i]
      enddo
      print*,'rho_global',rho8,ng_global
    endif; sync all
    rho8=rho8[1]; sync all
    ! convert to density contrast
    rho_c=rho_c/(rho8/ng_global/ng_global/ng_global)-1
    ! check normalization
    print*, minval(rho_c),maxval(rho_c),sum(rho_c*1d0)/ng/ng/ng; sync all
    if (head) print*,'Write delta_c into',output_name('delta_c')
    open(11,file=output_name('delta_c'),status='replace',access='stream')
    write(11) rho_c
    close(11); sync all



    ! neutrinos
    allocate(xp_nu(3,nplocal_nu),vp_nu(3,nplocal_nu))
    open(11,file=output_name('xp_nu'),status='old',action='read',access='stream')
    read(11) xp_nu; close(11)
    open(11,file=output_name('vp_nu'),status='old',action='read',access='stream')
    read(11) vp_nu; close(11)
    open(11,file=output_name('np_nu'),status='old',action='read',access='stream')
    read(11) rhoc; close(11)
    open(11,file=output_name('vc_nu'),status='old',action='read',access='stream')
    read(11) vfield; close(11)
    rho0=0; vel=0; nlast=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l
          pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp_nu(:,ip)+ishift_nu,izipx_nu)+rshift_nu)*x_resolution_nu
          pos1=pos1*real(ng)/real(nc) - 0.5
          idx1=floor(pos1)+1
          idx2=idx1+1
          dx1=idx1-pos1
          dx2=1-dx1
          ! density field
          rho0(idx1(1),idx1(2),idx1(3))=rho0(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p
          rho0(idx2(1),idx1(2),idx1(3))=rho0(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p
          rho0(idx1(1),idx2(2),idx1(3))=rho0(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p
          rho0(idx1(1),idx1(2),idx2(3))=rho0(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p
          rho0(idx1(1),idx2(2),idx2(3))=rho0(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p
          rho0(idx2(1),idx1(2),idx2(3))=rho0(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p
          rho0(idx2(1),idx2(2),idx1(3))=rho0(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p
          rho0(idx2(1),idx2(2),idx2(3))=rho0(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p
          ! velocity field
          vreal=tan((pi*real(vp_nu(:,ip)))/real(nvbin_nu-1)) / (sqrt(pi/2)/(sigma_vi_nu*vrel_boost))
          vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
          vel(:,idx1(1),idx1(2),idx1(3))=vel(:,idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*vreal
          vel(:,idx2(1),idx1(2),idx1(3))=vel(:,idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*vreal
          vel(:,idx1(1),idx2(2),idx1(3))=vel(:,idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*vreal
          vel(:,idx1(1),idx1(2),idx2(3))=vel(:,idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*vreal
          vel(:,idx1(1),idx2(2),idx2(3))=vel(:,idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*vreal
          vel(:,idx2(1),idx1(2),idx2(3))=vel(:,idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*vreal
          vel(:,idx2(1),idx2(2),idx1(3))=vel(:,idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*vreal
          vel(:,idx2(1),idx2(2),idx2(3))=vel(:,idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*vreal
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all
    deallocate(xp_nu,vp_nu)

    if (head) print*, 'Start sync from buffer regions'
    sync all
    rho0(1,:,:)=rho0(1,:,:)+rho0(ng+1,:,:)[image1d(inx,icy,icz)]
    rho0(ng,:,:)=rho0(ng,:,:)+rho0(0,:,:)[image1d(ipx,icy,icz)]; sync all
    rho0(:,1,:)=rho0(:,1,:)+rho0(:,ng+1,:)[image1d(icx,iny,icz)]
    rho0(:,ng,:)=rho0(:,ng,:)+rho0(:,0,:)[image1d(icx,ipy,icz)]; sync all
    rho0(:,:,1)=rho0(:,:,1)+rho0(:,:,ng+1)[image1d(icx,icy,inz)]
    rho0(:,:,ng)=rho0(:,:,ng)+rho0(:,:,0)[image1d(icx,icy,ipz)]; sync all
    rho_nu=rho0(1:ng,1:ng,1:ng)
    sync all
    vel(:,1,:,:)=vel(:,1,:,:)+vel(:,ng+1,:,:)[image1d(inx,icy,icz)]
    vel(:,ng,:,:)=vel(:,ng,:,:)+vel(:,0,:,:)[image1d(ipx,icy,icz)]; sync all
    vel(:,:,1,:)=vel(:,:,1,:)+vel(:,:,ng+1,:)[image1d(icx,iny,icz)]
    vel(:,:,ng,:)=vel(:,:,ng,:)+vel(:,:,0,:)[image1d(icx,ipy,icz)]; sync all
    vel(:,:,:,1)=vel(:,:,:,1)+vel(:,:,:,ng+1)[image1d(icx,icy,inz)]
    vel(:,:,:,ng)=vel(:,:,:,ng)+vel(:,:,:,0)[image1d(icx,icy,ipz)]; sync all
    vel_nu=vel(:,1:ng,1:ng,1:ng)

    vel_nu(1,:,:,:)=vel_nu(1,:,:,:)/rho_c
    vel_nu(2,:,:,:)=vel_nu(2,:,:,:)/rho_c
    vel_nu(3,:,:,:)=vel_nu(3,:,:,:)/rho_c
    !print*, 'check: min,max,sum of rho0 = '
    !print*, minval(rho_nu),maxval(rho_nu),sum(rho_nu*1d0)

    rho8=sum(rho_nu*1d0); sync all
    ! co_sum
    if (head) then
      do i=2,nn**3
        rho8=rho8+rho8[i]
      enddo
      print*,'rho_global',rho8,ng_global
    endif; sync all
    rho8=rho8[1]; sync all
    ! convert to density contrast
    rho_nu=rho_nu/(rho8/ng_global/ng_global/ng_global)-1
    ! check normalization
    print*, minval(rho_nu),maxval(rho_nu),sum(rho_nu*1d0)/ng/ng/ng; sync all

    if (head) print*,'Write delta_nu into',output_name('delta_nu')
    open(11,file=output_name('delta_nu'),status='replace',access='stream')
    write(11) rho_nu
    close(11); sync all

    ! power spectrum
    write(str_i,'(i6)') image
    write(str_z,'(f7.3)') z_checkpoint(cur_checkpoint)
    call cross_power(xi,rho_c,rho_nu)
    sync all
    if (head) then
      open(15,file=output_name('cicpower'),status='replace',access='stream')
      write(15) xi
      close(15)
    endif
    sync all

  enddo
  call destroy_penfft_plan
  sync all
  if (head) print*,'cic_velpower done'
  sync all
end

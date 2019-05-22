program convert_format
  implicit none

  integer,parameter :: np=256**3
  integer,parameter :: nc=128
  integer,parameter :: nnt=2
  integer,parameter :: nt=nc/nnt
  integer,parameter :: ncell=4
  integer,parameter :: nf=nc*ncell
  integer,parameter :: nft=nt*ncell
  real,parameter :: box=12.5 ! Mpc/h
  real(8),parameter :: x_resolution=1.0/(int(2,8)**(2*8))
  integer(8),parameter :: nvbin=int(2,8)**(8*2)
  real(8),parameter :: vrel_boost=2.5
  real,parameter :: pi=4*atan(1.)

  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt)
  integer(4) rholocal(nt,nt,nt,nnt,nnt,nnt)
  integer(8) cume(nt,nt,nt,nnt,nnt,nnt)
  real(4) vfield(3,nt,nt,nt,nnt,nnt,nnt)

  real xvp(6,np),xpg(3),xpt(3),xpc(3),vpg(3),vreal(3)
  integer(2) xp(3,np),vp(3,np)


  integer(4) itile(3),icg(3)
  integer(8) ip,idx

  type sim_header
    integer(8) nplocal,npglobal,nplocal_nu,npglobal_nu
    integer(8) izipx,izipv,izipx_nu,izipv_nu
    integer(8) image
    integer(8) nn,nnt,nt,ncell,ncb
    integer(8) timestep
    integer(8) cur_checkpoint
    integer(8) cur_halofind

    real a, t, tau
    real dt_pp, dt_fine, dt_coarse, dt_vmax, dt_vmax_nu
    real mass_p_cdm,mass_p_nu
    real box

    real h0
    real omega_m
    real omega_l
    real s8
    real vsim2phys
    real sigma_vres
    real sigma_vi
    real sigma_vi_nu
    real z_i,z_i_nu
    real vz_max
  endtype

  type(sim_header) sim


  print*,'convert format'
  print*,'nf =',nf

  sim%nplocal=np
  sim%nplocal_nu=0
  sim%npglobal=np
  sim%npglobal_nu=0
  sim%a=1./100
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
  sim%image=1
  sim%nn=1
  sim%nnt=nnt
  sim%nt=nt
  sim%ncell=ncell
  sim%ncb=6
  sim%izipx=2
  sim%izipv=2
  sim%izipx_nu=2
  sim%izipv_nu=2
  sim%h0=0.6727
  sim%omega_m=0.3156
  sim%omega_l=0.6844
  sim%s8=0.8
  sim%vsim2phys=(150./sim%a)*box*sqrt(sim%omega_m)/nf
  sim%z_i=99
  sim%z_i_nu=99
  sim%mass_p_cdm=real(nf**3)/real(np)
  print*,'sim%vsim2phys =',sim%vsim2phys

  open(11,file='IC.dat',status='old',action='read',access='stream')
  read(11) xvp
  close(11)
  xvp(1:3,:)=xvp(1:3,:)*nf/box/1000
  xvp(4:6,:)=xvp(4:6,:)/10 ! Gadget factor

  print*, maxval(xvp(:3,:))
  print*, minval(xvp(:3,:))


  sim%sigma_vi=sqrt(sum((xvp(4:6,:)*1d0/sim%vsim2phys)**2)/np)
  print*,'sim%sigma_vi =',sim%sigma_vi
  rhoc=0
  vfield=0
  do ip=1,np
    xpg=modulo(xvp(1:3,ip),real(nf))
    vpg=xvp(4:6,ip)/sim%vsim2phys
    itile=xpg/nft+1   ! print*,itile
    xpt=xpg-nft*(itile-1) ! print*,xpt
    icg=xpt/ncell+1  ! print*,icg
    rhoc(icg(1),icg(2),icg(3),itile(1),itile(2),itile(3))=rhoc(icg(1),icg(2),icg(3),itile(1),itile(2),itile(3))+1
    vfield(:,icg(1),icg(2),icg(3),itile(1),itile(2),itile(3))=vfield(:,icg(1),icg(2),icg(3),itile(1),itile(2),itile(3))+vpg
  enddo
  print*, 'sum(rhoc), maxval(rhoc), minval(rhoc)'
  print*, sum(rhoc), maxval(rhoc), minval(rhoc)
  vfield(1,:,:,:,:,:,:)=vfield(1,:,:,:,:,:,:)/rhoc
  vfield(2,:,:,:,:,:,:)=vfield(2,:,:,:,:,:,:)/rhoc
  vfield(3,:,:,:,:,:,:)=vfield(3,:,:,:,:,:,:)/rhoc

  cume=cumsum6(rhoc)

  rholocal=0
  xp=0
  vp=0
  do ip=1,np
    xpg=modulo(xvp(1:3,ip),real(nf))
    vpg=xvp(4:6,ip)/sim%vsim2phys
    itile=xpg/nft+1 !; print*,itile
    xpt=xpg-nft*(itile-1) !; print*,xpt
    icg=xpt/ncell+1 !; print*,icg
    xpc=xpt-ncell*(icg-1) !; print*,xpc
    rholocal(icg(1),icg(2),icg(3),itile(1),itile(2),itile(3))=rholocal(icg(1),icg(2),icg(3),itile(1),itile(2),itile(3))+1
    idx=cume(icg(1),icg(2),icg(3),itile(1),itile(2),itile(3)) &
       -rhoc(icg(1),icg(2),icg(3),itile(1),itile(2),itile(3)) &
       +rholocal(icg(1),icg(2),icg(3),itile(1),itile(2),itile(3));

    xp(:,idx)=floor(xpc/x_resolution,kind=8) !; print*, xp(:,idx); stop
    vreal=vpg-vfield(:,icg(1),icg(2),icg(3),itile(1),itile(2),itile(3))
    vp(:,idx)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sim%sigma_vi*vrel_boost)*vreal)/pi,kind=2) !; print*, vp(:,idx); stop
  enddo


  open(11,file='../../output/universe200/image1/99.000_info_1.bin',status='replace',access='stream')
  write(11) sim
  close(11)

  open(11,file='../../output/universe200/image1/99.000_xp_1.bin',status='replace',access='stream')
  open(12,file='../../output/universe200/image1/99.000_vp_1.bin',status='replace',access='stream')
  open(13,file='../../output/universe200/image1/99.000_np_1.bin',status='replace',access='stream')
  open(14,file='../../output/universe200/image1/99.000_vc_1.bin',status='replace',access='stream')

  write(11) xp
  write(12) vp
  write(13) rhoc
  write(14) vfield

  close(11)
  close(12)
  close(13)
  close(14)

contains

  function cumsum6(input)
      implicit none
      integer(4) input(nt,nt,nt,nnt,nnt,nnt)
      integer(8) cumsum6(nt,nt,nt,nnt,nnt,nnt)
      integer(8) nsum,i,j,k,itx,ity,itz
      nsum=0
      do itz=1,nnt
      do ity=1,nnt
      do itx=1,nnt
        do k=1,nt
        do j=1,nt
        do i=1,nt
          nsum=nsum+input(i,j,k,itx,ity,itz)
          cumsum6(i,j,k,itx,ity,itz)=nsum
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
  endfunction
end

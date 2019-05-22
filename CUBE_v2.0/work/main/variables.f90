module variables
  use omp_lib
  use parameters
  implicit none
  save

  ! parameters
  integer(8),parameter :: np_image=(nc*np_nc)**3*merge(2,1,body_centered_cubic) ! average number of particles per image
  integer(8),parameter :: np_image_max=np_image*(nte*1./nt)**3*image_buffer
  integer(8),parameter :: np_tile_max=np_image/nnt**3*(nte*1./nt)**3*tile_buffer
  integer(8),parameter ::  nseedmax=200
  integer(8),parameter :: unit8=1
  real,parameter :: dt_max=1
  real,parameter :: dt_scale=1
  real,parameter :: GG=1.0/6.0/pi

  logical neutrino_flag
  integer n_checkpoint_neu ! the checkpoint to switch on neutrinos

  ! variables
  integer(8) istep
  real dt[*],dt_old[*],dt_mid[*],dt_e
  real a[*],da[*],a_mid[*],tau[*],t[*] ! time step
  real f2_max_fine(nnt,nnt,nnt)[*],f2_max_pp[*],f2_max_coarse[*]

  integer(4) iseed(nseedmax), iseedsize, nth,ith
  integer(8) itx,ity,itz,ix,iy,iz,i_dim
  integer(8) i,j,k,l,ip,ipp,pp,nzero
  integer(8) nptile(nnt,nnt,nnt), npcheck
  real(8) xq(3),deltax(3),deltav(3),vreal(3)

  ! FFT plans
  integer stat
  integer(8) plan_fft_fine,plan_ifft_fine
  real vmax(3),vmax_nu(3),overhead_tile[*],overhead_image[*]
  real sigma_vi,sigma_vi_new,sigma_vi_nu,sigma_vi_new_nu
  !real vdisp(506,2),sigma_vi_old,sigma_vi
  real(4) svz(500,2),svr(100,2)
  real(8) sigma_vci,sigma_vfi,sigma_vres,sigma_vci_old,sigma_vfi_old,sigma_vres_old
  real(8) sum_c,sum_r,sum_s
  real(8) std_vsim_c[*],std_vsim_res[*],std_vsim[*]
  real(8) std_vsim_c_nu[*],std_vsim_res_nu[*],std_vsim_nu[*]
  ! n^3
  integer(izipx) xp(3,np_image_max)[*], xp_new(3,np_tile_max)
  integer(izipv) vp(3,np_image_max)[*], vp_new(3,np_tile_max)
#ifdef PID
    integer(4) pid(np_image_max)[*], pid_new(np_tile_max)
#endif
#ifdef AFIELD
    real(4) afield(3,np_image_max)
#endif

  real rho_f(nfe+2,nfe,nfe)
  real crho_f(nfe+2,nfe,nfe)
  real kern_f(nfe/2+1,nfe,nfe,3)
  real force_f(3,nfb:nfe-nfb+1,nfb:nfe-nfb+1,nfb:nfe-nfb+1)

  integer(4) rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
  real(4) vfield(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt) ! cannot have >7 dims

  integer(8),dimension(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt),codimension[*] :: idx_b_l,idx_b_r
  integer(8),dimension(nt,nt,nnt,nnt,nnt),codimension[*] :: ppl0,pplr,pprl,ppr0,ppl,ppr

  ! the following variables are introduced because
  ! gcc only allows <= 7 ranks in arrays
  real(4) vtransx(3,ncb,nt+2*ncb,nt+2*ncb,nnt,nnt)[*]
  real(4) vtransy(3,nt+2*ncb,ncb,nt+2*ncb,nnt,nnt)[*]
  real(4) vtransz(3,nt+2*ncb,nt+2*ncb,ncb,nnt,nnt)[*]

  ! coarse kernel arrays
  real ck(3,nc,nc,nc)
  real kern_c(nc*nn/2+1,nc,npen,3)
  real crho_c(nc*nn+2,nc,npen) !!! temp
  real force_c(3,0:nc+1,0:nc+1,0:nc+1)[*]

  character (10) :: img_s, z_s

  !equivalence(rhoce,rhoce1d)
  integer(8),parameter :: np_pp_max=np_tile_max ! can be optimized lower
  integer(8) npairs,itest1,nlast
  real(8) xvec1(3),xvec2(3),xvec21(3),rmag,force_pp(3),rcut,pcut,f_tot(3)
  integer(4) ivec1(3),np,ii,jj,kk,np1,np2,l1,l2
  integer(4) igx,igy,igz,lp,ip1,ip2
  integer(4) t1,t2,tt1,tt2,ttt1,ttt2,t_rate

  type type_header_halo_tab
    integer(4) Ngroups,TotNgroups,Nids
    integer(8) TotNids
    integer(4) NFiles
  endtype
  type(type_header_halo_tab) header_halo_tab[*]
  integer nhalo[*],nhalo_tot[*],n_peak[*],n_peak_real[*],n_search_fail[*]

contains

  subroutine spine_tile(rhoce,idx_ex_r,pp_l,pp_r,ppe_l,ppe_r)
    !! make a particle index (cumulative sumation) on tile
    !! used in update_particle, initial_conditions
    !! input:
    !! rhoce -- particle number density on tile, with 2x buffer depth
    !! output:
    !! idx_ex_r -- last extended index on extended right boundary
    !! ppe_r -- last extended index on physical right boundary
    !! pp_l -- first physical index on physical left boundary
    !! pp_r -- last physical index on physical right boundary
    !! ppe_l -- first extended index on physical left boundary
    use omp_lib
    implicit none
    integer(4),intent(in) :: rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
    integer(8),dimension(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb),intent(out) :: idx_ex_r
    integer(8),dimension(nt,nt),intent(out) :: pp_l,pp_r,ppe_l,ppe_r
    integer(8) nsum,np_phy
    integer igz,igy
    ! spine := yz-plane to record cumulative sum
    nsum=0
    do igz=1-2*ncb,nt+2*ncb
    do igy=1-2*ncb,nt+2*ncb
      nsum=nsum+sum(rhoce(:,igy,igz))
      idx_ex_r(igy,igz)=nsum ! right index
    enddo
    enddo
    !$omp paralleldo default(shared) private(igz,igy)
    do igz=1,nt
    do igy=1,nt
      ppe_r(igy,igz)=idx_ex_r(igy,igz)-sum(rhoce(nt+1:,igy,igz))
    enddo
    enddo
    !$omp endparalleldo
    nsum=0
    do igz=1,nt
    do igy=1,nt
      pp_l(igy,igz)=nsum+1
      np_phy=sum(rhoce(1:nt,igy,igz))
      nsum=nsum+np_phy
      pp_r(igy,igz)=nsum
      ppe_l(igy,igz)=ppe_r(igy,igz)-np_phy+1
    enddo
    enddo
  endsubroutine

  subroutine spine_image(rhoc,idx_b_l,idx_b_r,ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr)
    !! make a particle index (cumulative sumation) on image
    !! used in buffer_grid
    !! input:
    !! rhoc -- particle number density on image, with 1x buffer depth
    !! output:
    !! idx_b_l -- zeroth extended index on extended left boundary
    !! idx_b_r -- last extended index on extended right boundary
    !! ppl -- zeroth physical index on physical left boundary
    !! ppr -- last physical index on physical right boundary
    !! ppe_r0 -- last extended index on physical right boundary
    !! ppe_l0 -- zeroth extended index on physical left boundary
    !! ppe_rl -- last extended index on inner right boundary
    !! ppe_lr -- zeroth extended index on inner left boundary
    implicit none
    integer(4),intent(in) :: rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)
    integer(8),dimension(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt),intent(out) :: idx_b_l,idx_b_r
    integer(8),dimension(nt,nt,nnt,nnt,nnt),intent(out) :: ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr
    integer(8) nsum,nsum_p,np_phy,ihz,ihy,ihx,igz,igy,ctile_mass(nnt,nnt,nnt),ctile_mass_p(nnt,nnt,nnt)
    ! spine := yz-plane to record cumulative sum
    nsum=0;nsum_p=0
    do ihz=1,nnt ! sum cumulative tile mass first
    do ihy=1,nnt
    do ihx=1,nnt
      ctile_mass(ihx,ihy,ihz)=nsum
      ctile_mass_p(ihx,ihy,ihz)=nsum_p
      nsum=nsum+sum(rhoc(:,:,:,ihx,ihy,ihz))
      nsum_p=nsum_p+sum(rhoc(1:nt,1:nt,1:nt,ihx,ihy,ihz))
    enddo
    enddo
    enddo
    !$omp paralleldo default(shared) private(ihz,ihy,ihx,nsum,igz,igy)
    do ihz=1,nnt ! calculate extended spine cumulative index on both sides
    do ihy=1,nnt
    do ihx=1,nnt
      nsum=ctile_mass(ihx,ihy,ihz)
      do igz=1-ncb,nt+ncb
      do igy=1-ncb,nt+ncb
        idx_b_l(igy,igz,ihx,ihy,ihz)=nsum
        nsum=nsum+sum(rhoc(:,igy,igz,ihx,ihy,ihz))
        idx_b_r(igy,igz,ihx,ihy,ihz)=nsum
      enddo
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
    !$omp paralleldo default(shared) private(ihz,ihy,ihx,nsum_p,igz,igy,np_phy)
    do ihz=1,nnt ! calculate physical spine
    do ihy=1,nnt
    do ihx=1,nnt
      nsum_p=ctile_mass_p(ihx,ihy,ihz)
      do igz=1,nt
      do igy=1,nt
        ppl(igy,igz,ihx,ihy,ihz)=nsum_p
        np_phy=sum(rhoc(1:nt,igy,igz,ihx,ihy,ihz))
        nsum_p=nsum_p+np_phy
        ppr(igy,igz,ihx,ihy,ihz)=nsum_p

        ppe_r0(igy,igz,ihx,ihy,ihz)=idx_b_r(igy,igz,ihx,ihy,ihz)-sum(rhoc(nt+1:,igy,igz,ihx,ihy,ihz))
        ppe_rl(igy,igz,ihx,ihy,ihz)= ppe_r0(igy,igz,ihx,ihy,ihz)-sum(rhoc(nt-ncb+1:nt,igy,igz,ihx,ihy,ihz))

        ppe_l0(igy,igz,ihx,ihy,ihz)=idx_b_l(igy,igz,ihx,ihy,ihz)+sum(rhoc(:0,igy,igz,ihx,ihy,ihz))
        ppe_lr(igy,igz,ihx,ihy,ihz)=ppe_l0(igy,igz,ihx,ihy,ihz) +sum(rhoc(1:ncb,igy,igz,ihx,ihy,ihz))
      enddo
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
  endsubroutine

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
    print*,term_z,term_r
  endfunction

endmodule

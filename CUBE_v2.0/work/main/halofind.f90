#define HID
module halo_output
  implicit none

  type type_halo_catalog
    integer nhalo_tot,nhalo
    real den_odc
  endtype

  type type_halo_info
    real hpos(3)
    real mass_odc,radius_odc,v_disp
    real x_mean(3),v_mean(3),ang_mom(3),var_x(3),inertia(3,3)
    real q_mean(3),inertia_q(3,3),s_mean(3)
  endtype
endmodule

subroutine halofind
  use omp_lib
  use variables
  use halo_output
  implicit none
  save

  type(type_halo_info) halo_info

  integer,parameter :: nc_halo_max=128
  integer,parameter :: newbox_max=300
  integer,parameter :: nlist=5*(nc_halo_max+1)**3
  integer,parameter :: n_peak_max=5*nc_halo_max**3
  integer,parameter :: max_halo_np=5*(nc_halo_max+1)**3/3
  integer, parameter :: search_ratio = 2
  integer, parameter :: refine_ratio = 3
  logical,parameter :: NGP = .false. ! NGP mass assignment by default
  logical, parameter :: odc_vir = .false. ! use virial density as threshold
  logical physical_halo
  real,parameter :: den_peak_cutoff=100 ! lower density peak threshold for determining which peaks should be inspected as halos
  integer,parameter :: min_halo_particles=100 ! lower mass cutoff for cataloging a halo, in particles

  integer i0,j0,k0,l0,idx1(3),idx2(3),irtot
  integer(1) p_tag(np_image_max)
  integer(4) ilist_odc(max_halo_np),i_odc,GroupOffset,temp_halo(2,200000)
  integer idist(3,nlist),isortdist(nlist),idist_tmp(3,nlist),isortpeak(n_peak_max)
  integer isortpos_odc(max_halo_np)
  integer iloc,np_odc,np_search,nptemp,qtemp
  integer crbox(3,2),csbox(3,2),frbox(3,2),itile(3),newbox(3),imax(3)
  real den_odc,d1,d2,r1,r2,w1,w2,den1,den2,dx(3),dv(3),qpos(3),dx_mean(3),qpos_mean(3),inert(3,3)
  real finegrid(0:newbox_max+1,0:newbox_max+1,0:newbox_max+1)
  real xv_odc(6,max_halo_np),r_odc(max_halo_np)
  integer(4) oct(max_halo_np)
  real dx1(3),dx2(3),pos1(3),rr,rdist(nlist),xflat,den_peak(n_peak_max),ipeak(3,n_peak_max),m_mesh,denmax
  real hpos(3),mass_proxy,newgrid,r_proxy,rsearch,dr(3)
  real rhof(1-nfb:nft+nfb,1-nfb:nft+nfb,1-nfb:nft+nfb),halo_mesh_mass(n_peak_max)
#ifdef HID
  integer(2) hid(nf**3)
#endif

  real,parameter :: cen(3)=[30,370,370]
  real(4) los(3),sx(3)

  if (head) print*, ''
  if (head) print*, 'halofind'

  ! initialize_halofinder
  ! Loop through a box of length 2*nc_halo_max
  ! if cell is within sphere of radius = box length / 2
  ! include distince in rdist at entry ii
  ! ordered bottom left to top right
  print*, 'initialize_halofinder'
  ii=0
  do i0=-nc_halo_max,nc_halo_max
  do j0=-nc_halo_max,nc_halo_max
  do k0=-nc_halo_max,nc_halo_max
    !rr=norm2((/real(i0),real(j0),real(k0)/))
    rr=norm2([real :: i0,j0,k0])
    if (rr>nc_halo_max) cycle
    ii=ii+1
    idist(:,ii)=[i0,j0,k0]
    rdist(ii)=rr
  enddo
  enddo
  enddo
  irtot=ii
  if (irtot>nlist) stop 'irtot>nlist'
  ! sorts the rdist array from lowest to highest radial position
  ! from center of sphere saves rdist array position values in idist
  isortdist(:irtot)=(/ (i0,i0=1,irtot) /)
  call indexedsort(irtot,rdist,isortdist)
  idist(:,:irtot)=idist(:,isortdist(:irtot))


  ! Find halo candidates based on local overdensities for each tile
  ! Determine which candidates are to be considered halos and write their properties to file.
  if (odc_vir) then
    ! Determine delta_vir using equation (6) of Bryan et al. (1997) for a flat universe
    print*, 'scale_factor =',1./(1+z_halofind(cur_halofind))
    a=1./(1+z_halofind(cur_halofind))
    xflat=(omega_m/a**3)/(omega_m/a**3+omega_l)-1.
    den_odc=18.*pi**2+82.*xflat-39.*xflat**2
  else
    den_odc=200
  endif

  if (head) then
    print*, '  mass_p       =',sim%mass_p_cdm
    print*, '  nplocal      =',sim%nplocal
    print*, '  nplocal_nu   =',sim%nplocal_nu
    print*, '  np_image_max =',np_image_max
    print*, '  scale_factor,den_odc =',a,den_odc
  endif

  ! open halo file
  open(11,file=output_name_halo('halo'),status='replace',access='stream')
  open(12,file=output_name_halo('halo_pid'),status='replace',access='stream')
  !open(101,file=output_name_halo('group_tab'),status='replace',access='stream')
  !open(102,file=output_name_halo('group_ids'),status='replace',access='stream')
  ! Determine which candidates are to be considered halos and write their properties to file.
# ifdef HID
    open(21,file=output_name_halo('hid'),status='replace',access='stream')
    hid=0 !! initialize hid so that all particles do not belong to any halo
# endif
  nhalo=0
  GroupOffset=0
  n_search_fail=0 ! Ticker recording how many particle searches end up empty handed
  p_tag=0 ! Initialize so that no particles are yet part of a halo
  write(11) nhalo_tot,nhalo,den_odc
  ! write(12) ! there is no header for halo_pid
  header_halo_tab%Ngroups=0
  header_halo_tab%TotNgroups=0
  header_halo_tab%Nids=0
  header_halo_tab%TotNids=0
  !write(101) header_halo_tab
  !write(102) header_halo_tab,0

  do itz=nnt,1,-1
  do ity=nnt,1,-1
  do itx=nnt,1,-1

    if (head) print*, '  fine mass assignment on tile',int(itx,1),int(ity,1),int(itz,1)
    !call find_halo_candidates(tile, n_peak)
    !subroutine find_halo_candidates(tile, ic)
    rhof=0; isortpeak=0; den_peak=0; n_peak=0; n_peak_real=0;
    !$omp paralleldo default(shared) &
    !$omp& private(k0,j0,i0,np,nzero,l0,ip,pos1,idx1)
    do k0=2-ncb,nt+ncb-1
    do j0=2-ncb,nt+ncb-1
    do i0=2-ncb,nt+ncb-1
      np=rhoc(i0,j0,k0,itx,ity,itz)
      nzero=idx_b_r(j0,k0,itx,ity,itz)-sum(rhoc(i0:,j0,k0,itx,ity,itz))
      do l0=1,np
        ip=nzero+l0
        if (NGP) then
          pos1=ncell*([i0,j0,k0]-1) + ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          idx1=floor(pos1)+1
          rhof(idx1(1),idx1(2),idx1(3))=rhof(idx1(1),idx1(2),idx1(3))+sim%mass_p_cdm
        else
          pos1=ncell*([i0,j0,k0]-1) + ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution - 0.5
          idx1=floor(pos1)+1
          idx2=idx1+1
          dx1=idx1-pos1
          dx2=1-dx1
          rhof(idx1(1),idx1(2),idx1(3))=rhof(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*sim%mass_p_cdm
          rhof(idx2(1),idx1(2),idx1(3))=rhof(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*sim%mass_p_cdm
          rhof(idx1(1),idx2(2),idx1(3))=rhof(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*sim%mass_p_cdm
          rhof(idx1(1),idx1(2),idx2(3))=rhof(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*sim%mass_p_cdm
          rhof(idx1(1),idx2(2),idx2(3))=rhof(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*sim%mass_p_cdm
          rhof(idx2(1),idx1(2),idx2(3))=rhof(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*sim%mass_p_cdm
          rhof(idx2(1),idx2(2),idx1(3))=rhof(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*sim%mass_p_cdm
          rhof(idx2(1),idx2(2),idx2(3))=rhof(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*sim%mass_p_cdm
        endif
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo

    if (head) print*,'  find density peaks'
    do k0=1-ncell,nft+ncell  ! 1 more coarse cell layer
    do j0=1-ncell,nft+ncell
    do i0=1-ncell,nft+ncell
      denmax=maxval(rhof(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1))
      if (denmax==rhof(i0,j0,k0).and.denmax>den_peak_cutoff) then
        !print*,maxloc(rhof(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1))
        if (n_peak>n_peak_max-2) stop 'too many peaks'
        ! Find the fine mesh mass of this peak
        m_mesh=0
        do ii=1,irtot ! keep adding mass until mean density is less than den_odc
          ix=i0+idist(1,ii); if(ix<=1-nfb .or. ix>=nft+nfb-1) cycle
          iy=j0+idist(2,ii); if(iy<=1-nfb .or. iy>=nft+nfb-1) cycle
          iz=k0+idist(3,ii); if(iz<=1-nfb .or. iz>=nft+nfb-1) cycle
          ! skip the outerest layer, due to CIC mass assignment
          m_mesh=m_mesh+rhof(ix,iy,iz)
          if (rdist(ii+1)==rdist(ii)) cycle ! keep completing a shell
          if (ii>=19 .and. m_mesh/ii<den_odc) exit ! for low-mass peaks go upto sqrt(2)
        enddo
        ! Consider this a halo candidate if the fine mesh mass is large enough
        if (m_mesh>sim%mass_p_cdm*min_halo_particles/2.) then
          n_peak=n_peak+1
          n_peak_real=n_peak_real+merge(1,0,minval([i0,j0,k0])>=1 .and. maxval([i0,j0,k0])<=nft)
          ! NGP maximum
          !ipeak(:,n_peak)=[i0,j0,k0]-0.5
          ! parabolic interpolation maximum
          ipeak(1,n_peak)=(i0-0.5)+0.5*(rhof(i0+1,j0,k0)-rhof(i0-1,j0,k0))&
          /(-rhof(i0-1,j0,k0)+2*rhof(i0,j0,k0)-rhof(i0+1,j0,k0))
          ipeak(2,n_peak)=(j0-0.5)+0.5*(rhof(i0,j0+1,k0)-rhof(i0,j0-1,k0))&
          /(-rhof(i0,j0-1,k0)+2*rhof(i0,j0,k0)-rhof(i0,j0+1,k0))
          ipeak(3,n_peak)=(k0-0.5)+0.5*(rhof(i0,j0,k0+1)-rhof(i0,j0,k0-1))&
          /(-rhof(i0,j0,k0-1)+2*rhof(i0,j0,k0)-rhof(i0,j0,k0+1))

          den_peak(n_peak)=denmax
          halo_mesh_mass(n_peak)=m_mesh
        endif
      endif
    enddo
    enddo
    enddo
    if (head) print*, '  real/total peaks =',n_peak_real,n_peak
    ! sort density maxima
    isortpeak(:n_peak) = [(i0,i0=1,n_peak)]
    call indexedsort(n_peak,den_peak(:),isortpeak(:))
    ipeak(:,:n_peak)=ipeak(:,isortpeak(:n_peak))
    halo_mesh_mass(:n_peak)=halo_mesh_mass(isortpeak(:n_peak))

    !print*, n_peak
    !print*, 'mass_proxy'
    !print*,halo_mesh_mass(1:n_peak);
    !print*, 'denmax'
    !print*,den_peak(1:n_peak); stop

    if (head) print*,'  find halo particles'
    do iloc=n_peak,1,-1
      ! determine searching regions
      mass_proxy=halo_mesh_mass(iloc)
      hpos=ipeak(:,iloc)

      ! find_halo_particles
      ! Search for particles by looking at local particle distribution
      !call find_halo_particles(halo_vir, mass_proxy, hpos(:), r_odc, i_odc, 1)
      !call find_halo_particles(den_odc, mass_proxy, hpos(:), r_odc, i_odc, 2)

      ! Use the mass proxy to guess at the size of the halo
      ! the radius within which the refined density peak will be found
      r_proxy=(0.75/pi/den_odc*mass_proxy)**(1./3.)
      ! the (larger) radius to store particle positions determine which are part of the halo
      rsearch = min(search_ratio*r_proxy,real(ncb-1)*ncell)
      itile=[itx,ity,itz]
#ifdef analysis
      print*, 'den_peak =',den_peak(iloc)
      print*, 'mass_proxy =',mass_proxy
      print*, 'r_proxy =',r_proxy
      print*, 'rsearch =',rsearch
      print*, 'rho max at:',([itx,ity,itz]-1)*nft+hpos
#endif
      ! larger box for identifying particles
      csbox(:,1)=floor((hpos-rsearch)/ncell)+1
      csbox(:,2)=floor((hpos+rsearch)/ncell)+1
      if (minval(csbox)<1-ncb) then
        print*,csbox
        stop
      endif
      if (maxval(csbox)>nt+ncb) then
        print*,csbox
        stop
      endif
      ! smaller box for finding center
      crbox(:,1)=floor((hpos-r_proxy)/ncell)+1
      crbox(:,2)=floor((hpos+r_proxy)/ncell)+1
      frbox(:,1)=ncell*(crbox(:,1) - 1) ! Boundary of the refined region in fine mesh cell units, integer coordinates
      frbox(:,2)=ncell*crbox(:,2)
      newbox(1:3)=refine_ratio*(frbox(:,2)-frbox(:,1)) ! Number of extra refined cells in this region
      newgrid    = 1./refine_ratio ! and their spacing
      if (maxval(newbox)>newbox_max) stop 'newbox>newbox_max'

      ! find halo particles
      np_search=0;np_odc=0;finegrid=0;xv_odc=0
#     ifdef analysis
        open(23,file=output_name_halo('xvp_halo'),access='stream',status='replace')
#     ifdef PID
        open(24,file=output_name_halo('pid_halo'),access='stream',status='replace')
#     endif
        write(23) halo_info%hpos*0, halo_info%x_mean*0, halo_info%radius_odc*0, halo_info%mass_odc*0, i_odc*0
        write(23) hpos+([itx,ity,itz]-1)*nft
#     endif
      do k0=csbox(3,1),csbox(3,2)
      do j0=csbox(2,1),csbox(2,2)
      do i0=csbox(1,1),csbox(1,2)
        np=rhoc(i0,j0,k0,itx,ity,itz)
        nzero=idx_b_r(j0,k0,itx,ity,itz)-sum(rhoc(i0:,j0,k0,itx,ity,itz))
        do l0=1,np
          ip=nzero+l0
          pos1=ncell*([i0,j0,k0]-1)+ncell*(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          dr=pos1-hpos
          rr=norm2(dr)
          if (rr < rsearch) then
            if (p_tag(ip)==0) then
              np_odc=np_odc+1
              ilist_odc(np_odc)=ip
              xv_odc(1:3,np_odc)=pos1
              xv_odc(4:6,np_odc)=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost)) &
                                 +vfield(:,i0,j0,k0,itx,ity,itz)
              if (rr < r_proxy) then
                if (NGP) then
                  idx1=floor((pos1-frbox(:,1))/newgrid)+1
                  finegrid(idx1(1),idx1(2),idx1(3))=finegrid(idx1(1),idx1(2),idx1(3))+1
                else
                  pos1=pos1-0.5*newgrid
                  idx1=floor((pos1-frbox(:,1))/newgrid)+1
                  idx2=idx1+1
                  dx1=idx1-(pos1-frbox(:,1))/newgrid
                  dx2=1-dx1
!                  if (minval(idx1)<0 .or. minval(idx2)<0) then
!                    print*,'r',r_proxy,rsearch
!                    print*,'frbox',frbox
!                    print*,'csbox',csbox
!                    print*,'crbox',crbox
!                    stop
!                  endif
                  !print*, idx1
                  !print*, dx1
                  !stop
                  finegrid(idx1(1),idx1(2),idx1(3))=finegrid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
                  finegrid(idx2(1),idx1(2),idx1(3))=finegrid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
                  finegrid(idx1(1),idx2(2),idx1(3))=finegrid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
                  finegrid(idx1(1),idx1(2),idx2(3))=finegrid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
                  finegrid(idx1(1),idx2(2),idx2(3))=finegrid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
                  finegrid(idx2(1),idx1(2),idx2(3))=finegrid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
                  finegrid(idx2(1),idx2(2),idx1(3))=finegrid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
                  finegrid(idx2(1),idx2(2),idx2(3))=finegrid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)
                endif
              endif
            endif
          endif
        enddo
      enddo
      enddo
      enddo
      if(np_odc>max_halo_np) stop 'np_search>max_halo_np'
      ! Find refined mesh density maximum, w.r.t. tile
      !hpos=frbox(:,1)+newgrid*(maxloc(finegrid)-0.5) !
      !parabolic
        !{ debug
        !finegrid(80:82,80:82,80:82)=100000
        !finegrid(81:82,81:82,81:82)=110000
        !finegrid(81,81,81)=170500
        ! debug}
      !imax=maxloc(finegrid)-1
      imax=maxloc(finegrid(1:newbox_max,1:newbox_max,1:newbox_max))

      !print*,'imax',imax
!if(i0-1==-1) then
!  print*, 'maxloc',maxloc(finegrid)
!  stop
!endif
      i0=imax(1);j0=imax(2);k0=imax(3)
      hpos(1)=(i0-0.5)+0.5*(finegrid(i0+1,j0,k0)-finegrid(i0-1,j0,k0))&
                       /(-finegrid(i0-1,j0,k0)+2*finegrid(i0,j0,k0)-finegrid(i0+1,j0,k0))
      hpos(2)=(j0-0.5)+0.5*(finegrid(i0,j0+1,k0)-finegrid(i0,j0-1,k0))&
                       /(-finegrid(i0,j0-1,k0)+2*finegrid(i0,j0,k0)-finegrid(i0,j0+1,k0))
      hpos(3)=(k0-0.5)+0.5*(finegrid(i0,j0,k0+1)-finegrid(i0,j0,k0-1))&
                       /(-finegrid(i0,j0,k0-1)+2*finegrid(i0,j0,k0)-finegrid(i0,j0,k0-1))
      !print*,'hpos in finegrid',hpos
#     ifdef analysis
        open(55,file=output_name_halo('refine'),access='stream',status='replace')
        write(55) hpos,finegrid
        close(55)
#     endif
      hpos=frbox(:,1)+newgrid*hpos

      physical_halo=(minval(hpos)>0 .and. maxval(hpos)<nft)
      ! sort by radius -------------------------------------------------------
      r_odc(:np_odc)=norm2(xv_odc(:3,:np_odc)-spread(hpos,2,np_odc),1)
      oct(:np_odc)=1+sum(spread([1,2,4],2,np_odc)*merge(0,1,xv_odc(:3,:np_odc)-spread(hpos,2,np_odc)>0),1)
      isortpos_odc(:np_odc)=[(i0,i0=1,np_odc)]
      call indexedsort(np_odc,r_odc,isortpos_odc)
      xv_odc(:,:np_odc)=xv_odc(:,isortpos_odc(:np_odc))
      ilist_odc(:np_odc)=ilist_odc(isortpos_odc(:np_odc))
      oct(:np_odc)=oct(isortpos_odc(:np_odc))
#ifdef analysis
  xv_odc(:3,:np_odc)=xv_odc(:3,:np_odc)+spread(([itx,ity,itz]-1)*nft, dim=2, ncopies=np_odc)
  write(23) xv_odc(:,:np_odc)
  xv_odc(:3,:np_odc)=xv_odc(:3,:np_odc)-spread(([itx,ity,itz]-1)*nft, dim=2, ncopies=np_odc)
#ifdef PID
  write(24) pid(ilist_odc(:np_odc))
#endif
  !close(23)
  open(33,file=output_name_halo('r_odc'),access='stream',status='replace')
  write(33) r_odc(:np_odc)
  close(33)
  open(33,file=output_name_halo('oct_odc'),access='stream',status='replace')
  write(33) int(oct(:np_odc),1)
  close(33)
  print*, 'wrote',np_odc,'particles into file'
#endif
      ! calculate halo radius ------------------------------------------------
      i_odc=0
      do ii=2,np_odc
        den2=0.75/pi*ii*sim%mass_p_cdm/r_odc(ii)**3
        if (den2<den_odc) then
          den1=0.75/pi*(ii-1)*sim%mass_p_cdm/r_odc(ii-1)**3
          r2 = log10(r_odc(ii))
          r1 = log10(r_odc(ii-1))
          d2 = log10(den2)
          d1 = log10(den1)
          w1 = log10(den_odc) - d2
          w2 = d1 - log10(den_odc)
          halo_info%radius_odc = 10**((w1 * r1 + w2 * r2) / (d1 - d2))
          i_odc=ii-1
          exit
        endif
      enddo

      !if (i_odc==0) then ! check number of particles
      !  n_search_fail=n_search_fail+1
      !  print*,'  search fail:'
      !  print*, ipeak(:,iloc)+([itx,ity,itz]-1)*nft
      !  print*,'  r',r_proxy,rsearch
      !  !stop
      !endif

      if (i_odc >= min_halo_particles) then
        p_tag(ilist_odc(1:i_odc))=1
        if (physical_halo) then  ! if the maximum is in physical regions
#         ifdef analysis
            print*,'  it is a physical halo'
#         endif
          nhalo=nhalo+1
          halo_info%hpos=modulo(hpos+([itx,ity,itz]-1)*nft,real(nf_global))
          halo_info%mass_odc=real(i_odc)
          halo_info%x_mean=sum(xv_odc(1:3,:i_odc),2)/i_odc
          halo_info%var_x=sum((xv_odc(1:3,:i_odc)-spread(halo_info%x_mean,2,i_odc))**2,2)/(i_odc-1)
          halo_info%x_mean=modulo(halo_info%x_mean+([itx,ity,itz]-1)*nft,real(nf_global))
          halo_info%v_mean=sum(xv_odc(4:6,:i_odc),2)/i_odc
          halo_info%v_disp=sqrt(sum((xv_odc(4:6,:i_odc)-spread(halo_info%v_mean,2,i_odc))**2)/(i_odc-1))

          los=halo_info%x_mean-(cen/500.)*ng_global
          los=los/norm2(los)
          sx=sum(halo_info%v_mean*los)*los ! project velocity to line of sight
          sx=sx*sim%vsim2phys/sim%a/(100*h0) ! convert to km/h and multiply 1/aH, in Mpc
          sx=sx/(h0*box/nf_global) ! convert to find grid
          halo_info%s_mean=modulo(halo_info%x_mean+sx,real(nf_global))

          halo_info%ang_mom=0
          do ii=1,i_odc
            dx=xv_odc(1:3,ii)-(halo_info%x_mean-([itx,ity,itz]-1)*nft)
            dv=xv_odc(4:6,ii)
            halo_info%ang_mom(1)=halo_info%ang_mom(1)+dx(2)*dv(3)-dx(3)*dv(2)
            halo_info%ang_mom(2)=halo_info%ang_mom(2)+dx(3)*dv(1)-dx(1)*dv(3)
            halo_info%ang_mom(3)=halo_info%ang_mom(3)+dx(1)*dv(2)-dx(2)*dv(1)
            halo_info%inertia(1,1)=halo_info%inertia(1,1)+dx(1)*dx(1)
            halo_info%inertia(2,2)=halo_info%inertia(2,2)+dx(2)*dx(2)
            halo_info%inertia(3,3)=halo_info%inertia(3,3)+dx(3)*dx(3)
            halo_info%inertia(1,2)=halo_info%inertia(1,2)+dx(1)*dx(2)
            halo_info%inertia(2,3)=halo_info%inertia(2,3)+dx(2)*dx(3)
            halo_info%inertia(3,1)=halo_info%inertia(3,1)+dx(3)*dx(1)
          enddo
          halo_info%inertia(2,1)=halo_info%inertia(1,2)
          halo_info%inertia(3,2)=halo_info%inertia(2,3)
          halo_info%inertia(1,3)=halo_info%inertia(3,1)
          !print*,halo_info
#ifdef PID
          if (minval(pid(ilist_odc(:i_odc)))<1) then
            print*, "  warning: pids are not all positive integers"
            !print*, pid(ilist_odc(:i_odc))
            !stop
          endif
          write(12) i_odc,pid(ilist_odc(:i_odc)),0
#         ifdef HID
            hid(pid(ilist_odc(:i_odc)))=nhalo
#         endif
#endif
          !write(101) i_odc,GroupOffset
          temp_halo(1,nhalo)=i_odc
          temp_halo(2,nhalo)=GroupOffset
          !write(102) pid(ilist_odc(:i_odc))
          GroupOffset=GroupOffset+i_odc
          header_halo_tab%Nids=header_halo_tab%Nids+i_odc
#ifdef PID
          ! compute inertia in q-space
          halo_info%q_mean=0
          nptemp=0;
          do ii=1,i_odc
            qtemp=pid(ilist_odc(ii))-1
            if (qtemp<0) then
              cycle
            endif
            nptemp=nptemp+1
            qpos(3)=qtemp/nf_global**2
            qpos(2)=(qtemp-(qtemp/nf_global**2)*nf_global**2)/nf_global
            qpos(1)=modulo(qtemp,nf_global)
            qpos=qpos+0.5
            dx=qpos-halo_info%x_mean
            dx=modulo(dx+nf_global/2,real(nf_global))-nf_global/2
            dx_mean=dx_mean+dx
          enddo
          dx_mean=dx_mean/nptemp
          qpos_mean=halo_info%x_mean+dx_mean
          qpos_mean=modulo(qpos_mean,real(nf_global))
          halo_info%q_mean=qpos_mean

          do ii=1,i_odc
            qtemp=pid(ilist_odc(ii))-1
            if (qtemp<0) then
              cycle
            endif
            qpos(3)=qtemp/nf_global**2
            qpos(2)=(qtemp-(qtemp/nf_global**2)*nf_global**2)/nf_global
            qpos(1)=modulo(qtemp,nf_global)
            qpos=qpos+0.5
            dx=qpos+0.5-qpos_mean
            dx=modulo(dx+nf_global/2,real(nf_global))-nf_global/2
            ! initial inertia
            inert(1,1)=inert(1,1)+dx(1)**2
            inert(2,2)=inert(2,2)+dx(2)**2
            inert(3,3)=inert(3,3)+dx(3)**2
            inert(1,2)=inert(1,2)+dx(1)*dx(2)
            inert(2,3)=inert(2,3)+dx(2)*dx(3)
            inert(3,1)=inert(3,1)+dx(3)*dx(1)
          enddo
          inert(2,1)=inert(1,2); inert(3,2)=inert(2,3); inert(1,3)=inert(3,1)
          halo_info%inertia_q=inert
          if (nhalo==1) then
            print*,'i_odc',i_odc
            print*,'hpos',halo_info%hpos
            print*,'x_mean',halo_info%x_mean
            print*,'v_mean',halo_info%v_mean
            print*,'ang_mom',halo_info%ang_mom
            print*,'q_mean',halo_info%q_mean
            print*,'inertia_q',halo_info%inertia_q
            !stop
          endif
#endif
          write(11) halo_info ! if the maximum is in physical regions
        else
#         ifdef analysis
            print*,'  it is a halo in buffer region'
#         endif
        endif ! physical_halo
      else
#ifdef analysis
        print*,'  it is not a main halo'
#endif
      endif

#ifdef analysis
  rewind(23)
  write(23) halo_info%hpos, halo_info%x_mean, halo_info%radius_odc, halo_info%mass_odc, i_odc
  close(23)
# ifdef PID
    close(24)
# endif
  print*,'  enter anything to continue'
  read(*,*)
#endif


    enddo ! iloc
    print*,'  found nhalo =',nhalo
    print*,'  total particles in halo =',GroupOffset,header_halo_tab%Nids
  enddo
  enddo
  enddo ! end looping over tiles
  sync all

  nhalo_tot=0
  !header_halo_tab%TotNids=0
  do ii=1,nn**3
    nhalo_tot=nhalo_tot+nhalo[ii]
    !header_halo_tab%TotNids=header_halo_tab%TotNids+header_halo_tab[ii]%Nids
  enddo
  sync all
  if (head) then
    print*, '  found',nhalo_tot,'halos'
  endif
  rewind(11)
  write(11) nhalo_tot,nhalo
  close(11)
  close(12)
#ifdef HID
  write(21) hid
  close(21)
  print*,'  count(hid>0) =', count(hid>0,kind=4)
#endif

  !header_halo_tab%Ngroups=nhalo
  !header_halo_tab%TotNgroups=nhalo_tot
  !header_halo_tab%NFiles=nn**3
  !rewind(101)
  !write(101) header_halo_tab
  !write(101) temp_halo(1,:nhalo),temp_halo(2,:nhalo)
  !close(101)

  !rewind(102)
  !write(102) header_halo_tab,0
  !close(102)

  sync all

endsubroutine

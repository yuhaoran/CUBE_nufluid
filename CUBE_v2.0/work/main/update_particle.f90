module update_particle

contains

subroutine update_x
  use variables
  use neutrinos
  implicit none
  save
  call update_xp
endsubroutine

subroutine update_xp()
  use variables
  implicit none
  save

  integer(4) ilayer,nlayer
  integer(8) idx,np_prev,checkv0,checkv1,ig(3)
  integer(4) rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! double buffer tile
  real(4) vfield_new(3,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  real(8),parameter :: weight_v=0.1 ! how previous-step vfield is mostly weighted
  integer(8) idx_ex_r(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
  integer(8),dimension(nt,nt) :: pp_l,pp_r,ppe_l,ppe_r
  integer(4) rholocal(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb) ! count writing

  if (head) print*,'update_xp'
  call system_clock(t1,t_rate)
  dt_mid=(dt_old+dt)/2
  if (head) print*,'  dt_mid =',dt_mid
  overhead_tile=0
  np_prev=0

  nlayer=2*ceiling(dt_mid*sim%vz_max/ncell)+1
  if (head) print*,'  nlayer =',nlayer

  do itz=1,nnt ! loop over tile
  do ity=1,nnt
  do itx=1,nnt
    !if (head) print*,'  tile',int(itx,1),int(ity,1),int(itz,1)
    rhoce=0
    rholocal=0
    ! for empty coarse cells, use previous-step vfield
    vfield_new=0
    vfield_new(:,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb)=vfield(:,:,:,:,itx,ity,itz)*weight_v
    xp_new=0
    vp_new=0
#   ifdef PID
    pid_new=0
#   endif
    !if (head) print*,'    density loop'
    do ilayer=0,nlayer-1
      !$omp paralleldo default(shared) schedule(static,2)&
      !$omp& private(k,j,i,np,nzero,l,ip,xq,vreal,deltax,ig)
      do k=1-ncb+ilayer,nt+ncb,nlayer
      !do k=1-ncb,nt+ncb
      do j=1-ncb,nt+ncb
      do i=1-ncb,nt+ncb
        np=rhoc(i,j,k,itx,ity,itz)
        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
        do l=1,np
          ip=nzero+l
          xq=((/i,j,k/)-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
          vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
          deltax=(dt_mid*vreal)/ncell
          ig=ceiling(xq+deltax)
          rhoce(ig(1),ig(2),ig(3))=rhoce(ig(1),ig(2),ig(3))+1
          vfield_new(:,ig(1),ig(2),ig(3))=vfield_new(:,ig(1),ig(2),ig(3))+vreal
        enddo
      enddo
      enddo
      enddo ! k
      !$omp endparalleldo
    enddo ! ilayer
    ! vfield_new is kept the same as previous-step if the grid is empty.
    ! vfield_new is at most weight_v weighted by previous-step for non-empty grids.
    vfield_new(1,:,:,:)=vfield_new(1,:,:,:)/(rhoce+weight_v)
    vfield_new(2,:,:,:)=vfield_new(2,:,:,:)/(rhoce+weight_v)
    vfield_new(3,:,:,:)=vfield_new(3,:,:,:)/(rhoce+weight_v)

    call spine_tile(rhoce,idx_ex_r,pp_l,pp_r,ppe_l,ppe_r)
    overhead_tile=max(overhead_tile,idx_ex_r(nt+2*ncb,nt+2*ncb)/real(np_tile_max))

    if (idx_ex_r(nt+2*ncb,nt+2*ncb)>np_tile_max) then
      print*, '  error: too many particles in this tile+buffer'
      print*, '  ',idx_ex_r(nt+2*ncb,nt+2*ncb),'>',np_tile_max
      print*, '  on',this_image(), itx,ity,itz
      print*, '  please set tile_buffer larger'
      stop
    endif
    ! create x_new v_new for this local tile
    !if (head) print*,'    xvnew loop'
    do ilayer=0,nlayer-1
      !$omp paralleldo default(shared) schedule(static,2)&
      !$omp& private(k,j,i,np,nzero,l,ip,xq,vreal,deltax,ig,idx)
      do k=1-ncb+ilayer,nt+ncb,nlayer
      !do k=1-ncb,nt+ncb
      do j=1-ncb,nt+ncb
      do i=1-ncb,nt+ncb
        np=rhoc(i,j,k,itx,ity,itz)
        nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
        do l=1,np
          ip=nzero+l
          xq=((/i,j,k/)-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
          vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
          deltax=(dt_mid*vreal)/ncell
          ig=ceiling(xq+deltax)
          rholocal(ig(1),ig(2),ig(3))=rholocal(ig(1),ig(2),ig(3))+1
          idx=idx_ex_r(ig(2),ig(3))-sum(rhoce(ig(1):,ig(2),ig(3)))+rholocal(ig(1),ig(2),ig(3))
          xp_new(:,idx)=xp(:,ip)+nint(dt_mid*vreal/(x_resolution*ncell))
          vreal=vreal-vfield_new(:,ig(1),ig(2),ig(3))
          vp_new(:,idx)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
#         ifdef PID
            pid_new(idx)=pid(ip)
#         endif
        enddo
      enddo
      enddo
      enddo ! k
      !$omp endparalleldo
    enddo ! ilayer
    sync all
    ! delete particles
    !if (head) print*,'    delete_particle loop'
    !$omp paralleldo default(shared) schedule(static,2)&
    !$omp& private(k,j)
    do k=1,nt
    do j=1,nt
      xp(:,np_prev+pp_l(j,k):np_prev+pp_r(j,k))=xp_new(:,ppe_l(j,k):ppe_r(j,k))
      vp(:,np_prev+pp_l(j,k):np_prev+pp_r(j,k))=vp_new(:,ppe_l(j,k):ppe_r(j,k))
#     ifdef PID
        pid(np_prev+pp_l(j,k):np_prev+pp_r(j,k))=pid_new(ppe_l(j,k):ppe_r(j,k))
#     endif
    enddo
    enddo
    !$omp endparalleldo
    np_prev=np_prev+pp_r(nt,nt)
    ! update rhoc
    rhoc(:,:,:,itx,ity,itz)=0
    rhoc(1:nt,1:nt,1:nt,itx,ity,itz)=rhoce(1:nt,1:nt,1:nt)
    vfield(:,:,:,:,itx,ity,itz)=0
    vfield(:,1:nt,1:nt,1:nt,itx,ity,itz)=vfield_new(:,1:nt,1:nt,1:nt)
  enddo
  enddo
  enddo ! end looping over tiles
  sim%nplocal=np_prev
  xp(:,sim%nplocal+1:)=0
  vp(:,sim%nplocal+1:)=0
# ifdef PID
    pid(sim%nplocal+1:)=0
    print*, '  check PID range:',minval(pid(:sim%nplocal)),maxval(pid(:sim%nplocal))
    if (minval(pid(:sim%nplocal))<1) then
      print*, '  pid are not all positive'
      print*, '  number of 0s =',count(pid(:sim%nplocal)==0)
      print*, minloc(pid(:sim%nplocal))
      !stop
    endif
# endif
  !sim%nplocal=sum(rhoc(1:nt,1:nt,1:nt,:,:,:))
  !print*, iright,sum(rhoc(1:nt,1:nt,1:nt,:,:,:))
  !stop
  sync all

  ! calculate std of the velocity field
  !cum=cumsum6(rhoc)
  call spine_image(rhoc,idx_b_l,idx_b_r,ppl0,pplr,pprl,ppr0,ppl,ppr)
  sum_c=0; sum_r=0; sum_s=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    !$omp paralleldo default(shared) schedule(static,2)&
    !$omp& private(k,j,i,np,nzero,l,ip,vreal) &
    !$omp& reduction(+:sum_c,sum_r,sum_s)
    do k=1,nt
    do j=1,nt
    do i=1,nt
      sum_c=sum_c+sum(vfield(:,i,j,k,itx,ity,itz)**2)
      np=rhoc(i,j,k,itx,ity,itz)
      nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
      do l=1,np
        ip=nzero+l
        vreal=tan(pi*real(vp(:,ip))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
        sum_r=sum_r+sum(vreal**2)
        vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
        sum_s=sum_s+sum(vreal**2)
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all
  std_vsim_c=sum_c
  std_vsim_res=sum_r
  std_vsim=sum_s
  ! co_sum
  if (head) then
    do i=2,nn**3
      std_vsim_c=std_vsim_c+std_vsim_c[i]
      std_vsim_res=std_vsim_res+std_vsim_res[i]
      std_vsim=std_vsim+std_vsim[i]
    enddo
  endif
  sync all
  ! broadcast
  std_vsim_c=std_vsim_c[1]
  std_vsim_res=std_vsim_res[1]
  std_vsim=std_vsim[1]
  sync all

  ! divide
  std_vsim=sqrt(std_vsim/sim%npglobal)
  std_vsim_c=sqrt(std_vsim_c/nc/nc/nc/nn/nn/nn)
  std_vsim_res=sqrt(std_vsim_res/sim%npglobal)

  ! set sigma_vi_new according to particle statistics
  sigma_vi_new=std_vsim_res/sqrt(3.)
  if (head) then
    print*,'  std_vsim    ',std_vsim*sim%vsim2phys,'km/s'
    print*,'  std_vsim_c  ',std_vsim_c*sim%vsim2phys,'km/s'
    print*,'  std_vsim_res',std_vsim_res*sim%vsim2phys,'km/s'
    print*,'  sigma_vi    ',sigma_vi,'(simulation unit)'
    print*,'  sigma_vi_new',sigma_vi_new,'(simulation unit)'
    write(77) a-da,real([std_vsim,std_vsim_c,std_vsim_res])*sim%vsim2phys
  endif
  sync all

  do i=1,nn**3
    overhead_tile=max(overhead_tile,overhead_tile[i])
  enddo

  if (head) then
    print*, '  tile overhead',overhead_tile*100,'% full'
    print*, '  comsumed tile_buffer =',overhead_tile*tile_buffer,'/',tile_buffer
    print*, '  clean buffer of rhoc'
  endif

  if (head) then
    npcheck=0
    do i=1,nn**3
      npcheck=npcheck+sim[i]%nplocal
    enddo
    print*, '  npcheck,npglobal=', npcheck,sim%npglobal
    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
    print*, ''
  endif
  sync all
endsubroutine update_xp

endmodule

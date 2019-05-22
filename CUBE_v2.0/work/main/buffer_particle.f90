module buffer_particle_subroutines

contains

subroutine buffer_x
  use variables
  use neutrinos
  implicit none
  save
  call buffer_xvp(.true.,.false.)
endsubroutine

subroutine buffer_v
  use variables
  use neutrinos
  implicit none
  save
  call buffer_xvp(.false.,.true.)
endsubroutine

subroutine buffer_xvp(dox,dov)
  use variables
  implicit none
  save
  logical dox,dov
  integer(8),dimension(nt,nt,nnt,nnt) :: myzl,myzr
  integer(8),dimension(1-ncb:nt+ncb,nnt,nnt) :: mzl,mzr
  integer(8),dimension(nnt,nnt) :: ml,mr

  if (head) print*, 'buffer_xvp'
  call system_clock(t1,t_rate)
  myzl=pprl(:,:,nnt,:,:)[image1d(inx,icy,icz)]
  myzr=ppr0(:,:,nnt,:,:)[image1d(inx,icy,icz)]
  sync all
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,1 ! do only tile_x=1
    !!$omp paralleldo default(shared) private(iz,iy)
    do iz=1,nt
    do iy=1,nt
      if (dox) &
      xp(:,idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
      xp(:,myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(inx,icy,icz)]
      if (dov) then
        vp(:,idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
        vp(:,myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(inx,icy,icz)]
#ifdef PID
        pid(idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
        pid(myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(inx,icy,icz)]
#endif
      endif
    enddo
    enddo
    !!$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x- buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=2,nnt ! skip tile_x=1
    !$omp paralleldo default(shared) schedule(static,2) private(iz,iy)
    do iz=1,nt
    do iy=1,nt
      if (dox) &
      xp(:,idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
      xp(:,pprl(iy,iz,itx-1,ity,itz)+1:ppr0(iy,iz,itx-1,ity,itz))
      if (dov) then
        vp(:,idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
        vp(:,pprl(iy,iz,itx-1,ity,itz)+1:ppr0(iy,iz,itx-1,ity,itz))
#ifdef PID
        pid(idx_b_l(iy,iz,itx,ity,itz)+1:ppl0(iy,iz,itx,ity,itz))=&
        pid(pprl(iy,iz,itx-1,ity,itz)+1:ppr0(iy,iz,itx-1,ity,itz))
#endif
      endif
    enddo
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! sync x+ buffer with node on the right
  myzl=ppl0(:,:,1,:,:)[image1d(ipx,icy,icz)]
  myzr=pplr(:,:,1,:,:)[image1d(ipx,icy,icz)]
  sync all
  do itz=1,nnt
  do ity=1,nnt
  do itx=nnt,nnt ! do only tile_x=nnt
    !!$omp paralleldo default(shared) private(iz,iy)
    do iz=1,nt
    do iy=1,nt
      if (dox) &
      xp(:,ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
      xp(:,myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(ipx,icy,icz)]
      if (dov) then
        vp(:,ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
        vp(:,myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(ipx,icy,icz)]
#ifdef PID
        pid(ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
        pid(myzl(iy,iz,ity,itz)+1:myzr(iy,iz,ity,itz))[image1d(ipx,icy,icz)]
#endif
      endif
    enddo
    enddo
    !!$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! redistribute local x+ buffer
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt-1 ! skip tile_x=nnt
    !$omp paralleldo default(shared) schedule(static,2) private(iz,iy)
    do iz=1,nt
    do iy=1,nt
      if (dox) &
      xp(:,ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
      xp(:,ppl0(iy,iz,itx+1,ity,itz)+1:pplr(iy,iz,itx+1,ity,itz))
      if (dov) then
        vp(:,ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
        vp(:,ppl0(iy,iz,itx+1,ity,itz)+1:pplr(iy,iz,itx+1,ity,itz))
#ifdef PID
        pid(ppr0(iy,iz,itx,ity,itz)+1:idx_b_r(iy,iz,itx,ity,itz))=&
        pid(ppl0(iy,iz,itx+1,ity,itz)+1:pplr(iy,iz,itx+1,ity,itz))
#endif
      endif
    enddo
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! buffer y direction
  ! sync y-
  mzl=idx_b_l(nt-ncb+1,:,:,nnt,:)[image1d(icx,iny,icz)]
  mzr=idx_b_r(nt,:,:,nnt,:)[image1d(icx,iny,icz)]
  sync all
  do itz=1,nnt
  do ity=1,1 ! do only ity=1
  do itx=1,nnt
    !!$omp paralleldo default(shared) private(iz)
    do iz=1,nt
      if (dox) &
      xp(:,idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
      xp(:,mzl(iz,itx,itz)+1:mzr(iz,itx,itz))[image1d(icx,iny,icz)]
      if (dov) then
        vp(:,idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
        vp(:,mzl(iz,itx,itz)+1:mzr(iz,itx,itz))[image1d(icx,iny,icz)]
#ifdef PID
        pid(idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
        pid(mzl(iz,itx,itz)+1:mzr(iz,itx,itz))[image1d(icx,iny,icz)]
#endif
      endif
    enddo
    !!$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y-
  do itz=1,nnt
  do ity=2,nnt ! skip ity=1
  do itx=1,nnt
    !$omp paralleldo default(shared) schedule(static,2) private(iz)
    do iz=1,nt
      if (dox) &
      xp(:,idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
      xp(:,idx_b_l(nt-ncb+1,iz,itx,ity-1,itz)+1:idx_b_r(nt,iz,itx,ity-1,itz))
      if (dov) then
        vp(:,idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
        vp(:,idx_b_l(nt-ncb+1,iz,itx,ity-1,itz)+1:idx_b_r(nt,iz,itx,ity-1,itz))
#ifdef PID
        pid(idx_b_l(1-ncb,iz,itx,ity,itz)+1:idx_b_r(0,iz,itx,ity,itz))=&
        pid(idx_b_l(nt-ncb+1,iz,itx,ity-1,itz)+1:idx_b_r(nt,iz,itx,ity-1,itz))
#endif
      endif
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! sync y+
  mzl=idx_b_l(1,:,:,1,:)[image1d(icx,ipy,icz)]
  mzr=idx_b_r(ncb,:,:,1,:)[image1d(icx,ipy,icz)]
  sync all
  do itz=1,nnt
  do ity=nnt,nnt ! do only ity=nnt
  do itx=1,nnt
    !!$omp paralleldo default(shared) private(iz)
    do iz=1,nt
      if (dox) &
      xp(:,idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
      xp(:,mzl(iz,itx,itz)+1:mzr(iz,itx,itz))[image1d(icx,ipy,icz)]
      if (dov) then
        vp(:,idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
        vp(:,mzl(iz,itx,itz)+1:mzr(iz,itx,itz))[image1d(icx,ipy,icz)]
#ifdef PID
        pid(idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
        pid(mzl(iz,itx,itz)+1:mzr(iz,itx,itz))[image1d(icx,ipy,icz)]
#endif
      endif
    enddo
    !!$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! redistribute y+
  do itz=1,nnt
  do ity=1,nnt-1 ! skip ity=nnt
  do itx=1,nnt
    !$omp paralleldo default(shared) schedule(static,2) private(iz)
    do iz=1,nt
      if (dox) &
      xp(:,idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
      xp(:,idx_b_l(1,iz,itx,ity+1,itz)+1:idx_b_r(ncb,iz,itx,ity+1,itz))
      if (dov) then
        vp(:,idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
        vp(:,idx_b_l(1,iz,itx,ity+1,itz)+1:idx_b_r(ncb,iz,itx,ity+1,itz))
#ifdef PID
        pid(idx_b_l(nt+1,iz,itx,ity,itz)+1:idx_b_r(nt+ncb,iz,itx,ity,itz))=&
        pid(idx_b_l(1,iz,itx,ity+1,itz)+1:idx_b_r(ncb,iz,itx,ity+1,itz))
#endif
      endif
    enddo
    !$omp endparalleldo
  enddo
  enddo
  enddo
  sync all

  ! buffer z direction
  ! sync z-
  ml=idx_b_l(1-ncb,nt-ncb+1,:,:,nnt)[image1d(icx,icy,inz)]
  mr=idx_b_r(nt+ncb,nt,:,:,nnt)[image1d(icx,icy,inz)]
  sync all
  do itz=1,1 ! do only itz=1
  !!$omp paralleldo default(shared) private(ity,itx)
  do ity=1,nnt
  do itx=1,nnt
    if (dox) &
    xp(:,idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
    xp(:,ml(itx,ity)+1:mr(itx,ity))[image1d(icx,icy,inz)]
    if (dov) then
      vp(:,idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
      vp(:,ml(itx,ity)+1:mr(itx,ity))[image1d(icx,icy,inz)]
#ifdef PID
      pid(idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
      pid(ml(itx,ity)+1:mr(itx,ity))[image1d(icx,icy,inz)]
#endif
    endif
  enddo
  enddo
  !!$omp endparalleldo
  enddo
  sync all

  ! redistribute z-
  do itz=2,nnt ! skip itz=1
  !$omp paralleldo default(shared) schedule(static,2) private(ity,itx)
  do ity=1,nnt
  do itx=1,nnt
    if (dox) &
    xp(:,idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
    xp(:,idx_b_l(1-ncb,nt-ncb+1,itx,ity,itz-1)+1:idx_b_r(nt+ncb,nt,itx,ity,itz-1))
    if (dov) then
      vp(:,idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
      vp(:,idx_b_l(1-ncb,nt-ncb+1,itx,ity,itz-1)+1:idx_b_r(nt+ncb,nt,itx,ity,itz-1))
#ifdef PID
      pid(idx_b_l(1-ncb,1-ncb,itx,ity,itz)+1:idx_b_r(nt+ncb,0,itx,ity,itz))=&
      pid(idx_b_l(1-ncb,nt-ncb+1,itx,ity,itz-1)+1:idx_b_r(nt+ncb,nt,itx,ity,itz-1))
#endif
    endif
  enddo
  enddo
  !$omp endparalleldo
  enddo
  sync all

  ! sync z+
  ml=idx_b_l(1-ncb,1,:,:,1)[image1d(icx,icy,ipz)]
  mr=idx_b_r(nt+ncb,ncb,:,:,1)[image1d(icx,icy,ipz)]
  sync all
  do itz=nnt,nnt ! do only itz=nnt
  !!$omp paralleldo default(shared) private(ity,itx)
  do ity=1,nnt
  do itx=1,nnt
    if (dox) &
    xp(:,idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
    xp(:,ml(itx,ity)+1:mr(itx,ity))[image1d(icx,icy,ipz)]
    if (dov) then
      vp(:,idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
      vp(:,ml(itx,ity)+1:mr(itx,ity))[image1d(icx,icy,ipz)]
#ifdef PID
      pid(idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
      pid(ml(itx,ity)+1:mr(itx,ity))[image1d(icx,icy,ipz)]
#endif
    endif
  enddo
  enddo
  !!$omp endparalleldo
  enddo
  sync all

  ! redistribute z+
  do itz=1,nnt-1 ! skip itz=nnt
  !$omp paralleldo default(shared) schedule(static,2) private(ity,itx)
  do ity=1,nnt
  do itx=1,nnt
    if (dox) &
    xp(:,idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
    xp(:,idx_b_l(1-ncb,1,itx,ity,itz+1)+1:idx_b_r(nt+ncb,ncb,itx,ity,itz+1))
    if (dov) then
      vp(:,idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
      vp(:,idx_b_l(1-ncb,1,itx,ity,itz+1)+1:idx_b_r(nt+ncb,ncb,itx,ity,itz+1))
#ifdef PID
      pid(idx_b_l(1-ncb,nt+1,itx,ity,itz)+1:idx_b_r(nt+ncb,nt+ncb,itx,ity,itz))=&
      pid(idx_b_l(1-ncb,1,itx,ity,itz+1)+1:idx_b_r(nt+ncb,ncb,itx,ity,itz+1))
#endif
    endif
  enddo
  enddo
  !$omp endparalleldo
  enddo
  sync all
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
# ifdef FORCETEST
    print*,'xp =',xp(:,1:20)
# endif
endsubroutine buffer_xvp


endmodule

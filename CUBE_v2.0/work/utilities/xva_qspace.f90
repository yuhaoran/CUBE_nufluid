program xva_qspace
  use parameters
  implicit none
  save
  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim
  integer(8),parameter :: npnode=nf**3 ! only true for this project
  real,parameter :: density_buffer=1.2
  integer(8),parameter :: npmax=npnode*density_buffer
  integer(8) i,j,k,l,i_dim,iq(3),pid8,nplocal,itx,ity,itz,nlast,ip,np,ipp
  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt),rhoc0(nt,nt,nt,nnt,nnt,nnt)
  real(4) rho0(0:ng+1,0:ng+1,0:ng+1),rho_grid(0:ng+1,0:ng+1,0:ng+1)
  real dsp(3,0:ng+1,0:ng+1,0:ng+1),mass_p,pos0(3),pos1(3),dpos(3)
  real vel(3,0:ng+1,0:ng+1,0:ng+1)
  real acc(3,0:ng+1,0:ng+1,0:ng+1)
  real vfield(3,nt,nt,nt,nnt,nnt,nnt)

  integer idx,idx1(3),idx2(3)
  real dx1(3),dx2(3)
  integer(izipx) xp(3,npmax)
  integer(izipv) vp(3,npmax)
  integer(8)   pid(npmax)
  real(4) afield(3,npmax),vreal(3),sigma_vi

  call geometry

  if (head) then
    print*, 'Displacement field analysis on resolution:'
    print*, 'ng=',ng
  endif
  sync all

  if (head) then
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

  do cur_checkpoint= 1,n_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    open(12,file=output_name('zip2'),status='old',action='read',access='stream')
    read(12) sim
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, 'zip format incompatable'
      close(12)
      stop
    endif
    read(12) rhoc ! coarse grid density
    close(12)
    mass_p=sim%mass_p
    print*, 'mass_p =',mass_p
    nplocal=sim%nplocal
    print*, 'nplocal =',nplocal
    open(10,file=output_name('zip0'),status='old',action='read',access='stream')
    read(10) xp(:,:nplocal)
    close(10)
    open(10,file=output_name('zip1'),status='old',action='read',access='stream')
    read(10) vp(:,:nplocal)
    close(10)
    open(10,file=output_name('zipid'),status='old',action='read',access='stream')
    read(10) pid(:nplocal)
    close(10)
    open(10,file=output_name('vfield'),status='old',action='read',access='stream')
    read(10) vfield
    close(10)
    sigma_vi=sim%sigma_vi
    open(10,file=output_name('afield'),status='old',action='read',access='stream')
    read(10) afield(:,:nplocal)
    close(10)

    rho0=0
    dsp=0
    ipp=0
    nlast=0
    acc=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      if (head) print*, 'CIC interpolation on tile',int(itx,1),int(ity,1),int(itz,1)
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l
          ipp=ipp+1
          pid8=pid(ip)-1
          iq(3)=pid8/int(nf_global,4)**2
          iq(2)=(pid8-iq(3)*int(nf_global,4)**2)/int(nf_global,4)
          iq(1)=modulo(pid8,int(nf_global,4))

          pos0=iq+0.5
          pos1=nt*((/itx,ity,itz/)-1) + (/i,j,k/)-1 + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          pos1=real(ng)*((/icx,icy,icz/)-1) + pos1*real(ng)/real(nc)

          dpos=pos1-pos0
          dpos=modulo(dpos+ng*nn/2,real(ng*nn))-ng*nn/2
          dpos=dpos*real(nf)/real(ng)
          pos0=pos0-0.5 +0.5 ! dsp- "+0.5" assigns particle to the half-right grid
                         ! then the divergence will be done between current- and right-grid.
          idx1=floor(pos0)+1
          idx2=idx1+1
          dx1=idx1-pos0
          dx2=1-dx1
          ! displacement field
          dsp(:,idx1(1),idx1(2),idx1(3))=dsp(:,idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*dpos
          dsp(:,idx2(1),idx1(2),idx1(3))=dsp(:,idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*dpos
          dsp(:,idx1(1),idx2(2),idx1(3))=dsp(:,idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*dpos
          dsp(:,idx1(1),idx1(2),idx2(3))=dsp(:,idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*dpos
          dsp(:,idx1(1),idx2(2),idx2(3))=dsp(:,idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*dpos
          dsp(:,idx2(1),idx1(2),idx2(3))=dsp(:,idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*dpos
          dsp(:,idx2(1),idx2(2),idx1(3))=dsp(:,idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*dpos
          dsp(:,idx2(1),idx2(2),idx2(3))=dsp(:,idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*dpos
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
          ! acceleration field
          acc(:,idx1(1),idx1(2),idx1(3))=acc(:,idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*afield(:,ipp)
          acc(:,idx2(1),idx1(2),idx1(3))=acc(:,idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*afield(:,ipp)
          acc(:,idx1(1),idx2(2),idx1(3))=acc(:,idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*afield(:,ipp)
          acc(:,idx1(1),idx1(2),idx2(3))=acc(:,idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*afield(:,ipp)
          acc(:,idx1(1),idx2(2),idx2(3))=acc(:,idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*afield(:,ipp)
          acc(:,idx2(1),idx1(2),idx2(3))=acc(:,idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*afield(:,ipp)
          acc(:,idx2(1),idx2(2),idx1(3))=acc(:,idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*afield(:,ipp)
          acc(:,idx2(1),idx2(2),idx2(3))=acc(:,idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*afield(:,ipp)
          ! CIC number count
          rho0(idx1(1),idx1(2),idx1(3))=rho0(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)
          rho0(idx2(1),idx1(2),idx1(3))=rho0(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)
          rho0(idx1(1),idx2(2),idx1(3))=rho0(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)
          rho0(idx1(1),idx1(2),idx2(3))=rho0(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)
          rho0(idx1(1),idx2(2),idx2(3))=rho0(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)
          rho0(idx2(1),idx1(2),idx2(3))=rho0(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)
          rho0(idx2(1),idx2(2),idx1(3))=rho0(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)
          rho0(idx2(1),idx2(2),idx2(3))=rho0(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)
          ! CIC final density
          pos1=pos1-0.5
          idx1=floor(pos1)+1
          idx2=idx1+1
          dx1=idx1-pos1
          dx2=1-dx1
          rho_grid(idx1(1),idx1(2),idx1(3))=rho_grid(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p
          rho_grid(idx2(1),idx1(2),idx1(3))=rho_grid(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p
          rho_grid(idx1(1),idx2(2),idx1(3))=rho_grid(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p
          rho_grid(idx1(1),idx1(2),idx2(3))=rho_grid(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p
          rho_grid(idx1(1),idx2(2),idx2(3))=rho_grid(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p
          rho_grid(idx2(1),idx1(2),idx2(3))=rho_grid(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p
          rho_grid(idx2(1),idx2(2),idx1(3))=rho_grid(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p
          rho_grid(idx2(1),idx2(2),idx2(3))=rho_grid(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    print*, 'sum of rho_grid',sum(rho_grid*1d0)
    print*, 'Start sync from buffer regions'
    dsp(:,1,:,:)=dsp(:,1,:,:)+dsp(:,ng+1,:,:)
    dsp(:,ng,:,:)=dsp(:,ng,:,:)+dsp(:,0,:,:)
    dsp(:,:,1,:)=dsp(:,:,1,:)+dsp(:,:,ng+1,:)
    dsp(:,:,ng,:)=dsp(:,:,ng,:)+dsp(:,:,0,:)
    dsp(:,:,:,1)=dsp(:,:,:,1)+dsp(:,:,:,ng+1)
    dsp(:,:,:,ng)=dsp(:,:,:,ng)+dsp(:,:,:,0)

    vel(:,1,:,:)=vel(:,1,:,:)+vel(:,ng+1,:,:)
    vel(:,ng,:,:)=vel(:,ng,:,:)+vel(:,0,:,:)
    vel(:,:,1,:)=vel(:,:,1,:)+vel(:,:,ng+1,:)
    vel(:,:,ng,:)=vel(:,:,ng,:)+vel(:,:,0,:)
    vel(:,:,:,1)=vel(:,:,:,1)+vel(:,:,:,ng+1)
    vel(:,:,:,ng)=vel(:,:,:,ng)+vel(:,:,:,0)

    acc(:,1,:,:)=acc(:,1,:,:)+acc(:,ng+1,:,:)
    acc(:,ng,:,:)=acc(:,ng,:,:)+acc(:,0,:,:)
    acc(:,:,1,:)=acc(:,:,1,:)+acc(:,:,ng+1,:)
    acc(:,:,ng,:)=acc(:,:,ng,:)+acc(:,:,0,:)
    acc(:,:,:,1)=acc(:,:,:,1)+acc(:,:,:,ng+1)
    acc(:,:,:,ng)=acc(:,:,:,ng)+acc(:,:,:,0)

    rho0(1,:,:)=rho0(1,:,:)+rho0(ng+1,:,:)
    rho0(ng,:,:)=rho0(ng,:,:)+rho0(0,:,:)
    rho0(:,1,:)=rho0(:,1,:)+rho0(:,ng+1,:)
    rho0(:,ng,:)=rho0(:,ng,:)+rho0(:,0,:)
    rho0(:,:,1)=rho0(:,:,1)+rho0(:,:,ng+1)
    rho0(:,:,ng)=rho0(:,:,ng)+rho0(:,:,0)

    rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(ng+1,:,:)
    rho_grid(ng,:,:)=rho_grid(ng,:,:)+rho_grid(0,:,:)
    rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,ng+1,:)
    rho_grid(:,ng,:)=rho_grid(:,ng,:)+rho_grid(:,0,:)
    rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,ng+1)
    rho_grid(:,:,ng)=rho_grid(:,:,ng)+rho_grid(:,:,0)
    print*, 'check: min,max,sum of rho0 = '
    print*, minval(rho0(1:ng,1:ng,1:ng)),maxval(rho0(1:ng,1:ng,1:ng)),sum(rho0(1:ng,1:ng,1:ng))
    print*, 'check: min,max,sum of rho_grid = '
    print*, minval(rho_grid(1:ng,1:ng,1:ng)),maxval(rho_grid(1:ng,1:ng,1:ng)),sum(rho_grid(1:ng,1:ng,1:ng))

    do i_dim=1,3
      print*, 'dsp: dimension',int(i_dim,1),'min,max values ='
      print*, minval(dsp(i_dim,:,:,:)), maxval(dsp(i_dim,:,:,:))
    enddo

    if (head) print*,'Write into file'
    open(15,file=output_name('den'),status='replace',access='stream')
      write(15) rho_grid(1:ng,1:ng,1:ng)
    close(15)
    open(15,file=output_name('dsp'),status='replace',access='stream')
      write(15) dsp(1,1:ng,1:ng,1:ng)
      write(15) dsp(2,1:ng,1:ng,1:ng)
      write(15) dsp(3,1:ng,1:ng,1:ng)
    close(15)
    open(15,file=output_name('vel'),status='replace',access='stream')
      write(15) vel(1,1:ng,1:ng,1:ng)
      write(15) vel(2,1:ng,1:ng,1:ng)
      write(15) vel(3,1:ng,1:ng,1:ng)
    close(15)
    open(15,file=output_name('acc'),status='replace',access='stream')
      write(15) acc(1,1:ng,1:ng,1:ng)
      write(15) acc(2,1:ng,1:ng,1:ng)
      write(15) acc(3,1:ng,1:ng,1:ng)
    close(15)



  enddo
  print*,'xva_qspace done'
end

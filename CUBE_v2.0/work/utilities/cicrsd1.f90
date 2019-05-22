!! Include RSD effect of particles, and CIC interpolate into grid to get the redshift-space density field
!#define RSD_ELUCID
program cicrsd
  use parameters
  implicit none
  save
  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim
  real,parameter :: density_buffer=1.2
  integer,parameter :: ngrid=500
  real,parameter :: cen(3)=[30,370,370]
  integer(8) i,j,k,l,i_dim,iq(3),nplocal,nplocal_nu,itx,ity,itz
  integer(8) nlast,ip,np,idx1(3),idx2(3),pid8

  real(4) rho_grid(0:ngrid+1,0:ngrid+1,0:ngrid+1), dsp(3,0:ngrid+1,0:ngrid+1,0:ngrid+1)
  real(4) rho_c(ngrid,ngrid,ngrid),sigma_vi,zshift,los(3),sx(3)
  real(4) mass_p,pos0(3),dpos(3),pos1(3),dx1(3),dx2(3)
  real(8) rho8
  integer(izipx),allocatable :: xp(:,:)
  integer(izipv),allocatable :: vp(:,:)
  integer(4),allocatable :: pid(:)
  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt)
  real(4) vc(3,nt,nt,nt,nnt,nnt,nnt)

  character(20) str_z,str_i
  call geometry
    print*, 'cicrsd-ELUCID on resolution:'
    print*, 'ngrid=',ngrid
    print*, 'ng*nn=',ng*nn
    print*, 'cen=',cen
    print*, 'checkpoint at:'
    open(16,file='../main/z_checkpoint.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
      print*, z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
    print*,''

  do cur_checkpoint= n_checkpoint,n_checkpoint
    print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))

    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim
    close(11)
    sigma_vi=sim%sigma_vi
    !call print_header(sim); stop
    !mass_p=sim%mass_p
    mass_p=1.
    nplocal=sim%nplocal
    print*, 'mass_p =',mass_p
    print*, 'nplocal =',nplocal
    !cdm
    allocate(xp(3,nplocal),vp(3,nplocal),pid(sim%nplocal))
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp
    close(11)
    open(11,file=output_name('vp'),status='old',action='read',access='stream')
    read(11) vp
    close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc
    close(11)
    open(11,file=output_name('vc'),status='old',action='read',access='stream')
    read(11) vc
    close(11)
    open(14,file=output_name('id'),status='old',action='read',access='stream')
    read(14) pid
    close(14)
    print*,'check PID range:',minval(pid),maxval(pid)

    rho_grid=0
    dsp=0
    nlast=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l

          ! pos0
          pid8=pid(ip)-1
          if (pid8>0) then
            iq(3)=pid8/int(nf_global,4)**2
            iq(2)=(pid8-iq(3)*int(nf_global,4)**2)/int(nf_global,4)
            iq(1)=modulo(pid8,int(nf_global,4))
            pos0=(iq+0.5)*ngrid/nf_global
          endif
          ! pos1
          pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          pos1=pos1*real(ngrid)/real(nc)
#ifdef RSD_ELUCID
          sx=vc(:,i,j,k,itx,ity,itz)
          !print*,sx
          sx=sx+tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
          !print*,sx
          los=pos1-(cen/500.)*ngrid ! line of sight vector
          !print*,los
          los=los/norm2(los) ! line of sight unit vector
          !print*,los
          sx=sum(sx*los)*los ! project velocity to line of sight
          !print*,'sx',sx
          sx=sx*sim%vsim2phys/sim%a/(100*h0) ! convert to km/h and multiply 1/aH, in Mpc
          !print*,'sx',sx
          sx=sx/(h0*box/ngrid) ! convert to find grid
          !print*,'sx',sx
          !print*,'pos1',pos1
          pos1=pos1+sx
          !print*,'pos1',pos1
          !stop
          pos1=modulo(pos1,real(ngrid))
#endif
          dpos=pos0-pos1
          ! if pid==0, use pos0 of the previous particle
          dpos=modulo(dpos+ngrid/2,real(ngrid))-ngrid/2
          ! interpolation
          pos1=pos1 - 0.5
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
          dsp(:,idx1(1),idx1(2),idx1(3))=dsp(:,idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*dpos
          dsp(:,idx2(1),idx1(2),idx1(3))=dsp(:,idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*dpos
          dsp(:,idx1(1),idx2(2),idx1(3))=dsp(:,idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*dpos
          dsp(:,idx1(1),idx1(2),idx2(3))=dsp(:,idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*dpos
          dsp(:,idx1(1),idx2(2),idx2(3))=dsp(:,idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*dpos
          dsp(:,idx2(1),idx1(2),idx2(3))=dsp(:,idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*dpos
          dsp(:,idx2(1),idx2(2),idx1(3))=dsp(:,idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*dpos
          dsp(:,idx2(1),idx2(2),idx2(3))=dsp(:,idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*dpos
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    sync all
    deallocate(xp,vp,pid)
    print*, 'Start sync from buffer regions'
    sync all
    rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(ngrid+1,:,:); sync all
    rho_grid(ngrid,:,:)=rho_grid(ngrid,:,:)+rho_grid(0,:,:); sync all
    rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,ngrid+1,:); sync all
    rho_grid(:,ngrid,:)=rho_grid(:,ngrid,:)+rho_grid(:,0,:); sync all
    rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,ngrid+1); sync all
    rho_grid(:,:,ngrid)=rho_grid(:,:,ngrid)+rho_grid(:,:,0); sync all

    dsp(:,1,:,:)=dsp(:,1,:,:)+dsp(:,ngrid+1,:,:); sync all
    dsp(:,ngrid,:,:)=dsp(:,ngrid,:,:)+dsp(:,0,:,:); sync all
    dsp(:,:,1,:)=dsp(:,:,1,:)+dsp(:,:,ngrid+1,:); sync all
    dsp(:,:,ngrid,:)=dsp(:,:,ngrid,:)+dsp(:,:,0,:); sync all
    dsp(:,:,:,1)=dsp(:,:,:,1)+dsp(:,:,:,ngrid+1); sync all
    dsp(:,:,:,ngrid)=dsp(:,:,:,ngrid)+dsp(:,:,:,0); sync all

print*,dsp(1,1:10,1,1)
print*,''
print*,rho_grid(1:10,1,1)

sync all
    !rho_c=rho_grid(1:ng,1:ng,1:ng)
    do k=1,ngrid
    do j=1,ngrid
    do i=1,ngrid
      rho_c(i,j,k)=rho_grid(i,j,k)
      if (rho_grid(i,j,k)/=0)      dsp(1,i,j,k)=dsp(1,i,j,k)/rho_grid(i,j,k)
      if (rho_grid(i,j,k)/=0)      dsp(2,i,j,k)=dsp(2,i,j,k)/rho_grid(i,j,k)
      if (rho_grid(i,j,k)/=0)      dsp(3,i,j,k)=dsp(3,i,j,k)/rho_grid(i,j,k)
    enddo
    enddo
    enddo
sync all
!rho_grid=rho_grid+0.0001
!dsp(1,:,:,:)=dsp(1,:,:,:)/rho_grid
!dsp(2,:,:,:)=dsp(2,:,:,:)/rho_grid
!dsp(3,:,:,:)=dsp(3,:,:,:)/rho_grid
!rho_grid=rho_grid-0.0001

print*,''
print*,dsp(1,1:10,1,1)
!stop

    print*,rho_grid(1:2,1:2,1)
    print*,rho_grid(1:2,1:2,2)
    print*,rho_c(1:2,1:2,1)
    print*,rho_c(1:2,1:2,2)

    print*, 'check: min,max,sum of rho_grid = '
    print*, minval(rho_c),maxval(rho_c),sum(rho_c*1d0)

    rho8=sum(rho_c*1d0); sync all
    ! co_sum
    print*,'rho_global',rho8,ng_global
    ! convert to density contrast
    do i=1,ngrid
      rho_c(:,:,i)=rho_c(:,:,i)/(rho8/ng_global/ng_global/ng_global)-1
    enddo
    !rho_c=rho_c/(rho8/ng_global/ng_global/ng_global)-1
    ! check normalization
    print*,'min',minval(rho_c),'max',maxval(rho_c),'mean',sum(rho_c*1d0)/ngrid/ngrid/ngrid; sync all

#ifdef RSD_ELUCID
    open(11,file=output_name('delta_rsd'),status='replace',access='stream')
#else
    open(11,file=output_name('delta_c'),status='replace',access='stream')
#endif
    write(11) rho_c
    close(11)

#ifdef RSD_ELUCID
    open(11,file=output_name('idsp_rsd'),status='replace',access='stream')
#else
    open(11,file=output_name('idsp_c'),status='replace',access='stream')
#endif
    write(11) dsp(:,1:ngrid,1:ngrid,1:ngrid)
    close(11)

  enddo
  print*,'cicrsd done'
end

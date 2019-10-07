!! add -DNEUTRINOS to compute cross_power(delta_c,delta_nu)
!! otherwise compute cross_power(delta_L,delta_c)

!#define write_xreal
!#define analysis
!define RSD
!#define RSD_Emode
program cicpower
  use parameters
  use pencil_fft
  use powerspectrum
  use hydrodf
  implicit none
  save

  type(dfld) :: neu

  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim
  real,parameter :: density_buffer=1.2
  integer(8) i,j,k,l,i_dim,iq(3),nplocal,nplocal_nu,itx,ity,itz
  integer(8) nlast,ip,np,idx1(3),idx2(3)

  real(4) rho_grid(0:ng+1,0:ng+1,0:ng+1)[*]
  real(4) rho_c(ng,ng,ng),rho_nu(ng,ng,ng)
  real(4) mass_p,pos1(3),dx1(3),dx2(3)
  real(8) rho8[*]

  integer(izipx),allocatable :: xp(:,:)
  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt)

  real xi(10,nbin)[*]
  character(20) str_z,str_i
  call geometry
  if (head) then
    print*, 'cicpower on resolution:'
    print*, 'ng=',ng
    print*, 'ng*nn=',ng*nn
  endif
  sync all

  if (head) then
    print*, 'checkpoint at:'
    open(16,file='../main/z_checkpoint.txt',status='old')
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

  do cur_checkpoint= 2,n_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))

    !call particle_initialization

    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim
    close(11)
    !call print_header(sim); stop
    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
      print*, 'zip format incompatable'
      close(11)
      stop
    endif
    !mass_p=sim%mass_p
    mass_p=1.
    nplocal=sim%nplocal
    nplocal_nu=sim%nplocal_nu
    if (head) then
      print*, 'mass_p =',mass_p
      print*, 'nplocal =',nplocal
      print*, 'nplocal_nu =',nplocal_nu
    endif
    !cdm
    allocate(xp(3,nplocal))
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp
    close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc
    close(11)

#ifdef write_xreal
    open(12,file=output_name('xreal'),status='replace',access='stream')
#endif

    rho_grid=0
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
          pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution

#ifdef write_xreal
    write(12) pos1*real(ng)/real(nc)
#endif

          pos1=pos1*real(ng)/real(nc) - 0.5

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
    sync all
    deallocate(xp)

    if (head) print*, 'Start sync from buffer regions'
    sync all
    rho_grid(1,:,:)=rho_grid(1,:,:)+rho_grid(ng+1,:,:)[image1d(inx,icy,icz)]
    rho_grid(ng,:,:)=rho_grid(ng,:,:)+rho_grid(0,:,:)[image1d(ipx,icy,icz)]; sync all
    rho_grid(:,1,:)=rho_grid(:,1,:)+rho_grid(:,ng+1,:)[image1d(icx,iny,icz)]
    rho_grid(:,ng,:)=rho_grid(:,ng,:)+rho_grid(:,0,:)[image1d(icx,ipy,icz)]; sync all
    rho_grid(:,:,1)=rho_grid(:,:,1)+rho_grid(:,:,ng+1)[image1d(icx,icy,inz)]
    rho_grid(:,:,ng)=rho_grid(:,:,ng)+rho_grid(:,:,0)[image1d(icx,icy,ipz)]; sync all
    !rho_c=rho_grid(1:ng,1:ng,1:ng)
    do i=1,ng
    do j=1,ng
    do k=1,ng
      rho_c(k,j,i)=rho_grid(k,j,i)
    enddo
    enddo
    enddo
    print*,rho_grid(1:2,1:2,1)
    print*,rho_grid(1:2,1:2,2)
    print*,rho_c(1:2,1:2,1)
    print*,rho_c(1:2,1:2,2)
    print*, 'check: min,max,sum of rho_grid = '
    print*, minval(rho_c),maxval(rho_c),sum(rho_c*1d0)

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
    do i=1,ng
      rho_c(:,:,i)=rho_c(:,:,i)/(rho8/ng_global/ng_global/ng_global)-1
    enddo
    !rho_c=rho_c/(rho8/ng_global/ng_global/ng_global)-1
    ! check normalization
    print*,'min',minval(rho_c),'max',maxval(rho_c),'mean',sum(rho_c*1d0)/ng/ng/ng; sync all
    print*,'sample'
    print*,rho_c(1:2,1:2,1)
    print*,rho_c(1:2,1:2,2)



    if (head) print*,'Write delta_c into',output_name('delta_c')
    open(11,file=output_name('delta_c'),status='replace',access='stream')
    write(11) rho_c
    close(11); sync all
    open(11,file=output_name('delta_c_proj'),status='replace',access='stream')
    write(11) sum(rho_c(:,:,:),dim=3)/ng
    close(11); sync all

    ! power spectrum
    write(str_i,'(i6)') image
    write(str_z,'(f7.3)') z_checkpoint(cur_checkpoint)


  
  open(15,file=output_dir()//'delta_L'//output_suffix(),status='old',access='stream')
  read(15) rho_nu
  close(15)

  !Read in neutrinos
  if (sum(f_neu).gt.0) then
     if (head) write(*,*) 'computing neu density'
     call neu%restart(output_name('neu1'))
     !call neu%restart(output_dir()//'neu1'//output_suffix())
     do k=1,ng
     do j=1,ng
     do i=1,ng
        dx1(1)=(i-0.5)*hg_nf/ng
        dx1(2)=(j-0.5)*hg_nf/ng
        dx1(3)=(k-0.5)*hg_nf/ng
        rho_nu(i,j,k)=neu%density(dx1,int(nf/ng))
     end do
     end do
     end do
     rho_nu=rho_nu-1.d0
  else
     rho_nu=0
  end if

  open(17,file='./'//trim(adjustl(str_z))//'slice_x.bin',access='stream',status='replace')
  write(17) rho_c(ng/2,:,:)
  close(17)

  open(17,file='./'//trim(adjustl(str_z))//'slice_y.bin',access='stream',status='replace')
  write(17) rho_c(:,ng/2,:)
  close(17)

  open(17,file='./'//trim(adjustl(str_z))//'slice_z.bin',access='stream',status='replace')
  write(17) rho_c(:,:,ng/2)
  close(17)

  open(17,file='./'//trim(adjustl(str_z))//'nuslice_x.bin',access='stream',status='replace')
  write(17) rho_nu(ng/2,:,:)
  close(17)

  open(17,file='./'//trim(adjustl(str_z))//'nuslice_y.bin',access='stream',status='replace')
  write(17) rho_nu(:,ng/2,:)
  close(17)

  open(17,file='./'//trim(adjustl(str_z))//'nuslice_z.bin',access='stream',status='replace')
  write(17) rho_nu(:,:,ng/2)
  close(17)


    call cross_power(xi,rho_c,rho_nu)
  sync all
  if (head) then
    open(15,file=output_name('cicpower'),status='replace',access='stream')
    write(15) xi
    close(15)

    open(15,file='./'//trim(adjustl(str_z))//'cicpower.txt',status='replace')
    do k=1,nbin
       write(15,*) xi(:,k)
    end do
    close(15)

  endif
  sync all


#ifdef write_xreal
    close(12)
#endif

  enddo
  call destroy_penfft_plan
  sync all
  if (head) print*,'cicpower done'
  sync all
end

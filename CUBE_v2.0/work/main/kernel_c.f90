#define LRCKCORR
subroutine kernel_c
  use variables
  use pencil_fft
  implicit none
  save
  include 'fftw3.f'

  character(*),parameter :: dir_kern='../../kernels/'
  integer(8),parameter :: ncglobal=nc*nn
  integer(8) ig,jg,kg,itemp(3)
  real rx(3),r,ck_table(3,4,4,4),kx(3),kx_sin(3),kr

  if (head) print*, 'coarse kernel initialization'
  ! construct uncorrected force kernel
  ck=0
  do k=1,nc
  do j=1,nc
  do i=1,nc
    ig=i+nc*(icx-1)
    jg=j+nc*(icy-1)
    kg=k+nc*(icz-1)
    rx=ncell*(mod((/ig,jg,kg/)+ncglobal/2-1,ncglobal)-ncglobal/2)
    r=sqrt(rx(1)**2+rx(2)**2+rx(3)**2)
    ck(:,i,j,k)=merge(0.0,-rx/r**3,r==0)
  enddo
  enddo
  enddo

  ! Read in kernel for 2-level matching
  open(20,file=dir_kern//'wfxyzc.2.ascii',status='old')
  do k=1,4
  do j=1,4
  do i=1,4
    read(20,'(3i4,3e16.8)') itemp(1:3),ck_table(:,i,j,k)
  enddo
  enddo
  enddo
  close(20)

  ! copy corrections to other octants of the kernel
  ! need nc > 8
  if (icx==1 .and. icy==1 .and. icz==1) then
    ck(:,:4,:4,:4)=ck_table !1
  endif
  if (icx==nn .and. icy==1 .and. icz==1) then
    ck(1,nc-2:nc,:4,:4)=-ck_table(1,4:2:-1,:,:) !5 x
    ck(2:3,nc-2:nc,:4,:4)=ck_table(2:3,4:2:-1,:,:)
  endif
  if (icx==1 .and. icy==nn .and. icz==1) then
    ck(1:3:2,:4,nc-2:nc,:4)=ck_table(1:3:2,:,4:2:-1,:) !3 y
    ck(2,:4,nc-2:nc,:4)=-ck_table(2,:,4:2:-1,:)
  endif
  if (icx==1 .and. icy==1 .and. icz==nn) then
    ck(1:2,:4,:4,nc-2:nc)=ck_table(1:2,:,:,4:2:-1) !2 z
    ck(3,:4,:4,nc-2:nc)=-ck_table(3,:,:,4:2:-1)
  endif
  if (icx==1 .and. icy==nn .and. icz==nn) then
    ck(1,:4,nc-2:nc,nc-2:nc)=ck_table(1,:,4:2:-1,4:2:-1) !4 yz
    ck(2:3,:4,nc-2:nc,nc-2:nc)=-ck_table(2:3,:,4:2:-1,4:2:-1)
  endif
  if (icx==nn .and. icy==1 .and. icz==nn) then
    ck(1:3:2,nc-2:nc,:4,nc-2:nc)=-ck_table(1:3:2,4:2:-1,:,4:2:-1) !6 xz
    ck(2,nc-2:nc,:4,nc-2:nc)=ck_table(2,4:2:-1,:,4:2:-1)
  endif
  if (icx==nn .and. icy==nn .and. icz==1) then
    ck(:2,nc-2:nc,nc-2:nc,:4)=-ck_table(:2,4:2:-1,4:2:-1,:) !7 xy
    ck(3,nc-2:nc,nc-2:nc,:4)=ck_table(3,4:2:-1,4:2:-1,:)
  endif
  if (icx==nn .and. icy==nn .and. icz==nn) then
    ck(:,nc-2:nc,nc-2:nc,nc-2:nc)=-ck_table(:,4:2:-1,4:2:-1,4:2:-1) !8 xyz
  endif
  sync all
  if (head) print*, '  finished octants'

#ifdef LRCKCORR
  if (head) print*, '  LRCKCORR'
  force_c(:,1:nc,1:nc,1:nc)=ck ! back up ck
  ck=0
  do k=1,nc
  do j=1,nc
  do i=1,nc
    ig=i+nc*(icx-1)
    jg=j+nc*(icy-1)
    kg=k+nc*(icz-1)
    rx=ncell*(mod((/ig,jg,kg/)+ncglobal/2-1,ncglobal)-ncglobal/2)
    r=sqrt(rx(1)**2+rx(2)**2+rx(3)**2)
    ck(:,i,j,k)=merge(0.0,-rx/r**3,r==0)
  enddo
  enddo
  enddo

  do i_dim=1,3
    !if (head) print*,'correct dim',i_dim
    r3=force_c(i_dim,1:nc,1:nc,1:nc)
    call pencil_fft_forward
    kern_c(:,:,:,i_dim)=imag(cxyz)

    r3=ck(i_dim,:,:,:)
    call pencil_fft_forward

    do k=1,npen
    do j=1,ng
    do i=1,ng*nn/2+1
      kg=(nn*(icz-1)+icy-1)*npen+k
      jg=(icx-1)*ng+j
      ig=i
      kx=mod((/ig,jg,kg/)+ncglobal/2-1,ncglobal)-ncglobal/2
      kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
      kx_sin=2*sin(pi*kx/real(ncglobal))
      kern_c(i,j,k,i_dim)=merge(kern_c(i,j,k,i_dim), &
        kern_c(i,j,k,i_dim)*0.25*pi*kx_sin(i_dim)/sum(kx_sin**2)/imag(cxyz(i,j,k)), &
        (kr>8.0 .or. kx(i_dim)==0))
    enddo
    enddo
    enddo
  enddo
#else
  if (head) print*, '  without LRCKCORR'
  do i_dim=1,3
    r3=ck(i_dim,:,:,:)
    call pencil_fft_forward
    kern_c(:,:,:,i_dim)=imag(cxyz)
  enddo
#endif

  !if (head) print*, 'kernel_c done',kern_c(1,2,1,1)

  open(21,file=output_dir()//'kern_c'//output_suffix(),access='stream',status='replace')
  write(21) kern_c
  close(21)
  !print*, 'sum of kern_c', sum(kern_c)
  sync all
endsubroutine kernel_c

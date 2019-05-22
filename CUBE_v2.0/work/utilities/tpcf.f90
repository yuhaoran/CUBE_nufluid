! two point correlation function for vectors
program tpcf
!  use omp_lib
  use parameters
  use variables, only: spine_tile
  use iso_fortran_env, only : int64
  use halo_output, only: type_halo_catalog, type_halo_info
  implicit none
  save

  integer(4) t1,t2,tt1,tt2,ttt1,ttt2,t_rate,nhalo,ihalo,jhalo,imassbin,n_dist
  integer(8) kg,jg,ig,ii,jj
  real kr,kx,ky,kz,pow,r_filter,rm,masstemp
  real spin_u(3),spin_v(3),spin_r(3),spin_e(3,3),spin_vor(3)

  integer i,j,k,n_rsmall,n_ratio,nmassbin,itemp
  real t11,t22,t33,t12,t23,t31
  real tsmall(3,3),tlarge(3,3),torque(3,3)
  real dq,dx,corrq,corrx
  type(type_halo_info) halo_info
  type(type_halo_catalog) halo_catalog
  real,allocatable :: spin_x(:,:),spin_q(:,:),spin_t(:,:),xloc(:,:),qloc(:,:),imass(:),imass_info(:,:)
  real,allocatable :: corr_x(:,:,:),corr_q(:,:,:),corr_t(:,:,:),r_small(:),ratio_scale(:)
  real,allocatable :: tp(:,:)
  integer,allocatable :: isort_mass(:),i1(:),i2(:),ii1,ii2

  ! nc: coarse grid per image per dim
  ! nf: fine grid per image per dim, ng=nf
  ! nyquest: Nyquest frequency
  call omp_set_num_threads(ncore)
  call geometry
  open(16,file='../main/z_checkpoint.txt',status='old')
  do i=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  cur_checkpoint=n_checkpoint ! read halos first
  if (head) then
    print*, ''
    print*, 'program tpcf'
    print*, 'on',nn**3,' images and',ncore,' cores'
    print*, 'Resolution', ng*nn
    print*, 'at redshift=',z_checkpoint(cur_checkpoint)
    print*, 'Box size', box
    print*, 'output: ', opath
    print*, '-----------------------------------------'
  endif
  sync all
  ! open halo catalog
  print*,''
  print*,'read halo catalog'
  open(11,file=output_name('halo'),status='old',access='stream')
  open(12,file=output_name('halo_init_spin'),status='old',access='stream')
  open(13,file=output_name('halo_sort_index'),status='old',access='stream')
  read(11) halo_catalog
  nhalo=halo_catalog%nhalo
  allocate(isort_mass(nhalo),imass(nhalo))
  allocate(spin_x(3,nhalo),spin_q(3,nhalo),spin_t(3,nhalo),qloc(3,nhalo),xloc(3,nhalo))
  ! read halo q-pos & spin
  do ihalo=1,nhalo
    read(11) halo_info
    read(12) spin_q(:,ihalo),spin_t(:,ihalo),spin_x(:,ihalo),spin_u,spin_v,spin_r,spin_e(:,1:3)
    imass(ihalo)=halo_info%mass_odc
    qloc(:,ihalo)=halo_info%q_mean
    xloc(:,ihalo)=halo_info%x_mean
  enddo
  read(13) isort_mass ! read index by halo mass sort
  close(11)
  close(12)
  close(13)
  ! sort according to halo mass
  spin_x=spin_x(:,isort_mass)
  spin_q=spin_q(:,isort_mass)
  spin_t=spin_t(:,isort_mass)
  qloc=qloc(:,isort_mass)
  xloc=xloc(:,isort_mass)
  imass=imass(isort_mass)
  ! construct halo mass bins
  rm=2.0 ! mass bin ratio
  n_rsmall=50
  n_ratio=20
  n_dist=ceiling(ng*nn/2*sqrt(3.))+1
  nmassbin=ceiling(log(imass(1)/imass(nhalo))/log(rm))
  allocate(imass_info(4,nmassbin),i1(nmassbin),i2(nmassbin))
  allocate(corr_x(n_rsmall,n_ratio,nmassbin),corr_q(n_rsmall,n_ratio,nmassbin),corr_t(n_rsmall,n_ratio,nmassbin))
  allocate(r_small(n_rsmall),ratio_scale(n_ratio))
  allocate(tp(6,n_dist))
  print*,' N_halos_global =',halo_catalog%nhalo_tot
  print*,' N_halos_local  =',halo_catalog%nhalo
  print*,' Overdensity    =',halo_catalog%den_odc
  print*,' nmassbin =',nmassbin
  imassbin=nmassbin
  i2(nmassbin)=nhalo
  i1(1)=1
  itemp=1
  do ihalo=nhalo,1,-1
    if (imass(ihalo)<imass(nhalo)*rm**itemp .or. imassbin==1) then
      i1(imassbin)=ihalo ! still in the bin
    else
      imassbin=imassbin-1 ! assign to previous bin
      i2(imassbin)=ihalo
      itemp=itemp+1
    endif
  enddo

  tp=0
  ! correlation
  do ihalo=1,nhalo
  do jhalo=ihalo+1,nhalo
    dq=distance_periodic(qloc(:,ihalo),qloc(:,jhalo))
    dx=distance_periodic(xloc(:,ihalo),xloc(:,jhalo))
    corrq=ccc(spin_q(:,ihalo),spin_q(:,jhalo))
    corrx=ccc(spin_x(:,ihalo),spin_x(:,jhalo))
    tp(1,nint(dq)+1)=tp(1,nint(dq)+1)+1.0
    tp(2,nint(dx)+1)=tp(2,nint(dx)+1)+1.0
    tp(3,nint(dq)+1)=tp(3,nint(dq)+1)+corrq
    tp(4,nint(dq)+1)=tp(4,nint(dq)+1)+corrx
    tp(5,nint(dx)+1)=tp(5,nint(dx)+1)+corrq
    tp(6,nint(dx)+1)=tp(6,nint(dx)+1)+corrx
  enddo
  enddo

  do i=1,n_dist
    tp(3:4,i)=tp(3:4,i)/tp(1,i)
    tp(5:6,i)=tp(5:6,i)/tp(2,i)
  enddo

  open(11,file=output_name('tpcf'),status='replace',access='stream')
  write(11) int(6,kind=4), n_dist,  tp
  close(11)


contains

  real function ccc(vec1,vec2)
    implicit none
    real vec1(3),vec2(3)
    vec1=vec1/norm2(vec1)
    vec2=vec2/norm2(vec2)
    ccc=sum(vec1*vec2)
  endfunction

  real function distance_periodic(vec1,vec2)
    implicit none
    real vec1(3),vec2(3),dvec12(3)
    dvec12=vec2-vec1
    dvec12=modulo(dvec12+ng*nn/2,real(ng*nn))-ng*nn/2
    distance_periodic=norm2(dvec12)
  endfunction

endprogram

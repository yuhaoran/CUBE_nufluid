!! init.f90 written by Hy Trac
!! dist_init.f90 Parallelized: Hugh Merz Jun 4, 2005
!! Updated Jun 2 2006 -- minimal memory usage, DM only !!
!!-----------------------------------------------------!!
!! Compile with: mpif77 -fpp -g -w -O3 -axN dist_init.f90 -o dist_init  -L$MCKENZIE_FFTW_LIB_PATH -I$MCKENZIE_FFTW_INC_PATH -lsrfftw_mpi -lsrfftw -lsfftw_mpi -lsfftw -lm -ldl
!! On Ranger compile with:
!! mpif90 -Wl,-rpath, -L$TACC_MKL_LIB -I$TACC_MKL_INC -fpp -g -O2 -xW -shared-intel -DBINARY dist_init_dm.f90 -o dist_init -L$TACC_MKL_LIB -lmkl_em64t -openmp
!!mpif90 -shared-intel -mcmodel=medium -qopenmp -fpp -g -O3 -DVELTRANSFER -DNEUTRINOS dist_init_dmnu.f90 -I$SCINET_FFTW_INC -I$P3DFFT_INC -o dist_init_dmnu_nu -L$SCINET_FFTW_LIB -L$P3DFFT_LIB -lp3dfft -lfftw3f
program dist_init

  use omp_lib

#ifdef SFFTW3
  use, intrinsic :: iso_c_binding
#endif

  implicit none

  include 'mpif.h'
  include '../../parameters'
#ifdef SFFTW3
  include 'fftw3-mpi.f03'
  integer(C_INTPTR_T), parameter :: nct = nc
#endif

  integer, parameter  :: num_threads = cores*nested_threads
#ifdef NEUTRINOS
  logical, parameter :: generate_seeds=.false.
#else
  logical, parameter :: generate_seeds=.true.
#endif
  logical, parameter :: correct_kernel=.true.

  real, parameter :: ns = 0.9619
  real, parameter :: s8 = 0.8276
  real, parameter :: omegal=omega_l
  real, parameter :: omegam=1.0-omegal

#ifdef NEUTRINOS
  real, parameter :: redshift=z_i_nu
#else
  real, parameter :: redshift=z_i
#endif
  real, parameter :: scalefactor=1/(1+redshift)

#ifdef WDM
  real, parameter :: m_wdm = 1. ! keV
  real, parameter :: mu = 1.12
  real, parameter :: alpha = 0.049 * m_wdm**-1.11 * (omegam/0.25)**0.11 * (H0/100./0.7)**1.22
#endif

#ifdef BROKEN_POWER
  real, parameter :: dns = -0.05
  real, parameter :: k0  = 0.06
#endif

#ifdef NEUTRINOS
  integer, parameter :: nv=10000
  character(*), parameter :: cdfTable = 'CDFTable.txt'
  real(4), dimension(2,nv) :: cdf    !Col1 is v, col2 is cdf
#endif

  !! NEUTRINO sims may use this:
  integer, parameter      :: nk=1000
!fntf depends on redshift which could vary for neutrinos and dm
#ifdef NEUTRINOS
  character(*), parameter :: fntf = 'ith2_mnu0p05_z5_tk.dat'
#else
  character(*), parameter :: fntf = 'ith2_mnu0p05_z5_tk.dat'
#endif
#ifdef VELTRANSFER
#ifdef NEUTRINOS
  character(*), parameter :: vfntf = 'ith2_mnu0p05_z5_v_tk.dat'
#else
  character(*), parameter :: vfntf = 'ith2_mnu0p05_z5_v_tk.dat'
#endif
  real(4), parameter :: Vphys2sim = 1.0/(300. * sqrt(omega_m) * box * (1. + redshift) / 2. / nc)!(180.8892437/mass_neutrino)/(box*300.0*(omega_m)**0.5/2.0/nc)
#endif
  !! Transfer function file
!  integer, parameter      :: nk=922
!  character(*), parameter :: fntf='camb_jddefault_transfer_z0.dat'

  !! nc is the number of cells per box length
  integer, parameter :: hc=nc/2
  real, parameter    :: ncr=nc
  real, parameter    :: hcr=hc

  !! np is the number of particles
  !! np should be set to nc (1:1), hc (1:2), or qc (1:4)
!#ifdef NEUTRINOS
  integer, parameter :: np=nc
!#else
!  integer, parameter :: np=hc/ratio_nudm_dim
!#endif
  real, parameter    :: npr=np

  !! internal parallelization parameters
  integer(4), parameter :: nc_node_dim = nc/nodes_dim
  integer(4), parameter :: np_node_dim = np/nodes_dim
  integer(4), parameter :: nodes = nodes_dim * nodes_dim * nodes_dim
  integer(4), parameter :: nodes_slab = nodes_dim * nodes_dim

  !! parallelization variables
  integer(4), dimension(6) :: cart_neighbor
  integer(4), dimension(3) :: slab_coord, cart_coords
  integer(4) :: slab_rank, mpi_comm_cart, cart_rank, rank, ierr

  integer(4) :: np_local,wc_counter, count_i,count_f,count_r
  logical :: firstfftw
#ifdef SFFTW3
  type(C_PTR) :: plan, iplan
#else
  integer(8) :: plan, iplan
#endif

! :: simulation variables

  !! Other parameters
  real, parameter :: pi=3.141592654

  !! Power spectrum arrays
!  real, dimension(5,nk) :: tf   !cmbfast
  real, dimension(7,nk) :: tf    !CAMB
  real, dimension(2,nc) :: pkm,pkn

  !! For pencils decomposition:
#ifdef SLAB
  integer(4), dimension(0:nodes_dim-1,0:nodes_dim-1) :: slab_neighbor
  integer(4), parameter :: nc_slab = nc / nodes
#else
  integer(4), parameter   :: nodes_pen = nodes_dim
  integer(4), parameter   :: nc_pen = nc_node_dim / nodes_dim
  integer(4), parameter   :: dim_y = nodes_dim
  integer(4), parameter   :: dim_z = nodes_dim**2
  integer(4) :: pen_dims(2), istart(3), iend(3), isize(3), fstart(3), fend(3), fsize(3), mypadd
  integer(4), dimension(0:nodes_dim-1) :: pen_neighbor_to
  integer(4), dimension(0:nodes_dim-1) :: pen_neighbor_fm
#endif

  !! Fourier transform arrays
  real, dimension(nc_node_dim,nc_node_dim,nc_node_dim) :: cube
#ifdef SLAB
  real, dimension(nc+2,nc,nc_slab) :: slab, slab_work
  real, dimension(nc_node_dim,nc_node_dim,nc_slab,0:nodes_slab-1) :: recv_cube
#ifdef SFFTW3
  complex(4), dimension(nc+2,nc,nc_slab) :: slab_cmplx
#endif
#else
  real, dimension(nc, nc_node_dim, nc_pen+2) :: slab, slab_work
  real, dimension(nc_node_dim, nc_node_dim, nc_pen, 0:nodes_pen-1)      :: recv_cube
#endif
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1,0:nc_node_dim+1) :: phi
  real, dimension(0:nc_node_dim+1,0:nc_node_dim+1) :: phi_buf

  !! Particles arrays for subroutine dm
#ifdef ZIP
  real, parameter :: dbuf = 1.4
  integer, parameter :: mesh_scale = 4
  integer, parameter :: max_np = int(dbuf*np_node_dim)*np_node_dim**2
  integer, parameter :: np_buffer = max_np - np_node_dim**3
  real, dimension(6,np_node_dim,np_node_dim,np_node_dim) :: xvp_grid
  real, dimension(6,max_np) :: xvp
  real, dimension(6,np_buffer) :: xp_buf
  real, dimension(6*np_buffer) :: send_buf, recv_buf
  integer, parameter :: nm_node_dim = nc_node_dim / mesh_scale
  integer(4), dimension(nm_node_dim, nm_node_dim, nm_node_dim) :: hoc
  integer(4) :: ll(max_np)
  real(4), parameter :: dbuf_cell = 100
  integer(4), parameter :: max_cell_np = int(dbuf_cell*(hcr/ncr*mesh_scale)**3) + 1
  integer(4), dimension(nm_node_dim) :: rhoc_i4
  integer(1), dimension(3, max_cell_np, nm_node_dim) :: pos_i1
  integer(2), dimension(3, max_cell_np, nm_node_dim) :: vel_i2
#ifdef SLAB
  real, dimension(nc+2,nc,nc_slab) :: slab_cache
#else
  real, dimension(nc, nc_node_dim, nc_pen+2) :: slab_cache
#endif
#else
  real, dimension(6,np_node_dim,np_node_dim,num_threads) :: xvp
#endif

  !! Timing variables
  real(8) :: sec1, sec2

  !! Equivalence arrays to save memory
#ifdef ZIP
  equivalence (phi,slab_work,recv_cube,send_buf)
  equivalence (slab,cube,xp_buf,hoc)
  equivalence (slab_cache, ll, recv_buf)
  equivalence (xvp, xvp_grid)
#else
  equivalence (phi,slab_work,recv_cube)
#ifdef SFFTW3
  equivalence (slab_cmplx,slab)
#else
  equivalence (slab,cube)
#endif
#endif

  !! Common block
  common /rvar/ tf, pkm, pkn, phi_buf
  common / equiv1 / phi
#ifdef SFFTW3
  common / equiv2 / slab_cmplx
  common / cvar / cube
#else
  common / equiv2 / slab
#endif
#ifdef ZIP
  common / equiv3 / slab_cache
  common / zvar / rhoc_i4, pos_i1, vel_i2
#endif
  common / xvar / xvp

  call mpi_initialize

  sec1 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STARTING INITIAL CONDITIONS: ", sec1

  call omp_set_num_threads(num_threads)
  if (rank == 0) print*, 'Note: compile with openmp: num_threads=', num_threads

  if (rank == 0) call writeparams

  call initvar
  call transferfnc
  call noisemap
  call deltafield
  call potentialfield

#ifdef ZIP
# ifdef VELTRANSFER
    call dm_zip(0)
    call veltransfer
    call dm_zip(1)
# else
    call dm_zip
# endif
  call pass_particles
  call link_list
  call zip_checkpoint
#else !! ZIP
# ifdef VELTRANSFER
    call dm(0)
    call veltransfer
    call dm(1)
# else
    call dm
# endif
#endif

  if (rank == 0) call writepowerspectra

  call di_fftw(0)

  sec2 = mpi_wtime(ierr)
  if (rank == 0) write(*,*) "STOPPING INITIAL CONDITIONS: ", sec2
  if (rank == 0) write(*,*) "ELAPSED TIME: ", sec2-sec1
  call mpi_barrier(mpi_comm_world, ierr)

  call mpi_finalize(ierr)

contains

  subroutine mpi_initialize
    implicit none

    integer(4) :: i, j, nodes_returned
    integer(4) :: dims(3), ndim
    logical :: periodic(3), reorder

!! set up global mpi communicator

    call mpi_init(ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)
    if (nodes_returned /= nodes ) then
      write(*,*) 'dist_init compiled for a different number of nodes'
      write(*,*) 'mpirun nodes=',nodes_returned,'dist_init nodes=',nodes
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#ifdef SLAB
    if (mod(nc,nodes) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into slabs'
      write(*,*) 'nc=',nc,'nodes=',nodes,'mod(nc,nodes) != 0'
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#else
    if (mod(nc,nodes_dim**2) /= 0) then
      write(*,*) 'cannot evenly decompose mesh into pencils'
      write(*,*) 'nc=',nc, 'nodes_dim**2=',nodes_dim**2,'mod(nc,nodes_dim**2)=', mod(nc,nodes_dim**2)
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif
#endif
    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    if (rank==0) then
      write(*,*) 'dist_init running on',nodes,'nodes'
      write(*,*) 'using cubic distribution:',nodes_dim,'nodes per dimension'
      write(*,*) nc,'cells in mesh'
    endif

!! calculate coordinates within slab for cube processes

    slab_coord(3) = rank / nodes_slab
    slab_rank = rank - slab_coord(3) * nodes_slab
    slab_coord(2) = slab_rank / nodes_dim
    slab_coord(1) = slab_rank - slab_coord(2) * nodes_dim
    do j = 0, nodes_dim - 1
#ifdef SLAB
      do i = 0, nodes_dim - 1
        slab_neighbor(i,j) = i + j * nodes_dim + slab_coord(3) &
                           * nodes_slab
      enddo
#else
        pen_neighbor_to(j) = nodes_slab*slab_coord(3) + slab_coord(2) + j*nodes_dim
        pen_neighbor_fm(j) = nodes_slab*slab_coord(3) + j + nodes_dim*slab_coord(1)
#endif
    enddo

!! create cartesian communicator based on cubic decomposition

    dims(:) = nodes_dim
    periodic(:) = .true.
    reorder = .false.
    ndim = 3

    call mpi_cart_create(mpi_comm_world, ndim,dims, periodic, &
                       reorder, mpi_comm_cart, ierr)
    call mpi_comm_rank(mpi_comm_cart, cart_rank, ierr)
    call mpi_cart_coords(mpi_comm_cart, cart_rank, ndim,  &
                         cart_coords, ierr)

! cart_neighbor(1) -> down (negative z)
! cart_neighbor(2) -> up (positive z)
! cart_neighbor(3) -> back (negative y)
! cart_neighbor(4) -> front (positive y)
! cart_neighbor(5) -> left (negative x)
! cart_neighbor(6) -> right (positive x)

    do i = 0, ndim-1
      call mpi_cart_shift(mpi_comm_cart, i, 1, cart_neighbor(2*(i+1)-1), &
                          cart_neighbor(2*(i+1)), ierr)
    enddo

#ifdef DEBUG_LOW
  do i=0,nodes-1
    if (i==rank) write(*,'(8i5)') rank,cart_rank,cart_neighbor
    call mpi_barrier(mpi_comm_world,ierr)
  enddo
#endif

  end subroutine mpi_initialize

!-------------------------------------------------------------------!

#ifdef SLAB
subroutine pack_slab
!! pack cubic data into slab decomposition for fftw transform
    implicit none

    integer(4) :: i,j,k,i0,j0,i1,j1,k1
    integer(4) :: slab_slice,num_elements,tag,rtag

    integer(4), dimension(2*nodes_dim**2) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim**2) :: wait_status

    num_elements = nc_node_dim * nc_node_dim * nc_slab

!! swap data

    do j = 0, nodes_dim - 1
      do i = 0, nodes_dim - 1
        slab_slice = i + j * nodes_dim
        tag = rank**2
        rtag= slab_neighbor(i,j)**2
        call mpi_isend(cube(1,1,slab_slice*nc_slab + 1), num_elements, &
                       mpi_real, slab_neighbor(i,j), tag, mpi_comm_world, &
                       requests(slab_slice+1),ierr)
        call mpi_irecv(recv_cube(1,1,1,slab_slice), &
                       num_elements, mpi_real, slab_neighbor(i,j),rtag, &
                       mpi_comm_world, requests(slab_slice+1+nodes_dim**2), &
                       ierr)
      enddo
    enddo

    call mpi_waitall(2*nodes_dim**2, requests, wait_status, ierr)

!! place data in the slab

    do j = 0, nodes_dim - 1
      j0 = j * nc_node_dim + 1
      j1 = (j + 1) * nc_node_dim
      do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        slab_slice = i + j * nodes_dim
        slab(i0:i1,j0:j1,:) = recv_cube(:,:,:,slab_slice)
      enddo
    enddo

end subroutine pack_slab
#else
subroutine pack_pencils
    !
    ! Pack cubic data into pencils for p3dfft transform.
    !

    implicit none

    integer(4) :: i,j,k,i0,i1,k1
    integer(4) :: pen_slice,tag,rtag
    integer(8) :: num_elements_i8
    integer(4) :: num_elements
    integer(4), dimension(2*nodes_dim) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim) :: wait_status
    integer(4) nc_pen_break, breakup
    real(4) :: passGB

    !
    ! Ensure that send/recv buffers are no larger than 1 GB (really after 2 GB we get problems)
    !

    breakup = 1
    num_elements_i8 = int(nc_node_dim,kind=8) * nc_node_dim * nc_pen
    passGB = 4. * num_elements_i8 / 1024.**3
    if (passGB > 1.) then
        breakup = 2**ceiling(log(passGB)/log(2.))
    endif
    num_elements = num_elements_i8 / breakup

    !
    ! Send the data from cube to recv_cube
    !

    do k = 1, breakup
        nc_pen_break = nc_pen/breakup*(k-1)
        do j = 0, nodes_dim - 1
            pen_slice = j
            tag  = rank**2
            rtag = pen_neighbor_fm(j)**2
            call mpi_isend(cube(1,1, pen_slice*nc_pen + nc_pen_break + 1), num_elements, &
                           mpi_real, pen_neighbor_to(j), tag, mpi_comm_world, &
                           requests(pen_slice+1),ierr)
            call mpi_irecv(recv_cube(1,1,1+nc_pen_break,pen_slice), &
                           num_elements, mpi_real, pen_neighbor_fm(j),rtag, &
                           mpi_comm_world, requests(pen_slice+1+nodes_dim), &
                           ierr)
        enddo

        call mpi_waitall(2*nodes_dim, requests, wait_status, ierr)

    enddo

    !
    ! Place this data into the pencils (stored in the slab array)
    !

    do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        pen_slice = i

        do k = 1, nc_pen
            do j = 1, nc_node_dim
                slab(i0:i1,j,k) = recv_cube(:,j,k,pen_slice)
            enddo
        enddo

    enddo

end subroutine pack_pencils
#endif

!-------------------------------------------------------------------!

#ifdef SLAB
subroutine unpack_slab
!! unpack slab data into cubic decomposition following fftw transform
    implicit none

    integer(4) :: i,j,k,i0,j0,i1,j1,k1
    integer(4) :: slab_slice,num_elements,tag,rtag
    integer(4), dimension(2*nodes_dim**2) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim**2) :: wait_status

!! place data in the recv_cube buffer

    do j = 0, nodes_dim - 1
      j0 = j * nc_node_dim + 1
      j1 = (j + 1) * nc_node_dim
      do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        slab_slice = i + j * nodes_dim
        recv_cube(:,:,:,slab_slice) = slab(i0:i1,j0:j1,:)
      enddo
    enddo

    num_elements = nc_node_dim * nc_node_dim * nc_slab

!! swap data

   do j = 0, nodes_dim - 1
      do i = 0, nodes_dim - 1
        slab_slice = i + j * nodes_dim
        tag  = rank**2
        rtag = slab_neighbor(i,j)**2
        call mpi_isend(recv_cube(1,1,1,slab_slice), num_elements, &
                       mpi_real, slab_neighbor(i,j), tag, mpi_comm_world, &
                       requests(slab_slice+1),ierr)
        call mpi_irecv(cube(1,1,slab_slice * nc_slab +1), &
                       num_elements, mpi_real, slab_neighbor(i,j), rtag, &
                       mpi_comm_world, requests(slab_slice+1+nodes_dim**2), &
                       ierr)
      enddo
    enddo

    call mpi_waitall(2*nodes_dim**2,requests, wait_status, ierr)

end subroutine unpack_slab
#else
subroutine unpack_pencils
    !
    ! Unpack data from the pencils back into the cubic decompisition following
    ! p3dfft transform.
    !

    implicit none

    integer(4) :: i,j,k,i0,i1,k1
    integer(4) :: pen_slice,num_elements,tag,rtag
    integer(8) :: num_elements_i8
    integer(4), dimension(2*nodes_dim) :: requests
    integer(4), dimension(MPI_STATUS_SIZE,2*nodes_dim) :: wait_status
    integer(4) nc_pen_break, breakup
    real(4) :: passGB

    !
    ! Place data in the recv_cube buffer
    !

    do i = 0, nodes_dim - 1
        i0 = i * nc_node_dim + 1
        i1 = (i + 1) * nc_node_dim
        pen_slice = i
        do k = 1, nc_pen
            do j = 1, nc_node_dim
                recv_cube(:, j, k, pen_slice) = slab(i0:i1, j, k)
            enddo
        enddo
    enddo

    !
    ! Ensure that send/recv buffers are no larger than 1 GB (really after 2 GB we get problems)
    !

    breakup = 1
    num_elements_i8 = int(nc_node_dim,kind=8) * nc_node_dim * nc_pen
    passGB = 4. * num_elements_i8 / 1024.**3
    if (passGB > 1.) then
        breakup = 2**ceiling(log(passGB)/log(2.))
    endif
    num_elements = num_elements_i8 / breakup

    !
    ! Put this data back into cube
    !

    do k = 1, breakup
        nc_pen_break = nc_pen/breakup*(k-1)
        do j = 0, nodes_dim - 1
            pen_slice = j
            tag  = rank**2
            rtag = pen_neighbor_to(j)**2
            call mpi_isend(recv_cube(1,1,1+nc_pen_break,pen_slice), num_elements, &
                           mpi_real, pen_neighbor_fm(j), tag, mpi_comm_world, &
                           requests(pen_slice+1),ierr)
            call mpi_irecv(cube(1,1,pen_slice*nc_pen + nc_pen_break + 1), &
                           num_elements, mpi_real, pen_neighbor_to(j),rtag, &
                           mpi_comm_world, requests(pen_slice+1+nodes_dim), &
                           ierr)
        enddo

        call mpi_waitall(2*nodes_dim,requests, wait_status, ierr)

    enddo

end subroutine unpack_pencils
#endif

!-------------------------------------------------------------------!

subroutine di_fftw(command)
    !
    ! Calculate fftw transform using P3DFFT or FFTW slab decomposition (ifdef SLAB)
    ! 0 ends fftw subprogram, 1 starts forward fft, -1 starts backwards
    !

#ifndef SLAB
    use p3dfft
#endif

#ifdef SFFTW3
    use, intrinsic :: iso_c_binding
#endif

    implicit none
#ifdef SLAB
#ifndef SFFTW3
    include 'fftw_f77.i'
    integer(4), parameter :: order = FFTW_NORMAL_ORDER
#endif
#endif

    integer(4) :: i
    integer(4) :: command

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'starting fftw',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

! initialize plan variables for fftw

    if (firstfftw) then
#ifdef SLAB
#ifdef SFFTW3
      call fftwf_mpi_init
      plan  = fftwf_mpi_plan_dft_r2c_3d(nct, nct, nct, slab, slab_cmplx, mpi_comm_world, FFTW_ESTIMATE)
      iplan = fftwf_mpi_plan_dft_c2r_3d(nct, nct, nct, slab_cmplx, slab, mpi_comm_world, FFTW_ESTIMATE)
#else
      call rfftw3d_f77_mpi_create_plan(plan,mpi_comm_world,nc, nc,nc, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE)
      call rfftw3d_f77_mpi_create_plan(iplan,mpi_comm_world,nc, nc,nc, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE)
#endif
#else
        pen_dims = (/dim_y,dim_z/)
        call p3dfft_setup(pen_dims, nc, nc, nc, .true.)
        call p3dfft_get_dims(istart, iend, isize, 1, mypadd)
        call p3dfft_get_dims(fstart, fend, fsize, 2)
#endif
#ifdef DEBUG_LOW
      print *,'finished initialization of fftw',rank
#endif
      firstfftw=.false.
    endif

! giver

    if (command /= 0) then

!! call pack routine if we are going forward

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'starting pack',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

#ifdef SLAB
      if (command > 0) call pack_slab
#else
      if (command > 0) call pack_pencils
#endif

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'finished forward slab pack',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    if (command > 0) then
#ifdef SLAB
#ifdef SFFTW3
        call sfftw_execute_dft_r2c(plan, slab, slab_cmplx)
#else
        call rfftwnd_f77_mpi(plan,1,slab,slab_work,1,order)
#endif
#else
        call ftran_r2c(slab, slab, "fft")
#endif
    else
#ifdef SLAB
#ifdef SFFTW3
        call sfftw_execute_dft_c2r(iplan, slab_cmplx, slab)
#else
        call rfftwnd_f77_mpi(iplan,1,slab,slab_work,1,order)
#endif
#else
        call btran_c2r(slab, slab, "tff")
#endif
        slab=slab/(real(nc)*real(nc)*real(nc))
    endif

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (rank == i) print *,'finished fftw',rank
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

!! unpack the slab data

#ifdef SLAB
      if (command < 0) call unpack_slab
#else
      if (command < 0) call unpack_pencils
#endif

    else

! if command = 0 we delete the plans
#ifdef SLAB
#ifdef SFFTW3
        call sfftw_destroy_plan(plan)
        call sfftw_destroy_plan(iplan)
        call fftwf_mpi_cleanup
#else
        call rfftwnd_f77_mpi_destroy_plan(iplan)
        call rfftwnd_f77_mpi_destroy_plan(plan)
#endif
#else
        call p3dfft_clean
#endif

    endif

end subroutine di_fftw

!-------------------------------------------------------------------!

  subroutine writeparams
    implicit none

    real time1,time2
    call cpu_time(time1)

    write(*,*) 'nodes   ', nodes
    write(*,*) 'nc      ', nc
    write(*,*) 'np      ', np
    write(*,*)
    write(*,*) 'n_s      ',ns
    write(*,*) 'sigma_8  ',s8
    write(*,*) 'omega_m  ',omegam
    write(*,*) 'omega_l  ',omegal
    write(*,*)
    write(*,*) 'box      ',box
    write(*,*) 'redshift ',redshift
    write(*,*)

#ifdef MY_TRANSFER
    write(*,*) 'power index',power_index
#endif

#ifdef WDM
    write(*,*) "m_wdm = ", m_wdm
    write(*,*) "mu    = ", mu
    write(*,*) "alpha = ", alpha
#endif

#ifdef BROKEN_POWER
    write(*,*) "dns = ", dns
    write(*,*) "k0  = ", k0
#endif

#ifdef NEUTRINOS
    write(*,*) 'ratio_nudm_dim ', ratio_nudm_dim
    write(*,*) 'Vphys2sim ',Vphys2sim
#endif

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write params'
    return
  end subroutine writeparams

!-------------------------------------------------------------------!

  subroutine writepowerspectra
    implicit none
    integer      :: k
    real         :: kr
    character*180 :: fn

    real time1,time2
    call cpu_time(time1)

    !! Output power spectrum
    !! 1st column is k
    !! 2nd is matter p(k)
    !! 3rd is standard deviation
    !! 6th is noise p(k)
    !! 7th is noise standard deviation
#ifdef NEUTRINOS
    fn=output_path//'pk_nu.init'
#else
    fn=output_path//'pk.init'
#endif
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,hc+1
       kr=2*pi*(k-1)/box
       write(11,*) kr,pkm(:,k),pkn(:,k)
    enddo
    close(11)

    !! Output cmbfast power spectrum
#ifdef NEUTRINOS
    fn=output_path//'pk0_nu.init'
#else
    fn=output_path//'pk0.init'
#endif
    write(*,*) 'Writing ',fn
    open(11,file=fn,recl=500)
    do k=2,2*nc+1
       kr=2*pi*(k-1)/(2*box)
       write(11,*) kr,power(kr,1,2),power(kr,1,3)
    enddo
    close(11)

    call cpu_time(time2)
    time2=time2-time1
    write(*,"(f8.2,a)") time2,'  Called write power spectra'
    return
  end subroutine writepowerspectra

!!------------------------------------------------------------------!!

  subroutine transferfnc
    implicit none
    integer :: i,k
    real    :: kr,kmax,dummy
    real*8  :: v8

    real time1,time2
    call cpu_time(time1)

    if (rank == 0) then
      !! Get transfer function from CMBFAST or CAMB
      write(*,*) 'Reading ',fntf
      open(11,file=fntf)
!      read(11,*) tf
      do k=1,nk
         read(11,*) tf(1,k),tf(2,k),tf(3,k),tf(4,k),tf(5,k),tf(6,k),tf(7,k)
      end do
      close(11)

      !! Compute \Delta^2
      do k=1,nk
         kr     =tf(1,k)
         tf(2,k)=kr**(3+ns)*tf(2,k)**2/(2*pi**2)
         tf(3,k)=kr**(3+ns)*tf(3,k)**2/(2*pi**2)
         tf(6,k)=kr**(3+ns)*tf(6,k)**2/(2*pi**2)
      enddo

      !! Compute dk
      tf(4,1)=tf(1,2)/2
!print*, tf(4,1);stop
      do k=2,nk-1
         tf(4,k)=(tf(1,k+1)-tf(1,k-1))/2
      enddo
      tf(4,nk)=tf(1,nk)-tf(1,nk-1)
!print*, tf(4,nk);stop
      !! Compute variance in 8 h^-1 Mpc spheres
      v8=0
      kmax=2*pi*sqrt(3.)*hc/box
!print*,kmax;stop
      do k=1,nk
         if (tf(1,k) .gt. kmax) exit
         v8=v8+tf(2,k)*tophat(tf(1,k)*8)**2*tf(4,k)/tf(1,k)
      enddo
!print*,k,v8,s8,sum(tf*1d0);stop
      !! Normalize to \sigma_8
      !! Include growth factor
#ifdef NEUTRINOS
      !! Set dm transfer function = neutrino one
      !! This step makes sure that the neutrinos have the correct
      !! sigma8 normalization
      tf(2,:) = tf(6,:)*(s8**2/v8)*Dgrow(scalefactor)**2
#else
      tf(2:3,:)=tf(2:3,:)*(s8**2/v8)*Dgrow(scalefactor)**2
#endif
    endif
    call mpi_barrier(mpi_comm_world,ierr)

#ifdef DEBUG_LOW
    do i=0,nodes-1
      if (i==rank) print *,rank,'finished mpi_barrier',ierr
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call mpi_bcast(tf,7*nk,mpi_real,0,mpi_comm_world,ierr)

#ifdef DEBUG_LOW
    print *,rank,'finished mpi_bcast',ierr
#endif

    !! tf(1,i) stores k
    !! tf(2,i) stores \Delta^2_m
    !! tf(3,i) stores \Delta^2_b
    !! tf(4,i) stores dk

    call cpu_time(time2)
    time2=time2-time1
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called transfer fnc'
    return
  end subroutine transferfnc

!------------------------------------------------------------------!

  subroutine noisemap
    implicit none
!! for more than 1+1e6 processes this will have to be increased
    character(len=6) :: rank_s,image_s
    integer      :: i,j,k,seedsize,clock
    integer, dimension(8) :: values

#ifdef POWER_CHECK
    integer :: ig,jg,kg
#endif
    real         :: x,x1,x2
    character*180 :: fn
    integer(4),allocatable,dimension(:) :: iseed
    integer(4),allocatable,dimension(:) :: iseed_all

    real time1,time2
    call cpu_time(time1)

    write(rank_s,'(i6)') rank
    write(image_s,'(i6)') rank+1
    rank_s=adjustl(rank_s)
    image_s=adjustl(image_s)

#ifdef DEBUG_LOW
    print *,'rank:',rank,'starting noisemap'
#endif

    call random_seed()
#ifdef IBM
    call random_seed(generator=2)
#endif
    call random_seed(size=seedsize)

    allocate(iseed(seedsize))
    allocate(iseed_all(seedsize*nodes))

    call system_clock(count=clock)
    iseed = clock + 37 * (/ (i - 1, i = 1, 2) /)

    call date_and_time(values=values)
    if(rank==0) write(*,*) values(7:8), iseed

    call random_seed(put=values(7:8))
    call random_seed(put=iseed(1:seedsize))

    if (generate_seeds) then
#ifdef NODE_DEPENDANT
       if (rank == 0) then
         write(*,*) 'Generating seeds'
         do k=1,seedsize
            do j=0,rank
              call random_number(x)
            enddo
            iseed(k)=int(x*huge(0)+time1)+rank
         enddo
       endif
#else
       if (rank==0) then
         write(*,*) 'Generating seeds'
         do j=1,nodes
           do k=1,seedsize
             call random_number(x)
             iseed_all((j-1)*seedsize+k)=int(x*huge(0))
           enddo
         enddo
       endif
       call mpi_scatter(iseed_all,seedsize,mpi_integer,iseed,seedsize,mpi_integer,0,mpi_comm_world,ierr)
#endif
    else
       fn=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'seed'//rank_s(1:len_trim(rank_s))//'.init'
       if (rank == 0) print *, 'rank',rank,'Reading ',fn(1:len_trim(fn))
       open(11,file=fn)
       do k=1,seedsize
          read(11,*) i,iseed(k)
       enddo
       close(11)
    endif

    !! Generate random reals between 0 and 1
    call random_seed(put=iseed(1:seedsize))
    call random_number(cube)


open(11,file='../../CUBE/output/universe1/image'//image_s(1:len_trim(image_s))//'/'//&
'noise_'//image_s(1:len_trim(image_s))//'.bin',access='stream')
read(11) cube
print*, 'noise',rank,cube(1:4,1,1)
close(11)

print*, size(cube), 'read cube from CUBE'
call mpi_barrier(mpi_comm_world, ierr)

    fn=output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'seed'//rank_s(1:len_trim(rank_s))//'.init'
    if (rank == 0) print *, 'rank',rank,'Writing ',fn(1:len_trim(fn))
    open(11,file=fn)
    do k=1,seedsize
       write(11,*) k,iseed(k)
    enddo
    close(11)

    deallocate(iseed)
    deallocate(iseed_all)

    if (rank==0) then
       write(*,*) 'starting generation of Gaussian random numbers'
    end if

    !! Generate random Gaussian numbers
    !$omp parallel do default(shared) private(i,j,k,x1,x2)
    do k=1,nc_node_dim
       do j=1,nc_node_dim
          do i=1,nc_node_dim,2
             x1=2*pi*cube(i,j,k);! print*,x1
             !if (cube(i+1,j,k)<=0.0) cube(i+1,j,k)=0.0001
             x2=sqrt(-2*log(1-cube(i+1,j,k)));! print*,x2
             cube(i,j,k)=x2*cos(x1)
             cube(i+1,j,k)=x2*sin(x1)
!print*, cube(1:2,1,1); stop
          enddo
       enddo
    enddo
    !$omp end parallel do

    if (rank==0) then
       write(*,*) 'finished generation of Gaussian random numbers'
    end if

#ifdef POWER_CHECK
    do k=1,nc_node_dim
      kg=k+cart_coords(1)*nc_node_dim
      do j=1,nc_node_dim
        jg=j+cart_coords(2)*nc_node_dim
        do i=1,nc_node_dim
          ig=i+cart_coords(3)*nc_node_dim
          cube(i,j,k)=mod(ig,2)+mod(jg,4)+mod(kg,8)
        enddo
      enddo
    enddo
#endif

    if (rank==0) then
       write(*,*) 'starting forward FFT'
    end if

    !! Forward FFT white noise field
    call di_fftw(1)
    !! noise is now in slab array

    if (rank==0) then
       write(*,*) 'finished forward FFT'
    end if

    !! Generate noise spectrum
    call powerspectrum(pkn)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called noisemap'
    return
  end subroutine noisemap

!-----------------------------------------------------------------------!

  subroutine deltafield
    implicit none
    integer :: i,j,k,kg,mg,jg,ig,fstat
    real    :: kr,kx,ky,kz
    real    :: powb,powm
    real    :: d,dmin,dmax,dmint,dmaxt
    real*8  :: dsum,dvar,dsumt,dvart
    character(len=4) :: rank_string
    character(len=100) :: check_name

    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

    !! Interpolate \Delta^2 to grid
    !! Determine Fourier modes \delta(k)

#ifndef SLAB
    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
#endif

#ifdef SLAB
    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(kr,kx,ky,kz,powb,powm,kg)
    do k=1,nc_slab
       kg=k+nc_slab*rank
       if (kg .lt. hc+2) then
          kz=kg-1
       else
          kz=kg-1-nc
       endif
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          do i=1,nc+2,2
             kx=(i-1)/2
             kr=sqrt(kx**2+ky**2+kz**2)
             if (kr .eq. 0) then
                slab(i:i+1,j,k)=0.
             else
                powm=power(2*pi*kr/box,1,2)/(4*pi*kr**3)
                slab(i:i+1,j,k)=sqrt(powm*ncr**3)*slab(i:i+1,j,k)
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do
#else
    !! Cannot thread because of ind index
    do k = 1, nc_pen+mypadd
        do j = 1, nc_node_dim
            do i = 1, nc, 2
                kg = ind / dxy
                mg = ind - kg * dxy
                jg = mg / dx
                ig = mg - jg * dx
                kg = fstart(3) + kg
                jg = fstart(2) + jg
                ig = 2 * (fstart(1) + ig) - 1
                ind = ind + 1
                if (kg < hc+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc
                endif
                if (jg < hc+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc
                endif
                kx = (ig-1)/2
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if (kr .eq. 0) then
                    slab(i:i+1,j,k)=0.
                else
                    powm=power(2*pi*kr/box,1,2)/(4*pi*kr**3)
                    slab(i:i+1,j,k)=sqrt(powm*ncr**3)*slab(i:i+1,j,k)
                endif
            enddo
        enddo
    enddo
#endif

    !! Calculate matter delta power spectrum

    call powerspectrum(pkm)

    !! Calculate matter delta field statistics

    call di_fftw(-1)
  print*,'delta_L',cube(1:4,1,1)
  open(11,file='delta_L.dat',status='replace',access='stream')
  write(11) cube/Dgrow(scalefactor)
  close(11)
if (rank==0) then
  print*, 'omegam,omegal', omegam,omegal
  print*, 'Dgrow(scalefactor),scalefactor',Dgrow(scalefactor),scalefactor
endif
  !print*, 'delta_L',rank
  !print*, cube(1:4,1,1)
    dmin=0
    dmax=0
    dsum=0
    dvar=0

    !$omp parallel do default(shared) private(i,j,k,d) &
    !$omp& reduction(min:dmin) &
    !$omp& reduction(max:dmax) &
    !$omp& reduction(+:dsum,dvar)
    do k=1,nc_node_dim
       do j=1,nc_node_dim
          do i=1,nc_node_dim
             d=cube(i,j,k)
             dsum=dsum+d
             dvar=dvar+d*d
             dmin=min(dmin,d)
             dmax=max(dmax,d)
          enddo
       enddo
    enddo
    !$omp end parallel do

    call mpi_reduce(dsum,dsumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dvar,dvart,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    call mpi_reduce(dmin,dmint,1,mpi_real,mpi_min,0,mpi_comm_world,ierr)
    call mpi_reduce(dmax,dmaxt,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)

    if (rank==0) then
      dsum=dsumt/real(nc)**3
      dvar=sqrt(dvart/real(ncr)**3)
      write(*,*)
      write(*,*) 'Delta min    ',dmint
      write(*,*) 'Delta max    ',dmaxt
      write(*,*) 'Delta sum ',real(dsum)
      write(*,*) 'Delta var ',real(dvar)
      write(*,*)
    endif

#ifdef write_den
    if (rank == 0) then
        print *,'Writing density contrast to file'
    endif

    write(rank_string,'(i4)') rank
    rank_string = adjustl(rank_string)

#ifdef NEUTRINOS
    check_name = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'initdeltafield'// &
        rank_string(1:len_trim(rank_string))//'_nu.bin'
#else
    check_name = output_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'initdeltafield'// &
        rank_string(1:len_trim(rank_string))//'.bin'
#endif

    open(unit=21, file=check_name, status='replace', iostat=fstat, access='stream')

    if (fstat /= 0) then
        write(*,*) 'error opening density file'
        write(*,*) 'rank', rank, 'file:', check_name
        call mpi_abort(mpi_comm_world, ierr, ierr)
    endif

    do k = 1, nc_node_dim
        do j = 1, nc_node_dim

            write(21) cube(:, j, k)

        enddo
    enddo

#endif

    call di_fftw(1)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called delta field'
    return
  end subroutine deltafield

!!------------------------------------------------------------------!!

  subroutine potentialfield
    implicit none

    integer :: i,j,k,ioerr
    integer :: im,ip,ig,jm,jp,jg,km,kp,kg,mg
    real    :: r,x,y,z
    real    :: kr,ksq,kx,ky,kz
    real    :: phi8,phi8tot
    character*180 :: fn
!! for more than 1+1e6 processes this will have to be increased
    character(len=6) :: rank_s

    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

    if (correct_kernel) then
#ifdef ZIP
      slab_cache = slab
#else
      !! write delta to disk so we can build kernel
      if (rank == 0) write(*,*) 'Caching Delta on disk'
      write(rank_s,'(i6)') rank
      rank_s=adjustl(rank_s)
      fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'delta'//rank_s(1:len_trim(rank_s))

      open(11,file=fn,status='replace',iostat=ioerr, access='stream')

      if (ioerr /= 0) then
        print *,'error opening Delta cache file:',fn
        stop
      endif
#ifdef SLAB
      do k=1,nc_slab
#else
      do k=1,nc_pen+2
#endif
        write(11) slab(:,:,k)
      enddo
      close(11)
#endif
    else
      slab_work=slab
    endif

    !! Construct uncorrected potential kernel in Fourier space

#ifndef SLAB
    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
#endif

#ifdef SLAB
    !$omp parallel do default(shared) private(i,j,k,ip,kg) &
    !$omp& private(ksq,kx,ky,kz)
    do k=1,nc_slab
       kg=k+nc_slab*rank
       if (kg .lt. hc+2) then
          kz=kg-1
       else
          kz=kg-1-nc
       endif
       kz=2*sin(pi*kz/ncr)
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          ky=2*sin(pi*ky/ncr)
          do i=1,nc+2,2
             ip=i+1
             kx=(i-1)/2
             kx=2*sin(pi*kx/ncr)
             ksq=kx**2+ky**2+kz**2
             if (ksq .eq. 0) then
                slab(i:ip,j,k)=0
             else
                slab(i,j,k)=-4*pi/ksq
                slab(ip,j,k)=0
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do
#else
    !! Cannot thread because of ind index
    do k = 1, nc_pen+mypadd
        do j = 1, nc_node_dim
            do i = 1, nc, 2
                kg = ind / dxy
                mg = ind - kg * dxy
                jg = mg / dx
                ig = mg - jg * dx
                kg = fstart(3) + kg
                jg = fstart(2) + jg
                ig = 2 * (fstart(1) + ig) - 1
                ind = ind + 1
                if (kg < hc+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc
                endif
                kz=2*sin(pi*kz/ncr)
                if (jg < hc+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc
                endif
                ky=2*sin(pi*ky/ncr)
                kx = (ig-1)/2
                kx=2*sin(pi*kx/ncr)
                ksq=kx**2+ky**2+kz**2
                ip=i+1
                if (ksq .eq. 0) then
                    slab(i:ip,j,k)=0
                else
                    slab(i,j,k)=-4*pi/ksq
                    slab(ip,j,k)=0
                endif
            enddo
        enddo
    enddo
#endif
!print*,sum(slab*1d0),rank
    if (correct_kernel) then

       !! Inverse FFT potential kernel
       call di_fftw(-1)
print*, 'kernel',rank,sum(cube*1d0)
       phi8=0.0

       if (nc_node_dim < 9) then
         print *,'warning: mesh too small in potential kernel correction'
         call mpi_abort(mpi_comm_world,ierr,ierr)
       endif

       if (cart_coords(1) == 0 .and. cart_coords(2) == 0 .and. &
           cart_coords(3) == 0) phi8=cube(9,1,1)+cube(1,9,1)+cube(1,1,9)

       if (cart_coords(3) == nodes_dim-1 .and. cart_coords(2) == 0 .and. &
           cart_coords(1) == 0) phi8=phi8+cube(nc_node_dim-7,1,1)

       if (cart_coords(3) == 0 .and. cart_coords(2) == nodes_dim-1 .and. &
           cart_coords(1) == 0) phi8=phi8+cube(1,nc_node_dim-7,1)

       if (cart_coords(3) == 0 .and. cart_coords(2) == 0 .and. &
           cart_coords(1) == nodes_dim -1) phi8=phi8+cube(1,1,nc_node_dim-7)

       call mpi_reduce(phi8,phi8tot,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
       if (rank == 0) then
         phi8=phi8tot/6.0
         print*, 'phi8',phi8
       endif

       call mpi_bcast(phi8,1,mpi_real,0,mpi_comm_world,ierr)


       !! Construct Ewald potential kernel in real space
       !$omp parallel do default(shared) private(i,j,k) &
       !$omp& private(r,x,y,z,ig,jg,kg)
       do k=1,nc_node_dim
          kg=k+nc_node_dim*cart_coords(1)
          if (kg .lt. hc+2) then
             z=kg-1
          else
             z=kg-1-nc
          endif
          do j=1,nc_node_dim
             jg=j+nc_node_dim*cart_coords(2)
             if (jg .lt. hc+2) then
                y=jg-1
             else
                y=jg-1-nc
             endif
             do i=1,nc_node_dim
                ig=i+nc_node_dim*cart_coords(3)
                if (ig .lt. hc+2) then
                   x=ig-1
                else
                   x=ig-1-nc
                endif
                r=sqrt(x**2+y**2+z**2)
                if (r .gt. 8) then
                   cube(i,j,k)=cube(i,j,k)-(phi8+1/8.)
                else
                   if (r .eq. 0) then
                      cube(i,j,k)=-2.5
                   else
                      cube(i,j,k)=-1/r
                   endif
                endif
             enddo
          enddo
       enddo
       !$omp end parallel do

#ifdef DEBUG_KC
!       print *,'rank',rank,'phi8=',phi8
!       if (rank == 0) then
!         do i=1,16
!           print *,rank,cube(i,i,i)
!         enddo
!       endif
       print *,'rank',rank,'sum kern',sum(cube)
#endif
print*,'ewarld kernel',rank,sum(cube*1d0)

       !! Forward FFT potential kernel
       call di_fftw(1)


#ifdef DEBUG_KC
      phi8=0.0
      phi8tot=0.0
#ifdef SLAB
      do k=1,nc_slab
        do j=1,nc
          do i=1,nc+2
#else
      do k=1,nc_pen+2
        do j=1,nc_node_dim
          do i=1,nc
#endif
            phi8=phi8+slab(i,j,k)
          enddo
        enddo
      enddo
      call mpi_reduce(phi8,phi8tot,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
      if (rank==0) print *,'total slab=',phi8tot
      phi8=0.0
      phi8tot=0.0
#ifdef SLAB
      do k=1,nc_slab
        do j=1,nc
          do i=1,nc+2
#else
      do k=1,nc_pen+2
        do j=1,nc_node_dim
          do i=1,nc
#endif
            phi8=phi8+delta(i,j,k)
          enddo
        enddo
      enddo
      call mpi_reduce(phi8,phi8tot,1,mpi_real,mpi_sum,0,mpi_comm_world,ierr)
      if (rank==0) print *,'total delta=',phi8tot
#endif

#ifdef ZIP
      slab_work = slab_cache
#else
      open(11,file=fn,status='old',iostat=ioerr,access='stream')
      if (ioerr /= 0) then
        print *,'error opening Delta cache file:',fn
        stop
      endif
#ifdef SLAB
      do k=1,nc_slab
#else
      do k=1,nc_pen+2
#endif
        read(11) slab_work(:,:,k)
      enddo
      close(11)
#endif

    endif

    !! Complex multiply density field with potential kernel
    !$omp parallel do default(shared) private(i,j,k)
#ifdef SLAB
    do k=1,nc_slab
       do j=1,nc
          do i=1,nc+2,2
#else
    do k=1,nc_pen+2
       do j=1,nc_node_dim
          do i=1,nc,2
#endif
             slab(i:i+1,j,k)=slab_work(i:i+1,j,k)*slab(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    !! Inverse FFT potential field
    call di_fftw(-1)

    !! put cube in phi
    phi=0.0
    phi(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=cube

print*, 'phi',rank
print*, cube(1:4,1,1)

    call mesh_buffer

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called potential field'
    return
  end subroutine potentialfield

!!------------------------------------------------------------------!!

#ifdef VELTRANSFER
subroutine veltransfer
    !
    ! Replacing Fourier density field delta(k) with delta(k)*T_velocity(k)/T_density(k)
    ! where T_velocity(k) is the neutrino velocity transfer function.
    ! Store results in phi and pass along to dm subroutine with COMMAND == 1.
    !

    implicit none

    integer :: ii, k, j, i, l, kg, jg, ig, mg, ioerr
    real(4) :: kz, ky, kx, kr, interpVTF,interpMTF, kphys, w1, w2
    integer :: dx, dxy, ind
    real, dimension(7,nk) :: vtf

#ifdef NEUTRINOS
    integer, parameter :: nucol = 6 !6 for neutrinos, 2 for dm
#else
    integer, parameter :: nucol = 2
#endif
    character*180 :: fn
    character(len=6) :: rank_s

    real time1,time2
    call cpu_time(time1)

    !
    ! Read neutrino density transfer function
    !

    if (rank==0) then
        write(*,*) 'Reading ',fntf
        open(11,file=fntf)
        do k=1,nk
            read(11,*) tf(1,k),tf(2,k),tf(3,k),tf(4,k),tf(5,k),tf(6,k),tf(7,k)
        end do
        close(11)
    endif

    call mpi_bcast(tf, 7*nk, mpi_real, 0, mpi_comm_world, ierr)

    !
    ! Read neutrino velocity transfer function
    !

    if (rank == 0) then
        write(*,*) 'Reading ',vfntf
        open(11,file=vfntf)
        do k=1,nk
            read(11,*) vtf(1,k),vtf(2,k),vtf(3,k),vtf(4,k),vtf(5,k),vtf(6,k),vtf(7,k)
        end do
        close(11)
        !! Put vtf into simulation units
        do j=2,7
            do k=1, nk
                vtf(j,k) = vtf(j,k) * Vphys2sim
            enddo
        enddo
    endif

    call mpi_bcast(vtf, 7*nk, mpi_real, 0, mpi_comm_world, ierr)

    !
    ! Transform slab from real space to Fourier space (slab was last modified in potentialfield
    ! where it stored the Zel'dovich potential field delta(k)/k^2 in real space).
    !

    call di_fftw(1)

    !
    ! Multiply slab by T_vel(k)/T_den(k)*k^2 and then divide by k for finite difference.
    !

#ifndef SLAB
    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
#endif
#ifdef SLAB
    !$omp parallel default(shared) private(k, j, i, kg, kz, ky,kx, kphys, kr, interpMTF, interpVTF, w1, w2, ii)
    !$omp do schedule(dynamic)
    do k = 1, nc_slab
        kg=k+nc_slab*rank
        if (kg .lt. hc+2) then
            kz=kg-1
        else
            kz=kg-1-nc
        endif
        do j = 1, nc
            if (j .lt. hc+2) then
                ky=j-1
            else
                ky=j-1-nc
            endif
            do i = 1, nc+2, 2
                kx=(i-1)/2

                kphys = 2*pi*sqrt(kx**2+ky**2+kz**2)/box
                kx = 2*sin(pi*kx/ncr)
                ky = 2*sin(pi*ky/ncr)
                kz = 2*sin(pi*kz/ncr)

#else
    !$omp parallel default(shared) private(k, j, i, kg, mg, jg, ig, ind, kz, ky, kx, kphys, kr, interpMTF, interpVTF, w1, w2, ii)
    !$omp do schedule(dynamic)
    do k = 1, nc_pen+mypadd
        ind = (k-1)*nc_node_dim*nc/2
        do j = 1, nc_node_dim
            do i = 1, nc, 2
                kg = ind / dxy
                mg = ind - kg * dxy
                jg = mg / dx
                ig = mg - jg * dx
                kg = fstart(3) + kg
                jg = fstart(2) + jg
                ig = 2 * (fstart(1) + ig) - 1
                ind = ind + 1
                if (kg < hc+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc
                endif
                if (jg < hc+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc
                endif
                kx = (ig-1)/2

                kphys = 2*pi*sqrt(kx**2+ky**2+kz**2)/box
                kx = 2*sin(pi*kx/ncr)
                ky = 2*sin(pi*ky/ncr)
                kz = 2*sin(pi*kz/ncr)

#endif
                kr = sqrt(kx**2+ky**2+kz**2)

                !! Interpolate vt(kphys) and t(kphys)
                !! Funciton call turned out to be expensive ... put interpolation in here
                !interpVTF = linear_interpolate(kphys, vtf(1,:), vtf(nucol,:), nk)
                !interpMTF = linear_interpolate(kphys, tf(1, :), tf(nucol,:), nk)
                if (kphys <= vtf(1,1)) then
                    interpVTF = vtf(nucol, 1)
                else if (kphys >= vtf(1, nk)) then
                    interpVTF = vtf(nucol, nk)
                else
                    do ii = 1, nk-1
                        if (kphys <= vtf(1,ii+1)) then
                            w1 = vtf(1,ii+1) - kphys
                            w2 = kphys - vtf(1,ii)
                            interpVTF = (w1*vtf(nucol,ii)+w2*vtf(nucol,ii+1))/(vtf(1,ii+1)-vtf(1,ii))
                            exit
                        endif
                    enddo
                endif
                if (kphys <= tf(1,1)) then
                    interpMTF = tf(nucol, 1)
                else if (kphys >= tf(1, nk)) then
                    interpMTF = tf(nucol, nk)
                else
                    do ii = 1, nk-1
                        if (kphys <= tf(1,ii+1)) then
                            w1 = tf(1,ii+1) - kphys
                            w2 = kphys - tf(1,ii)
                            interpMTF = (w1*tf(nucol,ii)+w2*tf(nucol,ii+1))/(tf(1,ii+1)-tf(1,ii))
                            exit
                        endif
                    enddo
                endif

                !! Multiply slab by interpVTF/interpMTF*kr**2/kr
                slab(i:i+1, j, k) = slab(i:i+1, j, k) * interpVTF/interpMTF * kr * (-1.0) !* (-2.0) * pi / ncr
                if (kr .EQ. 0) slab(i:i+1,j,k) = 0.0
            enddo
        enddo
    enddo
    !$omp end do
    !$omp end parallel

    !
    ! Inverse FFT potential field
    !
    call di_fftw(-1)

    !
    ! Put cube in phi
    !

    phi=0.0
    phi(1:nc_node_dim,1:nc_node_dim,1:nc_node_dim)=cube
    call mesh_buffer

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called veltransfer'
    return

end subroutine veltransfer

!!------------------------------------------------------------------!!

real(4) function linear_interpolate(xi, x, y, n)
    !
    ! Linearlly interpolates (x, y) at (xi, yi).
    !

    implicit none

    real(4) :: xi
    integer(4) :: n
    real(4), dimension(n) :: x, y

    real :: w1, w2, yi
    integer :: i

    if (xi <= x(1)) then
        yi = y(1)
    else if (xi >= x(n)) then
        yi = y(n)
    else
        do i = 1, n-1
            if (xi >= x(i) .and. xi <= x(i+1)) then
                w1 = x(i+1) - xi
                w2 = xi - x(i)
                yi = (w1*y(i) + w2*y(i+1)) / (x(i+1)-x(i))
                exit
            endif
        enddo
    endif

    linear_interpolate = yi

end function linear_interpolate
#endif

!!------------------------------------------------------------------!!
#ifndef ZIP
    !! Dark matter data
    !! xvp(1:3) = x,y,z in grid units from 0 to nc
    !! xvp(4:6) = vx,vy,vz in km/s
#ifdef VELTRANSFER
  subroutine dm(COMMAND)
#else
  subroutine dm
#endif
    implicit none
    integer :: i,j,k,ioerr
    integer(8) :: pos_in,pos_out
    integer :: i1,j1,k1,lb,ub
    real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt,vf
    real*8  :: dsum,dvar,dsumt,dvart
    real, dimension(3) :: dis,x
    real*8, dimension(3) :: xav
    character*250 :: fn
!! for more than 1+1e6 processes this will have to be increased
    character(len=6) :: rank_s
    character(len=7) :: z_s
    integer(4), dimension(11) :: header_garbage
    integer(4) :: thread
    integer(4), parameter :: xsize1 = 4*3*np_node_dim**2
    integer(4), parameter :: xsize2 = 2*xsize1
#ifdef NEUTRINOS
    real(4), parameter :: mfac = 180.8892437/mass_neutrino/scalefactor/3.0**0.5
    real(4), parameter :: ffac = 50.2476/mass_neutrino/scalefactor !ckT/m =50.25
#endif

#ifdef VELTRANSFER
    integer :: COMMAND
#endif
#ifdef NEUTRINOS
    integer :: l
    real(4) :: rnum1,rnum2,rnum3,rnum4
#endif

    real time1,time2
    call cpu_time(time1)

    write(z_s,'(f7.3)') redshift
    z_s=adjustl(z_s)

    write(rank_s,'(i6)') rank
    rank_s=adjustl(rank_s)


!! Open input file
#ifdef VELTRANSFER
    if (COMMAND == 1) then
       fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'delta'//rank_s(1:len_trim(rank_s))
       open(unit=31,file=fn,status='old',iostat=ioerr,access='stream')
    endif
#endif

!! Open output file
#ifdef VELTRANSFER
    if (COMMAND == 0) then !! Temp file
       if (rank == 0) write(*,*) 'Caching xvp on disk'
       fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//'delta'//rank_s(1:len_trim(rank_s))
       open(unit=11,file=fn,status='replace',iostat=ioerr,access='stream')
    else !! xv file ! COMMAND==1
# ifdef NEUTRINOS
       fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'_nu.dat'
# else
       fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'
# endif
       open(unit=11,file=fn,status='replace',iostat=ioerr,access='stream')
    endif ! COMMAND
#else
 !! xv file
#ifdef NEUTRINOS
    fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'_nu.dat'
#else
    fn=scratch_path//'/node'//rank_s(1:len_trim(rank_s))//'/'//z_s(1:len_trim(z_s))//'xv'//rank_s(1:len_trim(rank_s))//'.dat'
#endif
    open(unit=11,file=fn,status='replace',iostat=ioerr,access='stream')
#endif

#ifdef NEUTRINOS
    !! Read cdfTable
    if (rank == 0) then
        write(*,*) 'Reading ',cdfTable
        open(12, file=cdfTable)
        do k = 1, nv
        read(12,*) cdf(1,k), cdf(2,k)
        enddo
        close(12)
    endif
    call mpi_bcast(cdf, 2*nv, mpi_real, 0, mpi_comm_world, ierr)
#endif

    !! Write the header in checkpoint format with zeros for all variables (except np_local)
    np_local=np_node_dim**3
    header_garbage(:) = 0

#ifdef VELTRANSFER
    if (COMMAND==1) then
#endif
       !write(11) np_local, header_garbage
       write(11) np_local,scalefactor,0.0,-3./sqrt(scalefactor),0,1000.,1000.,1000.,1,1,1,1.
#ifdef VELTRANSFER
    endif
#endif
    !! Displace particles
    !! Finite-difference potential to get displacement field
    vf=vfactor(scalefactor)

#ifndef NEUTRINOS
    !$omp parallel default(shared) private(k,k1,j,j1,i,i1,pos_in,pos_out,thread)
#else
    !$omp parallel default(shared) private(k,k1,j,j1,i,i1,pos_in,pos_out,l,rnum1,rnum2,rnum3,rnum4,thread)
#endif
    thread = 1
    thread = omp_get_thread_num() + 1
#ifdef VELTRANSFER
    if (COMMAND == 0) then
        !! First time we call dm We only need the positions
        !$omp do schedule(dynamic)
        do k=1,np_node_dim
            k1=(nc/np)*(k-1)+1
            pos_out = 1 + xsize1*(k-1)
            do j=1,np_node_dim
                j1=(nc/np)*(j-1)+1
                do i=1,np_node_dim
                    i1=(nc/np)*(i-1)+1

                    xvp(1,i,j,thread)=(phi(i1-1,j1,k1)-phi(i1+1,j1,k1))/2./(4.*pi)+(i1-0.5)
                    xvp(2,i,j,thread)=(phi(i1,j1-1,k1)-phi(i1,j1+1,k1))/2./(4.*pi)+(j1-0.5)
                    xvp(3,i,j,thread)=(phi(i1,j1,k1-1)-phi(i1,j1,k1+1))/2./(4.*pi)+(k1-0.5)

                enddo
            enddo
            write(unit=11,pos=pos_out) xvp(1:3,:,:,thread) !! save in temp file
        enddo
        !$omp end do
    else
        !! Second time we call dm we want velocity
# ifdef NEUTRINOS
        !$omp do schedule(static,1) ordered
# else
        !$omp do schedule(dynamic)
# endif
        do k=1,np_node_dim
            k1=(nc/np)*(k-1)+1
            pos_out = 1 + sizeof(np_local) + sizeof(header_garbage) + xsize2*(k-1)
            pos_in = 1 + xsize1*(k-1)
            read(unit=31,pos=pos_in) xvp(1:3,:,:,thread) !! catch from temp file
# ifdef NEUTRINOS
            !$omp ordered
            call random_number(xvp(4:6,:,:,thread))
            !$omp end ordered
# endif
            do j=1,np_node_dim
                j1=(nc/np)*(j-1)+1
                do i=1,np_node_dim
                    i1=(nc/np)*(i-1)+1

# ifdef NEUTRINOS
                    rnum1 = xvp(4,i,j,thread)
                    rnum2 = xvp(5,i,j,thread)
                    rnum3 = xvp(6,i,j,thread)
# endif

                    xvp(4,i,j,thread)=(phi(i1-1,j1,k1)-phi(i1+1,j1,k1))/2./(4.*pi)
                    xvp(5,i,j,thread)=(phi(i1,j1-1,k1)-phi(i1,j1+1,k1))/2./(4.*pi)
                    xvp(6,i,j,thread)=(phi(i1,j1,k1-1)-phi(i1,j1,k1+1))/2./(4.*pi)

# ifdef NEUTRINOS
                    !! Convert to fermi-dirac
                    do l=1,nv
                       if (cdf(2,l).GT.rnum1) then
                          rnum4 = cdf(2,l)-rnum1
                          if ( (l.NE.1) .AND. (rnum1-cdf(2,l-1) .LT. rnum4) ) then
                             rnum1 = cdf(1,l-1)
                          else
                             rnum1 = cdf(1,l)
                          endif
                          exit
                       endif
                    enddo

                    !! Convert to fermi-dirac with units rnum1<-rnum1*c*(kT/m)
                    rnum1 = rnum1 * ffac
                    rnum2 = rnum2*2.0-1.0 !convert to range 1,-1
                    rnum3 = rnum3*2.0*pi  !convert to range 0,2pi

                    xvp(4,i,j,thread)=xvp(4,i,j,thread) + rnum1*(1.0-rnum2**2)**0.5*cos(rnum3)*Vphys2sim
                    xvp(5,i,j,thread)=xvp(5,i,j,thread) + rnum1*(1.0-rnum2**2)**0.5*sin(rnum3)*Vphys2sim
                    xvp(6,i,j,thread)=xvp(6,i,j,thread) + rnum1*rnum2*Vphys2sim
# endif
                enddo
            enddo
            write(unit=11,pos=pos_out) xvp(:,:,:,thread)
        enddo
        !$omp end do
    endif


#else
! if not define VELTRANSFER


# ifdef NEUTRINOS
    !$omp do schedule(static,1) ordered
# else
    !$omp do schedule(dynamic)
# endif
    do k=1,np_node_dim
        k1=(nc/np)*(k-1)+1
        pos_out = 1 + sizeof(np_local) + sizeof(header_garbage) + xsize2*(k-1)
# ifdef NEUTRINOS
        !$omp ordered
        call random_number(xvp(4:6,:,:,thread))
        !$omp end ordered
# endif
        do j=1,np_node_dim
            j1=(nc/np)*(j-1)+1
            do i=1,np_node_dim
                i1=(nc/np)*(i-1)+1

# ifdef NEUTRINOS
                 rnum1 = xvp(4,i,j,thread)
                 rnum2 = xvp(5,i,j,thread)
                 rnum3 = xvp(6,i,j,thread)
# endif

                 xvp(4,i,j,thread)=(phi(i1-1,j1,k1)-phi(i1+1,j1,k1))/2./(4.*pi)
                 xvp(5,i,j,thread)=(phi(i1,j1-1,k1)-phi(i1,j1+1,k1))/2./(4.*pi)
                 xvp(6,i,j,thread)=(phi(i1,j1,k1-1)-phi(i1,j1,k1+1))/2./(4.*pi)
                 xvp(1,i,j,thread)=xvp(4,i,j,thread)+(i1-0.5)
                 xvp(2,i,j,thread)=xvp(5,i,j,thread)+(j1-0.5)
                 xvp(3,i,j,thread)=xvp(6,i,j,thread)+(k1-0.5)
                 xvp(4,i,j,thread)=xvp(4,i,j,thread)*vf
                 xvp(5,i,j,thread)=xvp(5,i,j,thread)*vf
                 xvp(6,i,j,thread)=xvp(6,i,j,thread)*vf
# ifdef NEUTRINOS
                 !! Convert to fermi-dirac
                 do l=1,nv
                    if (cdf(2,l).GT.rnum1) then
                       rnum4 = cdf(2,l)-rnum1
                       if ( (l.NE.1) .AND. (rnum1-cdf(2,l-1) .LT. rnum4) ) then
                          rnum1 = cdf(1,l-1)
                       else
                          rnum1 = cdf(1,l)
                       endif
                       exit
                    endif
                 enddo

                 !! Convert to fermi-dirac with units rnum1<-rnum1*c*(kT/m)
                 rnum1 = rnum1 * ffac
                 rnum2 = rnum2*2.0-1.0 !convert to range 1,-1
                 rnum3 = rnum3*2.0*pi  !convert to range 0,2pi
                 xvp(4,i,j,thread)=xvp(4,i,j,thread) + rnum1*(1.0-rnum2**2)**0.5*cos(rnum3)*Vphys2sim
                 xvp(5,i,j,thread)=xvp(5,i,j,thread) + rnum1*(1.0-rnum2**2)**0.5*sin(rnum3)*Vphys2sim
                 xvp(6,i,j,thread)=xvp(6,i,j,thread) + rnum1*rnum2*Vphys2sim
# endif

            enddo
        enddo
        write(unit=11,pos=pos_out) xvp(:,:,:,thread)
    enddo
    !$omp end do
#endif
    !$omp end parallel

    !! Close I/O files
    close(11)
#ifdef VELTRANSFER
    if (command == 1) close(31)
#endif
    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return
  end subroutine dm
!!------------------------------------------------------------------!!
#else
#ifdef VELTRANSFER
  subroutine dm_zip(COMMAND)
#else
  subroutine dm_zip
#endif
    implicit none
    integer :: i,j,k,ioerr
    integer(8) :: pos_out
    integer :: i1,j1,k1,lb,ub
    real    :: d,dmin,dmax,sum_dm,sum_dm_local,dmint,dmaxt,vf
    real*8  :: dsum,dvar,dsumt,dvart
    real, dimension(3) :: dis,x
    real*8, dimension(3) :: xav
    character*250 :: fn
!! for more than 1+1e6 processes this will have to be increased
    character(len=6) :: rank_s
    character(len=7) :: z_s
    integer(4), dimension(11) :: header_garbage
    integer(4), parameter :: xsize1 = 4*3*np_node_dim**2
    integer(4), parameter :: xsize2 = 2*xsize1
#ifdef NEUTRINOS
    real(4), parameter :: mfac = 180.8892437/mass_neutrino/scalefactor/3.0**0.5
    real(4), parameter :: ffac = 50.2476/mass_neutrino/scalefactor !ckT/m =50.25
#endif

#ifdef VELTRANSFER
    integer :: COMMAND
#endif
#ifdef NEUTRINOS
    integer :: l
    real(4) :: rnum1,rnum2,rnum3,rnum4
#endif

    real time1,time2
    call cpu_time(time1)

#ifdef NEUTRINOS
    !! Read cdfTable
    if (rank == 0) then
        write(*,*) 'Reading ',cdfTable
        open(12, file=cdfTable)
        do k = 1, nv
        read(12,*) cdf(1,k), cdf(2,k)
        enddo
        close(12)
    endif
    call mpi_bcast(cdf, 2*nv, mpi_real, 0, mpi_comm_world, ierr)
#endif

    !! Displace particles
    !! Finite-difference potential to get displacement field
    vf=vfactor(scalefactor)

#ifndef NEUTRINOS
    !$omp parallel default(shared) private(k,k1,j,j1,i,i1)
#else
    !$omp parallel default(shared) private(k,k1,j,j1,i,i1,l,rnum1,rnum2,rnum3,rnum4)
#endif
#ifdef VELTRANSFER
    if (COMMAND == 0) then
        !! First time we call dm We only need the positions
        !$omp do schedule(dynamic)
        do k=1,np_node_dim
            k1=(nc/np)*(k-1)+1
            do j=1,np_node_dim
                j1=(nc/np)*(j-1)+1
                do i=1,np_node_dim
                    i1=(nc/np)*(i-1)+1
                    xvp_grid(1,i,j,k)=(phi(i1-1,j1,k1)-phi(i1+1,j1,k1))/2./(4.*pi)+(i1-0.5)
                    xvp_grid(2,i,j,k)=(phi(i1,j1-1,k1)-phi(i1,j1+1,k1))/2./(4.*pi)+(j1-0.5)
                    xvp_grid(3,i,j,k)=(phi(i1,j1,k1-1)-phi(i1,j1,k1+1))/2./(4.*pi)+(k1-0.5)
                enddo
            enddo
        enddo
        !$omp end do
    else
        !! Second time we call dm we want velocity
#ifdef NEUTRINOS
        !$omp do schedule(static,1) ordered
#else
        !$omp do schedule(dynamic)
#endif
        do k=1,np_node_dim
            k1=(nc/np)*(k-1)+1
#ifdef NEUTRINOS
            !$omp ordered
            call random_number(xvp_grid(4:6,:,:,k))
            !$omp end ordered
#endif
            do j=1,np_node_dim
                j1=(nc/np)*(j-1)+1
                do i=1,np_node_dim
                    i1=(nc/np)*(i-1)+1
#ifdef NEUTRINOS
                    rnum1 = xvp_grid(4,i,j,k)
                    rnum2 = xvp_grid(5,i,j,k)
                    rnum3 = xvp_grid(6,i,j,k)
#endif
                    xvp_grid(4,i,j,k)=(phi(i1-1,j1,k1)-phi(i1+1,j1,k1))/2./(4.*pi)
                    xvp_grid(5,i,j,k)=(phi(i1,j1-1,k1)-phi(i1,j1+1,k1))/2./(4.*pi)
                    xvp_grid(6,i,j,k)=(phi(i1,j1,k1-1)-phi(i1,j1,k1+1))/2./(4.*pi)
#ifdef NEUTRINOS
                    !! Convert to fermi-dirac
                    do l=1,nv
                       if (cdf(2,l).GT.rnum1) then
                          rnum4 = cdf(2,l)-rnum1
                          if ( (l.NE.1) .AND. (rnum1-cdf(2,l-1) .LT. rnum4) ) then
                             rnum1 = cdf(1,l-1)
                          else
                             rnum1 = cdf(1,l)
                          endif
                          exit
                       endif
                    enddo

                    !! Convert to fermi-dirac with units rnum1<-rnum1*c*(kT/m)
                    rnum1 = rnum1 * ffac
                    rnum2 = rnum2*2.0-1.0 !convert to range 1,-1
                    rnum3 = rnum3*2.0*pi  !convert to range 0,2pi

                    xvp_grid(4,i,j,k)=xvp_grid(4,i,j,k) + rnum1*(1.0-rnum2**2)**0.5*cos(rnum3)*Vphys2sim
                    xvp_grid(5,i,j,k)=xvp_grid(5,i,j,k) + rnum1*(1.0-rnum2**2)**0.5*sin(rnum3)*Vphys2sim
                    xvp_grid(6,i,j,k)=xvp_grid(6,i,j,k) + rnum1*rnum2*Vphys2sim
#endif
                enddo
            enddo
        enddo
        !$omp end do
    endif
#else
#ifdef NEUTRINOS
    !$omp do schedule(static,1) ordered
#else
    !$omp do schedule(dynamic)
#endif
    do k=1,np_node_dim
        k1=(nc/np)*(k-1)+1
        pos_out = 1 + sizeof(np_local) + sizeof(header_garbage) + xsize2*(k-1)
#ifdef NEUTRINOS
        !$omp ordered
        call random_number(xvp_grid(4:6,:,:,k))
        !$omp end ordered
#endif
        do j=1,np_node_dim
            j1=(nc/np)*(j-1)+1
            do i=1,np_node_dim
                i1=(nc/np)*(i-1)+1
#ifdef NEUTRINOS
                 rnum1 = xvp_grid(4,i,j,k)
                 rnum2 = xvp_grid(5,i,j,k)
                 rnum3 = xvp_grid(6,i,j,k)
#endif
                 xvp_grid(4,i,j,k)=(phi(i1-1,j1,k1)-phi(i1+1,j1,k1))/2./(4.*pi)
                 xvp_grid(5,i,j,k)=(phi(i1,j1-1,k1)-phi(i1,j1+1,k1))/2./(4.*pi)
                 xvp_grid(6,i,j,k)=(phi(i1,j1,k1-1)-phi(i1,j1,k1+1))/2./(4.*pi)
                 xvp_grid(1,i,j,k)=xvp_grid(4,i,j,k)+(i1-0.5)
                 xvp_grid(2,i,j,k)=xvp_grid(5,i,j,k)+(j1-0.5)
                 xvp_grid(3,i,j,k)=xvp_grid(6,i,j,k)+(k1-0.5)
                 xvp_grid(4,i,j,k)=xvp_grid(4,i,j,k)*vf
                 xvp_grid(5,i,j,k)=xvp_grid(5,i,j,k)*vf
                 xvp_grid(6,i,j,k)=xvp_grid(6,i,j,k)*vf
#ifdef NEUTRINOS
                 !! Convert to fermi-dirac
                 do l=1,nv
                    if (cdf(2,l).GT.rnum1) then
                       rnum4 = cdf(2,l)-rnum1
                       if ( (l.NE.1) .AND. (rnum1-cdf(2,l-1) .LT. rnum4) ) then
                          rnum1 = cdf(1,l-1)
                       else
                          rnum1 = cdf(1,l)
                       endif
                       exit
                    endif
                 enddo

                 !! Convert to fermi-dirac with units rnum1<-rnum1*c*(kT/m)
                 rnum1 = rnum1 * ffac
                 rnum2 = rnum2*2.0-1.0 !convert to range 1,-1
                 rnum3 = rnum3*2.0*pi  !convert to range 0,2pi
                 xvp_grid(4,i,j,k)=xvp_grid(4,i,j,k) + rnum1*(1.0-rnum2**2)**0.5*cos(rnum3)*Vphys2sim
                 xvp_grid(5,i,j,k)=xvp_grid(5,i,j,k) + rnum1*(1.0-rnum2**2)**0.5*sin(rnum3)*Vphys2sim
                 xvp_grid(6,i,j,k)=xvp_grid(6,i,j,k) + rnum1*rnum2*Vphys2sim
#endif
            enddo
        enddo
    enddo
    !$omp end do
#endif
    !$omp end parallel

    call cpu_time(time2)
    time2=(time2-time1)

    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called dm'
    return

end subroutine dm_zip

!!--------------------------------------------------------------!!

subroutine pass_particles
    !
    ! Pass particles that may have been displaced out of their local volume to
    ! their corresponding node.
    !

    implicit none

    integer i,pp,np_buf,np_exit,npo,npi
    integer(8) :: np_final
    real x(6),lb,ub
    integer, dimension(mpi_status_size) :: status,sstatus,rstatus
    integer :: tag,srequest,rrequest,sierr,rierr
    real(4), parameter :: eps = 1.0e-03

    real(8) time1, time2
    time1 = mpi_wtime(ierr)

    np_local = np_node_dim**3

    !
    ! Identify particles within the buffer
    !

    lb = 0.
    ub = real(nc_node_dim)

    np_buf = 0
    pp = 1

    do

        if (pp > np_local) exit

        x = xvp(:, pp)

        !! See if it lies within the buffer
        if (x(1) < lb .or. x(1) >= ub .or. x(2) < lb .or. x(2) >= ub .or. &
            x(3) < lb .or. x(3) >= ub ) then

            !! Make sure we aren't exceeding the maximum
            np_buf = np_buf + 1

            if (np_buf > np_buffer) then
                print *, rank, 'np_buffer =', np_buffer, 'exceeded - np_buf =', np_buf
                call mpi_abort(mpi_comm_world, ierr, ierr)
            endif

            xp_buf(:, np_buf) = xvp(:, pp)
            xvp(:, pp)        = xvp(:, np_local)
            np_local          = np_local - 1

            cycle

        endif

        pp = pp + 1

    enddo

    call mpi_reduce(np_buf, np_exit, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)

    if (rank == 0) print *,'Total exiting particles = ',np_exit

    !
    ! Pass +x
    !

    !! Find particles that need to be passed

    tag = 11
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(1, pp) >= ub) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

#ifdef DEBUG
    do i = 0, nodes-1
      if (rank == i) print *, rank, 'np_out=', npo
      call mpi_barrier(mpi_comm_world, ierr)
    enddo
#endif

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(6), &
                              tag, cart_neighbor(5), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf + pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(1, np_buf + pp) = max(xp_buf(1, np_buf+pp) - ub, lb)
    enddo

    pp = 1

    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf + pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local + 1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi - 1
        cycle
      endif
      pp = pp + 1
    enddo

    np_buf = np_buf + npi

    !
    ! Pass -x
    !

    tag = 12
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(1, pp) < lb) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(5), &
                              tag, cart_neighbor(6), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(5), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(6), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf+pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(1, np_buf+pp) = min(xp_buf(1,np_buf+pp) + ub, ub-eps)
    enddo

    pp = 1
    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local + 1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi - 1
        cycle
      endif
      pp = pp + 1
    enddo

    np_buf = np_buf + npi

    !
    ! Pass +y
    !

    tag = 13
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(2, pp) >= ub) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(4), &
                              tag, cart_neighbor(3), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf+pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(2, np_buf+pp) = max(xp_buf(2, np_buf+pp)-ub, lb)
    enddo

    pp = 1
    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local + 1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi = npi-1
        cycle
      endif
      pp = pp + 1
    enddo

    np_buf = np_buf + npi

    !
    ! Pass -y
    !

    tag = 14
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(2,pp) < lb) then
        npo = npo+1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi, 1, mpi_integer, cart_neighbor(3), &
                              tag, cart_neighbor(4), tag, mpi_comm_world, &
                              status, ierr)

    call mpi_isend(send_buf, npo*6, mpi_real, cart_neighbor(3), &
                   tag, mpi_comm_world, srequest, sierr)
    call mpi_irecv(recv_buf, npi*6, mpi_real, cart_neighbor(4), &
                   tag, mpi_comm_world, rrequest, rierr)
    call mpi_wait(srequest, sstatus, sierr)
    call mpi_wait(rrequest, rstatus, rierr)

    do pp = 1, npi
      xp_buf(:, np_buf+pp) = recv_buf((pp-1)*6+1:pp*6)
      xp_buf(2, np_buf+pp) = min(xp_buf(2, np_buf+pp)+ub, ub-eps)
    enddo

    pp = 1
    do
      if (pp > npi) exit
      x = xp_buf(:, np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local = np_local+1
        xvp(:, np_local) = x
        xp_buf(:, np_buf+pp) = xp_buf(:, np_buf+npi)
        npi=npi-1
        cycle
      endif
      pp = pp + 1
    enddo

    np_buf = np_buf + npi

    !
    ! Pass +z
    !

    tag = 15
    npo = 0
    pp  = 1
    do
      if (pp > np_buf) exit
      if (xp_buf(3, pp) >= ub) then
        npo = npo + 1
        send_buf((npo-1)*6+1:npo*6) = xp_buf(:, pp)
        xp_buf(:, pp) = xp_buf(:, np_buf)
        np_buf = np_buf - 1
        cycle
      endif
      pp = pp + 1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(2), &
                              tag,cart_neighbor(1),tag,mpi_comm_world, &
                              status,ierr)

    call mpi_isend(send_buf,npo*6,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*6,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*6+1:pp*6)
      xp_buf(3,np_buf+pp)=max(xp_buf(3,np_buf+pp)-ub,lb)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle
      endif
      pp=pp+1
    enddo

    np_buf=np_buf+npi

    !
    ! Pass -z
    !

    tag=16
    npo=0
    pp=1
    do
      if (pp > np_buf) exit
      if (xp_buf(3,pp) < lb) then
        npo=npo+1
        send_buf((npo-1)*6+1:npo*6)=xp_buf(:,pp)
        xp_buf(:,pp)=xp_buf(:,np_buf)
        np_buf=np_buf-1
        cycle
      endif
      pp=pp+1
    enddo

    npi = npo

    call mpi_sendrecv_replace(npi,1,mpi_integer,cart_neighbor(1), &
                              tag,cart_neighbor(2),tag,mpi_comm_world, &
                              status,ierr)

    call mpi_isend(send_buf,npo*6,mpi_real,cart_neighbor(1), &
                   tag,mpi_comm_world,srequest,sierr)
    call mpi_irecv(recv_buf,npi*6,mpi_real,cart_neighbor(2), &
                   tag,mpi_comm_world,rrequest,rierr)
    call mpi_wait(srequest,sstatus,sierr)
    call mpi_wait(rrequest,rstatus,rierr)

    do pp=1,npi
      xp_buf(:,np_buf+pp)=recv_buf((pp-1)*6+1:pp*6)
      xp_buf(3,np_buf+pp)=min(xp_buf(3,np_buf+pp)+ub,ub-eps)
    enddo

    pp=1
    do
      if (pp > npi) exit
      x=xp_buf(:,np_buf+pp)
      if (x(1) >= lb .and. x(1) < ub .and. x(2) >= lb .and. x(2) < ub .and. &
          x(3) >= lb .and. x(3) < ub ) then
        np_local=np_local+1
        xvp(:,np_local)=x
        xp_buf(:,np_buf+pp)=xp_buf(:,np_buf+npi)
        npi=npi-1
        cycle
      endif
      pp=pp+1
    enddo

    np_buf=np_buf+npi

  call mpi_reduce(np_buf,np_exit,1,mpi_integer,mpi_sum,0, mpi_comm_world,ierr)
  call mpi_reduce(int(np_local,kind=8),np_final,1,mpi_integer8,mpi_sum,0,  mpi_comm_world,ierr)

  if (rank == 0) then
    if (np_final /= int(np,kind=8)**3) then
        print *,'ERROR: total number of particles incorrect after passing'
    endif
  endif

  !! Check for particles out of bounds
  do i=1,np_local
    if (xvp(1,i) < 0 .or. xvp(1,i) >= nc_node_dim .or. &
         xvp(2,i) < 0 .or. xvp(2,i) >= nc_node_dim .or. &
         xvp(3,i) < 0 .or. xvp(3,i) >= nc_node_dim) then
       print *,'particle out of bounds',rank,i,xvp(:3,i),nc_node_dim
    endif
  enddo

    time2 = mpi_wtime(ierr)
    if (rank == 0) write(*,*) 'total buffered particles =', np_exit
    if (rank == 0) write(*,*) 'total particles          =', np_final
    if (rank == 0) write(*, *) "Finished pass_particles ... elapsed time = ", time2-time1

    return

end subroutine pass_particles

!!--------------------------------------------------------------!!

subroutine link_list
    !
    ! Construct a link list on the scale nm_node_dim which is coarsed by a
    ! factor of mesh_scale w.r.t. to nc_node_dim.
    !

    implicit none

    integer(4) :: i, j, k, pp
    real(8) :: sec1, sec2

    sec1 = mpi_wtime(ierr)

    ll(:) = 0
    do k = 1, nm_node_dim
        hoc(:, :, k) = 0
    enddo

    pp = 1
    do
        if (pp > np_local) exit
        i = floor(xvp(1, pp)/mesh_scale) + 1
        j = floor(xvp(2, pp)/mesh_scale) + 1
        k = floor(xvp(3, pp)/mesh_scale) + 1
        if (i < 1 .or. i > nm_node_dim .or. &
            j < 1 .or. j > nm_node_dim .or. &
            k < 1 .or. k > nm_node_dim) then
            write(*,*) "WARNING: LinkList Particle Deleted: ", xvp(1:3,pp)
            xvp(:,pp) = xvp(:,np_local)
            np_local = np_local - 1
            cycle
        else
            ll(pp) = hoc(i, j, k)
            hoc(i, j, k) = pp
        endif
        pp = pp + 1
    enddo

    sec2 = mpi_wtime(ierr)
    if (rank == 0) write(*,*) "Finished link_list ... elapsed time = ", sec2 - sec1

end subroutine link_list

!!--------------------------------------------------------------!!

subroutine zip_checkpoint
    !
    ! Outputs particle positions and velocities using the zip checkpoint
    ! technique.
    !

    implicit none

    character (len=200) :: fdm_zip0,fdm_zip1,fdm_zip2,fdm_zip3
    character (len=6) :: rank_s
    character (len=7) :: z_s
    integer(4), dimension(11) :: header_garbage
    real(4), parameter :: v_resolution = 32767.499
    integer(4) :: thread, ind
    real(4) :: vmax, vmax_local, v_r2i
    real(4) :: shake_offset(3)
    integer :: k, j, i, pp
    integer :: fstat

    real(8) :: sec1, sec2
    sec1 = mpi_wtime(ierr)

    vmax_local = maxval(abs(xvp(4:6,1:np_local)))
    call mpi_allreduce(vmax_local, vmax, 1, mpi_real, mpi_max, mpi_comm_world, ierr)
    v_r2i = v_resolution/vmax

    write(z_s,'(f7.3)') redshift
    z_s=adjustl(z_s)

    write(rank_s,'(i6)') rank
    rank_s=adjustl(rank_s)

#ifdef NEUTRINOS
    fdm_zip0=scratch_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip0_'//trim(rank_s)//'_nu.dat'
    fdm_zip1=scratch_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip1_'//trim(rank_s)//'_nu.dat'
    fdm_zip2=scratch_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip2_'//trim(rank_s)//'_nu.dat'
    fdm_zip3=scratch_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip3_'//trim(rank_s)//'_nu.dat'
#else
    fdm_zip0=scratch_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip0_'//trim(rank_s)//'.dat'
    fdm_zip1=scratch_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip1_'//trim(rank_s)//'.dat'
    fdm_zip2=scratch_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip2_'//trim(rank_s)//'.dat'
    fdm_zip3=scratch_path//'/node'//trim(rank_s)//'/'//trim(z_s)//'zip3_'//trim(rank_s)//'.dat'
#endif
    open(unit=10, file=fdm_zip0, status="replace", iostat=fstat, access="stream", buffered='yes')
    open(unit=11, file=fdm_zip1, status="replace", iostat=fstat, access="stream", buffered='yes')
    open(unit=12, file=fdm_zip2, status="replace", iostat=fstat, access="stream", buffered='yes')
    open(unit=13, file=fdm_zip3, status="replace", iostat=fstat, access="stream", buffered='yes')

    header_garbage(:) = 0
    shake_offset(:) = 0
    !write(10) np_local, header_garbage, v_r2i, shake_offset
    !write(11) np_local, header_garbage, v_r2i, shake_offset
    write(10) np_local,scalefactor,0.0,-3./sqrt(scalefactor),0,1000.0,1000.0,1000.0,1,1,1,1.,v_r2i,0.,0.,0.
    write(11) np_local,scalefactor,0.0,-3./sqrt(scalefactor),0,1000.0,1000.0,1000.0,1,1,1,1.,v_r2i,0.,0.,0.
    do k = 1, nm_node_dim
        do j = 1, nm_node_dim
            !$omp parallel default(shared) private(i, ind, pp)
            !$omp do schedule(static)
            do i = 1, nm_node_dim
                ind = 1
                rhoc_i4(i) = 0
                pp = hoc(i, j, k)
                do while (pp > 0)
                    rhoc_i4(i) = rhoc_i4(i) + 1
                    pos_i1(:, ind, i) = int(mod(xvp(1:3,pp)/mesh_scale,1.)*256,kind=1)
                    vel_i2(:, ind, i) = int(xvp(4:6,pp)*v_r2i,kind=2)
                    ind = ind + 1
                    pp = ll(pp)
                enddo
                if (ind > max_cell_np+1) then
                    write(*,*) "ERROR: max_cell_np too small !!", rank, ind, max_cell_np
                    call mpi_abort(mpi_comm_world, ierr, ierr)
                endif
            enddo
            !$omp end do
            !$omp end parallel
            write(12) int(min(rhoc_i4, 255), kind=1)
            write(13) pack(rhoc_i4, rhoc_i4>254)
            do i = 1, nm_node_dim
                write(10) pos_i1(:, 1:rhoc_i4(i), i)
                write(11) vel_i2(:, 1:rhoc_i4(i), i)
            enddo
        enddo
    enddo

    close(10)
    close(11)
    close(12)
    close(13)

    sec2 = mpi_wtime(ierr)
    if (rank == 0) then
        write(*,*) "vmax, v_r2i = ", vmax, v_r2i
        write(*,*) "Finished zip_checkpoint ... elapsed time = ", sec2 - sec1
    endif

end subroutine zip_checkpoint
#endif
!!--------------------------------------------------------------!!

  subroutine powerspectrum(pk)
    implicit none
    real, dimension(2,nc)       :: pk

    integer :: i,j,k,kg,mg,jg,ig
    integer :: k1,k2
    real    :: kr,kx,ky,kz,w1,w2,pow
#ifdef SLAB
    real, dimension(3,nc,nc_slab) :: pkt
#else
    real, dimension(3,nc,nc_pen+2) :: pkt
#endif
    real, dimension(3,nc) :: pktsum
    integer :: ind, dx, dxy

    real time1,time2
    call cpu_time(time1)

#ifndef SLAB
    dx  = fsize(1)
    dxy = dx * fsize(2)
    ind = 0
#endif

    pkt=0.0
    pktsum=0.0

    !! Compute power spectrum
#ifdef SLAB
    !$omp parallel do default(shared) private(i,j,k) &
    !$omp& private(kr,kx,ky,kz,kg,k1,k2,w1,w2,pow)
    do k=1,nc_slab
       kg=k+nc_slab*rank
       if (kg .lt. hc+2) then
          kz=kg-1
       else
          kz=kg-1-nc
       endif
       do j=1,nc
          if (j .lt. hc+2) then
             ky=j-1
          else
             ky=j-1-nc
          endif
          do i=1,nc+2,2
             kx=(i-1)/2
             kr=sqrt(kx**2+ky**2+kz**2)
             if (kr .ne. 0) then
                k1=ceiling(kr)
                k2=k1+1
                w1=k1-kr
                w2=1-w1
                pow=sum((slab(i:i+1,j,k)/ncr**3)**2)
                pkt(1,k1,k)=pkt(1,k1,k)+w1*pow
                pkt(2,k1,k)=pkt(2,k1,k)+w1*pow**2
                pkt(3,k1,k)=pkt(3,k1,k)+w1
                pkt(1,k2,k)=pkt(1,k2,k)+w2*pow
                pkt(2,k2,k)=pkt(2,k2,k)+w2*pow**2
                pkt(3,k2,k)=pkt(3,k2,k)+w2
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do
#else
    !! Cannot thread because of ind index
    do k = 1, nc_pen+mypadd
        do j = 1, nc_node_dim
            do i = 1, nc, 2
                kg = ind / dxy
                mg = ind - kg * dxy
                jg = mg / dx
                ig = mg - jg * dx
                kg = fstart(3) + kg
                jg = fstart(2) + jg
                ig = 2 * (fstart(1) + ig) - 1
                ind = ind + 1
                if (kg < hc+2) then
                    kz=kg-1
                else
                    kz=kg-1-nc
                endif
                if (jg < hc+2) then
                    ky=jg-1
                else
                    ky=jg-1-nc
                endif
                kx = (ig-1)/2
                kr=sqrt(real(kx**2+ky**2+kz**2))
                if (kr .ne. 0) then
                    k1=ceiling(kr)
                    k2=k1+1
                    w1=k1-kr
                    w2=1-w1
                    pow=sum((slab(i:i+1,j,k)/ncr**3)**2)
                    pkt(1,k1,k)=pkt(1,k1,k)+w1*pow
                    pkt(2,k1,k)=pkt(2,k1,k)+w1*pow**2
                    pkt(3,k1,k)=pkt(3,k1,k)+w1
                    pkt(1,k2,k)=pkt(1,k2,k)+w2*pow
                    pkt(2,k2,k)=pkt(2,k2,k)+w2*pow**2
                    pkt(3,k2,k)=pkt(3,k2,k)+w2
                endif
            enddo
        enddo
    enddo
#endif

    !! Merge power spectrum from threads
#ifdef SLAB
    do k=2,nc_slab
#else
    do k=2,nc_pen+2
#endif
       pkt(:,:,1)=pkt(:,:,1)+pkt(:,:,k)
    enddo

    !! Reduce to rank 0
    call mpi_reduce(pkt(:,:,1),pktsum,3*nc,mpi_real,mpi_sum,0,mpi_comm_world,ierr)

    !! Divide by weights
    !! pk(1,k) stores pk(k)
    !! pk(2,k) stores standard deviation
    if (rank == 0) then
      do k=1,nc
        if (pktsum(3,k) .eq. 0) then
          pk(:,k)=0
        else
          pk(:,k)=pktsum(1:2,k)/pktsum(3,k)
          pk(2,k)=sqrt(abs((pk(2,k)-pk(1,k)**2)/(pktsum(3,k)-1)))
          pk(1:2,k)=4.*pi*(real(k)-1.)**3*pk(1:2,k)
       endif
      enddo
    endif

    call mpi_bcast(pk,2*nc,mpi_real,0,mpi_comm_world,ierr)

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called power spectrum'
    return
  end subroutine powerspectrum

!!------------------------------------------------------------------!!

  subroutine initvar
    implicit none
    integer :: k

    real time1,time2
    call cpu_time(time1)

#ifdef DEBUG_LOW
    do k=0,nodes-1
      if (rank==k) print *,'rank:',rank,'starting initvar'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    !! Initialize large arrays in parallel
    !$omp parallel default(shared) private(k)
    !$omp do
    do k=0,nc_node_dim+1
       phi(:,:,k)=0
    enddo
    !$omp end do
    !$omp do
#ifdef SLAB
    do k=1,nc_slab
#else
    do k=1,nc_pen+2
#endif
       slab(:,:,k)=0
#ifdef ZIP
       slab_cache(:,:,k)=0
#endif
    enddo
    !$omp end do
    !$omp do
    do k=1,nk
       tf(:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,nc
       pkm(:,k)=0
    enddo
    !$omp end do
    !$omp do
    do k=1,nc
       pkn(:,k)=0
    enddo
    !$omp end do
    !$omp end parallel

    !! Initialize fftw so that it generates the plans!
    firstfftw=.true.

#ifdef DEBUG_LOW
    do k=0,nodes-1
      if (k==rank) print *,rank,'finished initvar'
      call mpi_barrier(mpi_comm_world,ierr)
    enddo
#endif

    call cpu_time(time2)
    time2=(time2-time1)
    if (rank == 0) write(*,"(f8.2,a)") time2,'  Called init var'
    return
  end subroutine initvar

!!------------------------------------------------------------------!!

subroutine mesh_buffer
  implicit none

  integer(4) :: buffer_size
  integer(4) :: tag
  integer(4) :: status(MPI_STATUS_SIZE)

#ifdef DEBUG_BUFFER_MESH
    write(*,*) 'calling mesh_buffer on rank ', rank
#endif

  buffer_size = (nc_node_dim + 2)**2

  if (rank==0) write(*,*) 'buffer_size =', buffer_size


  tag=64

!! send to node in -x

      phi_buf(:,:)=phi(1,:,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(5),tag,cart_neighbor(6), &
                              tag,mpi_comm_cart,status,ierr)

      phi(nc_node_dim+1,:,:)=phi(nc_node_dim+1,:,:)+phi_buf(:,:)

#ifdef DEBUG_BUFFER_MESH
    write(*,*) 'finished -x exchanged on node ', rank
#endif


!! send to node in +x

      phi_buf(:,:)=phi(nc_node_dim,:,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(6),tag,cart_neighbor(5), &
                              tag,mpi_comm_cart,status,ierr)

      phi(0,:,:)=phi(0,:,:)+phi_buf(:,:)

!! send to node in -y

      phi_buf(:,:)=phi(:,1,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(3),tag,cart_neighbor(4), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,nc_node_dim+1,:)=phi(:,nc_node_dim+1,:)+phi_buf(:,:)

!! send to node in +y

      phi_buf(:,:)=phi(:,nc_node_dim,:)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(4),tag,cart_neighbor(3), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,0,:)=phi(:,0,:)+phi_buf(:,:)

!! send to node in -z

      phi_buf(:,:)=phi(:,:,1)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(1),tag,cart_neighbor(2), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,:,nc_node_dim+1)=phi(:,:,nc_node_dim+1)+phi_buf(:,:)

!! send to node in +z

      phi_buf(:,:)=phi(:,:,nc_node_dim)

      call mpi_sendrecv_replace(phi_buf,buffer_size,mpi_real, &
                              cart_neighbor(2),tag,cart_neighbor(1), &
                              tag,mpi_comm_cart,status,ierr)

      phi(:,:,0)=phi(:,:,0)+phi_buf(:,:)

  end subroutine mesh_buffer

!!------------------------------------------------------------------!!

  function power(kr,ix,iy)
   implicit none
    real    :: kr
    integer :: ix,iy

    integer :: i,i1,i2
    real    :: x,y,x1,x2,y1,y2
    real    :: power

    i1=1
    i2=nk
    do while (i2-i1 .gt. 1)
       i=(i1+i2)/2
       if (kr .gt. tf(ix,i)) then
          i1=i
       else
          i2=i
       endif
    enddo

    x1=log(tf(ix,i1))
    y1=log(tf(iy,i1))
    x2=log(tf(ix,i2))
    y2=log(tf(iy,i2))
    x =log(kr)
    y =y1+(y2-y1)*(x-x1)/(x2-x1)
    power=exp(y)

    return
  end function power

!!------------------------------------------------------------------!!

  function Dgrow(a)
    implicit none
    real, parameter :: om=omegam
    real, parameter :: ol=omegal
    real a
    real Dgrow

    real g,ga,hsq,oma,ola

    hsq=om/a**3+(1-om-ol)/a**2+ol
    oma=om/(a**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))

    Dgrow=a*ga/g
    return
  end function Dgrow

!!------------------------------------------------------------------!!

  function vfactor(a)
    implicit none
    real :: a

    real :: H,km,lm
    real :: vfactor

    lm=omegal/omegam
    km=(1-omegam-omegal)/omegam
    H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
    vfactor=a**2*H

    return
  end function vfactor

!!------------------------------------------------------------------!!

  function tophat(x)
    implicit none
    real :: x,tophat

    if (x .ne. 0) then
       tophat=3*(sin(x)-cos(x)*x)/x**3
    else
       tophat=1
    endif

    return
  end function tophat

!!------------------------------------------------------------------!!

#ifdef WDM
function wdm_transfer(x)
    !
    ! Returns WDM transfer function using equations (5) and (6) of Schneider, Smith, and Reed (2013)
    !

    implicit none

    real :: x, wdm_transfer

    wdm_transfer = (1. + (alpha*x)**(2.*mu))**(-5./mu)

    return

end function wdm_transfer
#endif

!!------------------------------------------------------------------!!

#ifdef BROKEN_POWER
function broken_power_transfer(x)
    !
    ! Applies a broken power law at a small scale pivot (k0) as an alternative to WDM.
    !

    implicit none

    real :: x, broken_power_transfer

    if (x < k0) then
        broken_power_transfer = 1.
    else
        broken_power_transfer = (x/k0)**(dns/2.)
    endif

    return

end function broken_power_transfer
#endif

end program dist_init

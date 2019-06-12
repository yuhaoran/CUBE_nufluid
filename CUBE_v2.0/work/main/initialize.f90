subroutine initialize
  use variables
  use cubefft
  use pencil_fft
  implicit none
  save
  include 'fftw3.f'

  logical exist_neu_checkpoint

  !call system('hostname')
  if (this_image()==1) then
    print*, ''
    print*, 'CUBE run on',int(nn**3,2),'images  x',int(ncore,1),'cores'
    print*, ''
    print*, 'initialize'
    print*, '  call geometry'
  endif
  sync all

  call geometry

  if (head) print*, '  call create_cubefft_plan ng = ',nfe
  call system_clock(t1,t_rate)
  call create_cubefft_plan
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
  sync all

  if (head) print*, '  call create_penfft_plan nfe = ',ng
  call system_clock(t1,t_rate)
  call create_penfft_plan
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
  sync all

  call kernel_f
  call kernel_c

  dt_neu=dt_max
  neutrino_flag=sum(f_neu).gt.0

  t=0
  dt=0
  dt_old=0
  da=0
  cur_checkpoint=1 ! change for resuming checkpoints
  cur_halofind=1
  z_checkpoint=-0.5
  z_checkpoint=-0.5
  checkpoint_step=.false.
  halofind_step=.false.
  final_step=.false.

    open(16,file='z_checkpoint.txt',status='old')
    do i=1,nmax_redshift-1
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
    if (n_checkpoint==0) stop 'z_checkpoint.txt empty'
# ifdef HALOFIND
    open(16,file='z_halofind.txt',status='old')
    do i=1,nmax_redshift-1
      read(16,end=81,fmt='(f8.4)') z_halofind(i)
    enddo
    81 n_halofind=i-1
    close(16)
    if (n_halofind==0) stop 'z_halofind.txt empty'
    n_halofind=n_halofind[1]
    z_halofind(:)=z_halofind(:)[1]
# endif
  tau=-3/sqrt(1./(1+z_checkpoint(cur_checkpoint)))
  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]
  ! create output directories
  call system('mkdir -p '//opath//'/image'//image2str(image))

# ifdef HALOFIND
  if (head) then
    print*, ''
    print*, 'runtime halofind information'
    print*, '  ',z_checkpoint(1),'< CDM initial conditions'
    do i=1,n_halofind
      print*, '  ',z_halofind(i)
    enddo
  endif
  sync all
# endif

  n_checkpoint_neu=0
# ifdef NEUTRINOS
    if (z_i_nu>z_i) then
      if (head) then
        print*, '  error: z_i_nu>z_i'
        print*, '  z_i_nu, z_i =',z_i_nu,z_i
      endif
      stop
    elseif (z_i_nu==z_i) then
      neutrino_flag=.true.
    else !(neutrinos will evolve at a later redshift)
      neutrino_flag=.false.
      exist_neu_checkpoint=.false.
      do i=1,n_checkpoint
        exist_neu_checkpoint=exist_neu_checkpoint.or.(z_checkpoint(i)==z_i_nu)
        if (exist_neu_checkpoint) then
          n_checkpoint_neu=i
          exit
        endif
      enddo
      if (.not. exist_neu_checkpoint) then
        if (z_checkpoint(n_checkpoint)>z_i_nu) then
          if (head) then
            print*, '  warning: neutrinos will not be added till the end'
            print*, '  final checkpoint:',z_checkpoint(n_checkpoint)
            print*, '  z_i_nu:',z_i_nu
            print*, '  will add one more checkpoint.'
          endif
          z_checkpoint(n_checkpoint+1)=z_i_nu
          n_checkpoint=n_checkpoint+1
          n_checkpoint_neu=n_checkpoint
        else
          do i=1,n_checkpoint
            if (z_checkpoint(i)<z_i_nu) then
              print*,'z_checkpoint(i),z_i_nu=',z_checkpoint(i),z_i_nu
              z_checkpoint(i+1:n_checkpoint+1)=z_checkpoint(i:n_checkpoint)
              z_checkpoint(i)=z_i_nu
              n_checkpoint=n_checkpoint+1
              n_checkpoint_neu=i
              exit
            endif
          enddo
        endif
      endif
    endif
# else
    neutrino_flag=.false.
# endif
  if (head) then
    print*, ''
    print*, 'checkpoint information'
    print*, '  ',z_checkpoint(1),'< CDM initial conditions'
    do i=2,n_checkpoint
      print*, '  ',z_checkpoint(i),merge('< adding neutrinos','                  ',i==n_checkpoint_neu)
    enddo
  endif
  sync all

  print*,'OpenMP information'
  nth=omp_get_num_procs()
  print*,'  omp_get_num_procs =',nth
  nth=omp_get_thread_limit()
  print*,'  omp_get_thread_limit =',nth
  call omp_set_num_threads(ncore)
  nth=omp_get_max_threads()
  print*,'  omp_get_max_threads =',nth
  print*, ''
endsubroutine

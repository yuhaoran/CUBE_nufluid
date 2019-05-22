!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@cita.utoronto.ca   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#define record_Fc
program main
  use omp_lib
  use variables
  use buffer_grid_subroutines
  use buffer_particle_subroutines
  use update_particle
  use pp_force
  implicit none
  save

  call initialize
  call particle_initialization
  call buffer_grid
  call buffer_x
  call buffer_v
  cur_checkpoint=cur_checkpoint+1
  cur_halofind=cur_checkpoint+1
  if (head) open(77,file=output_dir()//'vinfo'//output_suffix(),access='stream',status='replace')

  if (head) print*, '---------- starting main loop ----------'
  DO istep=1,istep_max
    call system_clock(ttt1,t_rate)
    call timestep
    call update_x
    call buffer_grid
    call buffer_x
    if (Extended_pp_force) then
      call ext_pp_force
    endif
    call particle_mesh
    call buffer_v
    if (checkpoint_step .or. halofind_step) then
      dt_old=0
      call update_x
      if (checkpoint_step) then
        call checkpoint
        cur_checkpoint=cur_checkpoint+1
      endif
      call buffer_grid
      call buffer_x
      call buffer_v
      if (halofind_step) then
        call halofind
        cur_halofind=cur_halofind+1
      endif
      call print_header(sim)
      if (final_step) exit
      dt=0
    endif
    call system_clock(ttt2,t_rate)
    print*, 'total elapsed time =',real(ttt2-ttt1)/t_rate,'secs';
  ENDDO
  if (head) close(77)
  call finalize

contains

  function Dgrow(scale_factor)
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real scale_factor
    real Dgrow
    real g,ga,hsq,oma,ola
    hsq=om/scale_factor**3+(1-om-ol)/scale_factor**2+ol
    oma=om/(scale_factor**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
    Dgrow=scale_factor*ga/g
  end function Dgrow

endprogram

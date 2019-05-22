subroutine finalize
  use variables
  use cubefft
  use pencil_fft
  implicit none
  save
  !include 'fftw3.f'

  call destroy_cubefft_plan

  call destroy_penfft_plan

endsubroutine

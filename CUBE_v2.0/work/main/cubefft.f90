#define macbook
!  if defined macbook, then omp-fftw is disabled.
module cubefft
  contains
  subroutine create_cubefft_plan
    use,intrinsic :: ISO_C_BINDING
    use variables
    !use intrinsic :: ISO_C_BINDING
    implicit none
    save
    !include 'fftw3.f'
    include 'fftw3.f03'
    integer istat,icore

#ifndef macbook
    call sfftw_init_threads(istat)
    if(head) print*, 'sfftw_init_threads status',istat
    icore=omp_get_max_threads()
    if(head) print*, 'omp_get_max_threads() =',icore
    call sfftw_plan_with_nthreads(icore)

    !call sfftw_plan_with_nthreads(32)
#endif

    call sfftw_plan_dft_r2c_3d(plan_fft_fine,nfe,nfe,nfe,rho_f,rho_f,FFTW_MEASURE)
    call sfftw_plan_dft_c2r_3d(plan_ifft_fine,nfe,nfe,nfe,rho_f,rho_f,FFTW_MEASURE)
  endsubroutine create_cubefft_plan

  subroutine destroy_cubefft_plan
    use,intrinsic :: ISO_C_BINDING
    use variables
    implicit none
    save
    !include 'fftw3.f'
    include 'fftw3.f03'
    call sfftw_destroy_plan(plan_fft_fine)
    call sfftw_destroy_plan(plan_ifft_fine)
    !call fftw_cleanup_threads()
  endsubroutine destroy_cubefft_plan

endmodule

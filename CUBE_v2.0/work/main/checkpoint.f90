subroutine checkpoint
  use omp_lib
  use variables
  use neutrinos
  implicit none
  save

  character(100) fn10,fn11,fn12,fn13,fn14,fn15,fn21,fn22,fn23,fn24,fn25

  if (head) print*, 'checkpoint'

  fn10=output_name('info')
  fn11=output_name('xp')
  fn12=output_name('vp')
  fn13=output_name('np')
  fn14=output_name('vc')
  fn15=output_name('id')
  fn21=output_name('xp_nu')
  fn22=output_name('vp_nu')
  fn23=output_name('np_nu')
  fn24=output_name('vc_nu')
  fn25=output_name('id_nu')

  sim%a=a
  sim%t=t
  sim%tau=tau

  sim%cur_checkpoint=cur_checkpoint

  sim%vsim2phys=(1.5/a)*box*100.*sqrt(omega_m)/nf_global
  sim%sigma_vi=sigma_vi
  sim%sigma_vi_nu=sigma_vi_nu

  !$omp parallelsections default(shared)
  !$omp section
  open(10,file=fn10,status='replace',access='stream'); write(10) sim; close(10)
  !$omp section
  open(11,file=fn11,status='replace',access='stream'); write(11) xp(:,:sim%nplocal); close(11)
  !$omp section
  open(12,file=fn12,status='replace',access='stream'); write(12) vp(:,:sim%nplocal); close(12)
  !$omp section
  open(13,file=fn13,status='replace',access='stream'); write(13) rhoc(1:nt,1:nt,1:nt,:,:,:); close(13)
  !$omp section
  open(14,file=fn14,status='replace',access='stream'); write(14) vfield(:,1:nt,1:nt,1:nt,:,:,:); close(14)
# ifdef PID
    !$omp section
    open(15,file=fn15,status='replace',access='stream'); write(15) pid(:sim%nplocal); close(15)
    print*, 'check PID range: ',minval(pid(:sim%nplocal)),maxval(pid(:sim%nplocal))
    print*, 'check PID range:',minval(pid(:sim%nplocal)),maxval(pid(:sim%nplocal))
    if (minval(pid(:sim%nplocal))<1) then
      print*, 'pid are not all positive'
      
    endif
# endif
  !$omp endparallelsections
  print*,'  image',this_image(),'wrote',sim%nplocal,'CDM particles'

# ifdef NEUTRINOS
  if (neutrino_flag) then
    !$omp parallelsections default(shared)
    !$omp section
    open(21,file=fn21,status='replace',access='stream'); write(21) xp_nu(:,:sim%nplocal_nu); close(21)
    !$omp section
    open(22,file=fn22,status='replace',access='stream'); write(22) vp_nu(:,:sim%nplocal_nu); close(22)
    !$omp section
    open(23,file=fn23,status='replace',access='stream'); write(23) rhoc_nu(1:nt,1:nt,1:nt,:,:,:); close(23)
    !$omp section
    open(24,file=fn24,status='replace',access='stream'); write(24) vfield_nu(:,1:nt,1:nt,1:nt,:,:,:); close(24)
#   ifdef EID
      !$omp section
      open(25,file=fn25,status='replace',access='stream'); write(25) pid_nu(:sim%nplocal_nu); close(25)
#   endif
    !$omp endparallelsections
  endif
# endif

  if (neutrino_flag) then
    print*,'  image',this_image(),'wrote',sim%nplocal_nu,'neutrino particles'
  endif

  sync all



  if (cur_checkpoint==n_checkpoint_neu) then
    print*, cur_checkpoint, n_checkpoint_neu
    print*, 'turn neutrino_flag on'
    neutrino_flag=.true.
    !sim%mass_p_cdm=real((nf*nn)**3*f_cdm)/sim%npglobal
    !sim%mass_p_nu=real((nf*nn)**3*f_nu)/sim%npglobal_nu
  endif
  !npglobal=0
  !do i=1,nn**3
  !  npglobal=npglobal+sim%nplocal[i]
  !enddo
  sync all


endsubroutine

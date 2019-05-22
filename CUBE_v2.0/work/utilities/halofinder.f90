!#define write_snapshot
program halofinder
  use omp_lib
  use parameters
  use buffer_grid_subroutines
  use buffer_particle_subroutines
#ifdef write_snapshot
  use variables
#else
  use variables, only : xp,vp,rhoc,pid,vfield,i
#endif
  implicit none

#ifdef write_snapshot
  !integer itx,ity,itz,j,k,l,np
  !integer(8) ip,nzero
  !real sigma_vi,xq(3),vreal
  !sigma_vi=sim%sigma_vi
#endif

  call geometry
  call omp_set_num_threads(ncore)
  if (head) then
    print*, 'halofinder at:'
    open(16,file='../main/z_checkpoint.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_halofind(i)
      print*, z_halofind(i)
    enddo
    71 n_halofind=i-1
    close(16)
    print*,''
  endif

  sync all
  n_halofind=n_halofind[1]
  z_halofind(:)=z_halofind(:)[1]
  n_checkpoint=n_halofind
  z_checkpoint=z_halofind
  sync all

  do cur_halofind=n_checkpoint,n_checkpoint
    cur_checkpoint=cur_halofind
    xp=0; vp=0; rhoc=0; vfield=0
    call particle_initialization
#   ifdef write_snapshot
      call spine_image(rhoc,idx_b_l,idx_b_r,ppl0,pplr,pprl,ppr0,ppl,ppr)
      !print*,2.7755e11*50**3*0.32/288**3/0.67 ! TianNu / TianZero CDM particle mass in unit of solar mass
      print*, 'Writing traditional snapshot file'
      open(11,file=output_name_halo('snapshot'),status='replace',access='stream')
      write(11) 256
      write(11) 0,int(sim%nplocal,4),0,0,0,0
      print*,'  Npart=',int(sim%nplocal,4)
      write(11) 0d0,2.7755d11*box**3*omega_m/sim%nplocal/1e10,0d0,0d0,0d0,0d0
      print*,'  Massarr[1]=',2.7755d11*box**3*omega_m/sim%nplocal/1e10
      write(11) 1d0/(z_halofind(cur_halofind)+1) ! Time
      print*,'  Time=',1d0/(z_halofind(cur_halofind)+1)
      write(11) 1d0*z_halofind(cur_halofind) ! Redshift
      print*,'  Redshift=',1d0*z_halofind(cur_halofind)
      write(11) 0,0 ! FlagSfr, FlagFeedback
      write(11) 0,int(sim%npglobal,4),0,0,0,0
      print*,'  Nall=',int(sim%npglobal,4)
      write(11) 0 ! FlagCooling
      write(11) 1 ! NumFiles
      write(11) box*1000d0 ! BoxSize in kpc/h
      print*,'  BoxSize=',box*1d3
      write(11) 1d0*omega_m,1d0*omega_l,1d0*h0
      print*,'  omega_m,omega_l,h0=',1d0*omega_m,1d0*omega_l,1d0*h0
      write(11) 0,0
      write(11) 0,0,0,0,0,0
      write(11) 0 ! flag_entr_ics ! 49 single pricision numbers
      write(11) spread(0,1,15) ! another 15 empty entries
      write(11) 256

      write(11) int(sim%nplocal,4)*12
      ! write pos
      do itz=1,nnt ! loop over tile
      do ity=1,nnt
      do itx=1,nnt
        do k=1,nt
        do j=1,nt
        do i=1,nt
          np=rhoc(i,j,k,itx,ity,itz)
          nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
          do l=1,np
            ip=nzero+l
            xq=((/i,j,k/)-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
            write(11) real(xq/nc*box*1d3,4)
            !vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
            !vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
            !write(11) real(vreal*sim%vsim2phys,4)
          enddo
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      write(11) int(sim%nplocal,4)*12

      write(11) int(sim%nplocal,4)*12
      ! write vel
      do itz=1,nnt ! loop over tile
      do ity=1,nnt
      do itx=1,nnt
        do k=1,nt
        do j=1,nt
        do i=1,nt
          np=rhoc(i,j,k,itx,ity,itz)
          nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
          do l=1,np
            ip=nzero+l
            !xq=((/i,j,k/)-1d0) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
            !write(11) real(xq/nc*box*1d3,4)
            vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sigma_vi*vrel_boost))
            vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
            write(11) real(vreal*sim%vsim2phys,4)
          enddo
        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      write(11) int(sim%nplocal,4)*12

      ! write id
      write(11) int(sim%nplocal,4)*4
      write(11) pid(:sim%nplocal)
      write(11) int(sim%nplocal,4)*4
      close(11)
#   endif
    call buffer_grid
    call buffer_x
    call buffer_v
    call halofind
  enddo
  sync all
  if (head) print*, 'halofinder done'

end

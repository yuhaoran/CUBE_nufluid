!#define eig
!#define Voronoi
!#define NEUTRINOS

module halo_output
  implicit none

  type type_halo_catalog
    integer nhalo_tot,nhalo
    real den_odc
  endtype

  type type_halo_info
    real hpos(3)
    real mass_odc,radius_odc,v_disp
    real x_mean(3),v_mean(3),ang_mom(3),var_x(3),inertia(3,3)
    real q_mean(3),inertia_q(3,3),s_mean(3)
  endtype
endmodule

program ang_mom_corr
  use parameters
  use halo_output, only: type_halo_catalog, type_halo_info
  implicit none
  save
  ! nc: coarse grid per node per dim
  ! nf: fine grid per node per dim



  integer,parameter :: nc_halo_max=128
  integer,parameter :: newbox_max=30
  integer,parameter :: nlist=5*(nc_halo_max+1)**3
  integer,parameter :: n_peak_max=5*nc_halo_max**3
  integer,parameter :: max_halo_np=5*(nc_halo_max+1)**3/3
  integer(8) i,j,k,l,ntemp
  integer(4) i0,ihalo,nhalo,n_compute,ip,np,nptemp,ncheck
  character(20) str_z,str_i
  integer pid(max_halo_np),ipos(3),ieig
  real qpos(3,max_halo_np),qpos_mean(3),vq(3),dx(3),vf,scale_factor,dx_mean(3)
  real spin_q(3),spin_t(3),spin_x(3),spin_u(3),spin_v(3),spin_r(3),spin_e(3,3),spin_vor(3)
  real inertia(3,3),tide(3,3),tidu(3,3),torque(3,3),torquu(3,3),torqur(3,3),torque_c,torque_u
  real force_cdm(3),force_neu(3),torque_v(3,3)
  real phi(0:nf+1,0:nf+1,0:nf+1)[*]
  real phu(0:nf+1,0:nf+1,0:nf+1)[*]

  type(type_halo_info) halo_info
  type(type_halo_catalog) halo_catalog
  integer(4),allocatable :: isort_mass(:)
  real,allocatable :: inert_v(:,:,:),inert_temp(:,:,:),Tc3(:,:,:,:),halonp(:),theta(:,:) !!qx,qt,tx,qu,uv,vr,ur,ff,ve1,ve2,ve3

  call geometry
  if (head) then
    print*, 'angular momemtum correlation on resolution:'
    print*, 'ng=',ng
    print*, 'ng*nn=',ng*nn
  endif
  sync all

  if (head) then
    print*, 'use halo catalog at:'
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

  cur_checkpoint=1 ! read phi for the initial condition
  open(21,file=output_name('phi1'),status='old',access='stream')
  !open(21,file='../../output/universe16/image1/50.000_phi1_1.bin',access='stream')
  !read(21) phi(1:nf,1:nf,1:nf)
  !open(21,file='../../output/universe2/image1/5.000_phi_z5_1.bin',status='old',access='stream')
  !open(21,file='../../output/universe3/image1/0.000_phiE_1.bin',status='old',access='stream')
  read(21) phi(1:nf,1:nf,1:nf)
  close(21)

  write(str_i,'(i6)') image
  !print*, output_name('phi1_nu')
  !open(21,file=output_name('phi1_nu'),status='old',access='stream')
#ifdef NEUTRINOS
  open(21,file=opath//'image'//trim(adjustl(str_i))//'/5.000_tf5.000_phi1_nu'//output_suffix(),status='old',access='stream')
    !open(21,file='../../output/universe16/image1/50.000_phi1_1.bin',status='old',access='stream')
  read(21) phu(1:nf,1:nf,1:nf)
  close(21)
#endif
  vf=vfactor(1/(1+z_checkpoint(cur_checkpoint)))

  if (head) print*, '  buffer phi'
  phi(0,:,:)=phi(nf,:,:)[image1d(inx,icy,icz)]
  phi(nf+1,:,:)=phi(1,:,:)[image1d(ipx,icy,icz)]
  sync all
  phi(:,0,:)=phi(:,nf,:)[image1d(icx,iny,icz)]
  phi(:,nf+1,:)=phi(:,1,:)[image1d(icx,ipy,icz)]
  sync all
  phi(:,:,0)=phi(:,:,nf)[image1d(icx,icy,inz)]
  phi(:,:,nf+1)=phi(:,:,1)[image1d(icx,icy,ipz)]
  sync all
  phu(0,:,:)=phu(nf,:,:)[image1d(inx,icy,icz)]
  phu(nf+1,:,:)=phu(1,:,:)[image1d(ipx,icy,icz)]
  sync all
  phu(:,0,:)=phu(:,nf,:)[image1d(icx,iny,icz)]
  phu(:,nf+1,:)=phu(:,1,:)[image1d(icx,ipy,icz)]
  sync all
  phu(:,:,0)=phu(:,:,nf)[image1d(icx,icy,inz)]
  phu(:,:,nf+1)=phu(:,:,1)[image1d(icx,icy,ipz)]
  sync all





  do cur_checkpoint= n_checkpoint,n_checkpoint
    if (head) print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    scale_factor=1/(1+z_checkpoint(cur_checkpoint))
    ! read halo catalog and their pid
    open(11,file=output_name('halo'),status='old',access='stream')
    read(11) halo_catalog
    print*,' N_halos_global =',halo_catalog%nhalo_tot
    print*,' N_halos_local  =',halo_catalog%nhalo
    print*,' Overdensity    =',halo_catalog%den_odc
    nhalo=halo_catalog%nhalo

    open(12,file=output_name('halo_pid'),status='old',access='stream') ! for computing protohalo regions
    ! read halo_init_spin file for writing
    open(13,file=output_name('halo_init_spin'),status='replace',access='stream') ! to write
    open(14,file=output_name('halo_init_Tc'),status='replace',access='stream') ! to write
    open(16,file=output_name('halo_init_I'),status='replace',access='stream') ! to write
#ifdef Voronoi
    open(19,file=output_name('protohalo_vor_I'),status='old',access='stream')
    allocate(inert_v(3,3,nhalo),inert_temp(3,3,nhalo))
    read(19) inert_temp
    close(19)
    !print*,'read ntemp'; read(*,*) ntemp
    ntemp=0
    ! randomize
    do ihalo=1,nhalo
      inert_v(:,:,ihalo)=inert_temp(:,:,modulo(ihalo-1+ntemp,nhalo)+1)
    enddo
    !open(19,file='0.000_halo_vor_I_1_res150.bin',status='old',access='stream')
#endif
    write(14) nhalo ! write number of halos as header
    write(16) nhalo
#   ifdef eig
      !print*, 'read',output_name('halo_init_Tc3')
      !open(15,file=output_name('halo_init_Tc3'),status='old',access='stream')
      print*, 'read',output_name('protohalo_vor_I3')
      open(15,file=output_name('protohalo_vor_I3'),status='old',access='stream')
      allocate(Tc3(3,3,3,nhalo))
      read(15) Tc3
      close(15)
#   endif

    allocate(theta(0:12,nhalo),halonp(nhalo),isort_mass(nhalo))
    !! calculate theta: qx,qt,tx,qu,uv,vr,ur,ff,ve1,ve2,ve3,voronoi
    do ihalo=1,nhalo
      pid=0; np=0
      read(11) halo_info
      read(12) np ! PID
      read(12) pid(1:np)
      read(12) ncheck
      halonp(ihalo)=real(np) ! for sorting output
      if (ncheck/=0) stop "halo_pid file error"
      if (minval(pid(1:np))<0) then
        print*, "  warning: pids are not all positive numbers"
        print*,ihalo,np
        print*,minloc(pid(1:np))
        print*,minval(pid(1:np))
        print*,'negative pid:',pid(minloc(pid(1:np)))
        print*,pid(1:np)
        stop
      endif
      pid(:np)=pid(:np)-1
      dx_mean=0; nptemp=0
      do ip=1,np
        if (pid(ip)==-1) then
          cycle
        endif
        nptemp=nptemp+1
        qpos(3,ip)=pid(ip)/nf_global**2
        qpos(2,ip)=(pid(ip)-(pid(ip)/nf_global**2)*nf_global**2)/nf_global
        qpos(1,ip)=modulo(pid(ip),nf_global)
        qpos(:,ip)=qpos(:,ip)+0.5
        dx=qpos(:,ip)-halo_info%x_mean
        dx=modulo(dx+nf_global/2,real(nf_global))-nf_global/2
        dx_mean=dx_mean+dx
      enddo
      dx_mean=dx_mean/nptemp
      qpos_mean=halo_info%x_mean+dx_mean
      qpos_mean=modulo(qpos_mean,real(nf_global))
      if(ihalo==1) then
        print*,'check halo position 1'
        print*,halo_info%x_mean
        print*,dx_mean
        print*,halo_info%x_mean+dx_mean
        print*,qpos_mean
      endif
      spin_q=0; spin_t=0; spin_x=0; spin_u=0; spin_v=0; spin_r=0;
      force_cdm=0; force_neu=0
      torque_c=0; torque_u=0
      dx_mean=0
      inertia=0
      tide=0
      tidu=0

      do ip=1,np
        if (pid(ip)==-1) then
          cycle
        endif
        ipos(3)=pid(ip)/nf_global**2
        ipos(2)=(pid(ip)-ipos(3)*nf_global**2)/nf_global
        ipos(1)=modulo(pid(ip),nf_global)
        ipos=ipos+1

        dx=ipos-0.5-qpos_mean
        dx=modulo(dx+nf_global/2,real(nf_global))-nf_global/2
        dx_mean=dx_mean+dx
        !CDM ===================================================================
        vq(1)=phi(ipos(1)+1,ipos(2),ipos(3))-phi(ipos(1)-1,ipos(2),ipos(3))
        vq(2)=phi(ipos(1),ipos(2)+1,ipos(3))-phi(ipos(1),ipos(2)-1,ipos(3))
        vq(3)=phi(ipos(1),ipos(2),ipos(3)+1)-phi(ipos(1),ipos(2),ipos(3)-1)
        vq=-vq/(8*pi)/Dgrow(1/(1+z_checkpoint(1)))*vf
        spin_q(1)=spin_q(1)+dx(2)*vq(3)-dx(3)*vq(2)
        spin_q(2)=spin_q(2)+dx(3)*vq(1)-dx(1)*vq(3)
        spin_q(3)=spin_q(3)+dx(1)*vq(2)-dx(2)*vq(1)
        force_cdm=force_cdm+vq
        !neutrinos
        vq(1)=phu(ipos(1)+1,ipos(2),ipos(3))-phu(ipos(1)-1,ipos(2),ipos(3))
        vq(2)=phu(ipos(1),ipos(2)+1,ipos(3))-phu(ipos(1),ipos(2)-1,ipos(3))
        vq(3)=phu(ipos(1),ipos(2),ipos(3)+1)-phu(ipos(1),ipos(2),ipos(3)-1)
        vq=-vq/(8*pi)/Dgrow(1/(1+z_checkpoint(1)))*vf
        spin_u(1)=spin_u(1)+dx(2)*vq(3)-dx(3)*vq(2)
        spin_u(2)=spin_u(2)+dx(3)*vq(1)-dx(1)*vq(3)
        spin_u(3)=spin_u(3)+dx(1)*vq(2)-dx(2)*vq(1)
        force_neu=force_neu+vq

        ! initial inertia
        inertia(1,1)=inertia(1,1)+dx(1)**2
        inertia(2,2)=inertia(2,2)+dx(2)**2
        inertia(3,3)=inertia(3,3)+dx(3)**2
        inertia(1,2)=inertia(1,2)+dx(1)*dx(2)
        inertia(2,3)=inertia(2,3)+dx(2)*dx(3)
        inertia(3,1)=inertia(3,1)+dx(3)*dx(1)

        ! initial tidal tensor
!ipos=ceiling(qpos_mean)
        tide(1,1)=tide(1,1)+phi(ipos(1)+1,ipos(2),ipos(3))-2*phi(ipos(1),ipos(2),ipos(3))+phi(ipos(1)-1,ipos(2),ipos(3))
        tide(2,2)=tide(2,2)+phi(ipos(1),ipos(2)+1,ipos(3))-2*phi(ipos(1),ipos(2),ipos(3))+phi(ipos(1),ipos(2)-1,ipos(3))
        tide(3,3)=tide(3,3)+phi(ipos(1),ipos(2),ipos(3)+1)-2*phi(ipos(1),ipos(2),ipos(3))+phi(ipos(1),ipos(2),ipos(3)-1)
        tide(1,2)=tide(1,2)+(phi(ipos(1)+1,ipos(2)+1,ipos(3))+phi(ipos(1)-1,ipos(2)-1,ipos(3))&
                            -phi(ipos(1)+1,ipos(2)-1,ipos(3))-phi(ipos(1)-1,ipos(2)+1,ipos(3)))/4
        tide(2,3)=tide(2,3)+(phi(ipos(1),ipos(2)+1,ipos(3)+1)+phi(ipos(1),ipos(2)-1,ipos(3)-1)&
                            -phi(ipos(1),ipos(2)+1,ipos(3)-1)-phi(ipos(1),ipos(2)-1,ipos(3)+1))/4
        tide(3,1)=tide(3,1)+(phi(ipos(1)+1,ipos(2),ipos(3)+1)+phi(ipos(1)-1,ipos(2),ipos(3)-1)&
                            -phi(ipos(1)+1,ipos(2),ipos(3)-1)-phi(ipos(1)-1,ipos(2),ipos(3)+1))/4

        tidu(1,1)=tidu(1,1)+phu(ipos(1)+1,ipos(2),ipos(3))-2*phu(ipos(1),ipos(2),ipos(3))+phu(ipos(1)-1,ipos(2),ipos(3))
        tidu(2,2)=tidu(2,2)+phu(ipos(1),ipos(2)+1,ipos(3))-2*phu(ipos(1),ipos(2),ipos(3))+phu(ipos(1),ipos(2)-1,ipos(3))
        tidu(3,3)=tidu(3,3)+phu(ipos(1),ipos(2),ipos(3)+1)-2*phu(ipos(1),ipos(2),ipos(3))+phu(ipos(1),ipos(2),ipos(3)-1)
        tidu(1,2)=tidu(1,2)+(phu(ipos(1)+1,ipos(2)+1,ipos(3))+phu(ipos(1)-1,ipos(2)-1,ipos(3))&
                            -phu(ipos(1)+1,ipos(2)-1,ipos(3))-phu(ipos(1)-1,ipos(2)+1,ipos(3)))/4
        tidu(2,3)=tidu(2,3)+(phu(ipos(1),ipos(2)+1,ipos(3)+1)+phu(ipos(1),ipos(2)-1,ipos(3)-1)&
                            -phu(ipos(1),ipos(2)+1,ipos(3)-1)-phu(ipos(1),ipos(2)-1,ipos(3)+1))/4
        tidu(3,1)=tidu(3,1)+(phu(ipos(1)+1,ipos(2),ipos(3)+1)+phu(ipos(1)-1,ipos(2),ipos(3)-1)&
                            -phu(ipos(1)+1,ipos(2),ipos(3)-1)-phu(ipos(1)-1,ipos(2),ipos(3)+1))/4


        spin_x=halo_info%ang_mom
      enddo
      dx_mean=dx_mean/nptemp ! check mean
      if(maxval(abs(dx_mean))>20) then
        stop "dx incorrect"
      endif

      inertia(2,1)=inertia(1,2); inertia(3,2)=inertia(2,3); inertia(1,3)=inertia(3,1)
      tide(2,1)=tide(1,2); tide(3,2)=tide(2,3); tide(1,3)=tide(3,1)
      tidu(2,1)=tidu(1,2); tidu(3,2)=tidu(2,3); tidu(1,3)=tidu(3,1)
      tide=tide/nptemp
      tidu=tidu/nptemp
      write(16) inertia,tidu
      write(14) tide ! write T_c into file

#     ifdef Voronoi
      !read(19) tide
      !torque_v=matmul(inert_v,tide)
      torque_v=matmul(inert_v(:,:,ihalo),tide)
      spin_vor(1)=-torque_v(2,3)+torque_v(3,2); spin_vor(2)=-torque_v(3,1)+torque_v(1,3); spin_vor(3)=-torque_v(1,2)+torque_v(2,1)
#     endif
      torque=matmul(inertia,tide)
      torquu=matmul(inertia,tidu)
      torqur=matmul(tide,tidu)
      spin_t(1)=-torque(2,3)+torque(3,2); spin_t(2)=-torque(3,1)+torque(1,3); spin_t(3)=-torque(1,2)+torque(2,1)
      spin_v(1)=-torquu(2,3)+torquu(3,2); spin_v(2)=-torquu(3,1)+torquu(1,3); spin_v(3)=-torquu(1,2)+torquu(2,1)
      spin_r(1)=-torqur(2,3)+torqur(3,2); spin_r(2)=-torqur(3,1)+torqur(1,3); spin_r(3)=-torqur(1,2)+torqur(2,1)
      !print*, 'e Tc Tnu =',spin_r
      !print*,tide
      !print*,sum(Tc3(:,:,:,ihalo),3)
#ifdef eig
      do ieig=1,3
        !torqur=matmul(Tc3(:,:,ieig,ihalo),tidu)
        torqur=matmul(Tc3(:,:,ieig,ihalo),tide)
        !torqur=matmul(sum(Tc3(:,:,:,ihalo),3),tidu)
        spin_e(1,ieig)=-torqur(2,3)+torqur(3,2)
        spin_e(2,ieig)=-torqur(3,1)+torqur(1,3)
        spin_e(3,ieig)=-torqur(1,2)+torqur(2,1)
      enddo
      torqur=matmul(sum(Tc3(:,:,:,ihalo),3),tidu)
      spin_r(1)=-torqur(2,3)+torqur(3,2)
      spin_r(2)=-torqur(3,1)+torqur(1,3)
      spin_r(3)=-torqur(1,2)+torqur(2,1)
      !print*,'e sum(Tc_eig) Tnu',spin_r
      !print*,'sum(e Tc_eig Tnu)',sum(spin_e,2)
      !print*,'corr in spin_e'
      !print*,sum(spin_e(:,1)*spin_e(:,2))/sqrt(sum(spin_e(:,1)**2))/sqrt(sum(spin_e(:,2)**2))
      !print*,sum(spin_e(:,2)*spin_e(:,3))/sqrt(sum(spin_e(:,2)**2))/sqrt(sum(spin_e(:,3)**2))
      !print*,sum(spin_e(:,3)*spin_e(:,1))/sqrt(sum(spin_e(:,3)**2))/sqrt(sum(spin_e(:,1)**2))
      !print*,'corr of spin_e - spin_r'
      !print*,sum(spin_e(:,1)*spin_r)/sqrt(sum(spin_e(:,1)**2))/sqrt(sum(spin_r**2))
      !print*,sum(spin_e(:,2)*spin_r)/sqrt(sum(spin_e(:,2)**2))/sqrt(sum(spin_r**2))
      !print*,sum(spin_e(:,3)*spin_r)/sqrt(sum(spin_e(:,3)**2))/sqrt(sum(spin_r**2))
#endif

      write(13) spin_q,spin_t,spin_x,spin_u,spin_v,spin_r,spin_e(:,1:3)
      theta(0,ihalo)=halonp(ihalo)
      theta(1,ihalo)=sum(spin_q*spin_x)/sqrt(sum(spin_q**2))/sqrt(sum(spin_x**2))
      theta(2,ihalo)=sum(spin_q*spin_t)/sqrt(sum(spin_q**2))/sqrt(sum(spin_t**2))
      theta(3,ihalo)=sum(spin_t*spin_x)/sqrt(sum(spin_t**2))/sqrt(sum(spin_x**2))
      theta(4,ihalo)=sum(spin_q*spin_u)/sqrt(sum(spin_q**2))/sqrt(sum(spin_u**2))

      theta(5,ihalo)=sum(spin_u*spin_v)/sqrt(sum(spin_u**2))/sqrt(sum(spin_v**2))
      theta(6,ihalo)=sum(spin_v*spin_r)/sqrt(sum(spin_v**2))/sqrt(sum(spin_r**2))
      theta(7,ihalo)=sum(spin_u*spin_r)/sqrt(sum(spin_u**2))/sqrt(sum(spin_r**2))

#ifdef eig
      do ieig=1,3
        theta(8+ieig,ihalo)=sum(spin_q*spin_e(:,ieig))/sqrt(sum(spin_q**2))/sqrt(sum(spin_e(:,ieig)**2))
      enddo
      !print*, 'theta_vr =',theta_vr(ihalo)
      !print*, 'theta_ve =',theta_ve(:,ihalo)
      !print*, sum(spin_v*sum(spin_e(:,:2),2))/sqrt(sum(spin_v**2))/sqrt(sum(sum(spin_e(:,:2),2)**2))
      !read(*,*)
#endif

      torque_c=torque_c+sqrt(sum(spin_q**2))
      torque_u=torque_u+sqrt(sum(spin_u**2))
      theta(8,ihalo)=sum(force_cdm*force_neu)/sqrt(sum(force_cdm**2))/sqrt(sum(force_neu**2))
#ifdef Voronoi
      theta(12,ihalo)=sum(spin_vor*spin_q)/sqrt(sum(spin_vor**2))/sqrt(sum(spin_q**2))
#endif
      if (ihalo==1) then
        print*,'halo',ihalo
        print*,'inertia'
        print*,inertia
        print*,'inert_v'
        print*,inert_v(:,:,ihalo)
        print*,'tide'
        print*,tide
        print*,'spin'
        print*,spin_t
        print*,'spin_vor'
        print*,spin_vor
      endif

    enddo

    print*,'mean qx correlation =',sum(theta(1,:nhalo))/nhalo
    print*,'mean qt correlation =',sum(theta(2,:nhalo))/nhalo
    print*,'mean tx correlation =',sum(theta(3,:nhalo))/nhalo
    print*,'mean vor correlation =',sum(theta(12,:nhalo))/nhalo
    print*,''
    print*,'mean qu correlation =',sum(theta(4,:nhalo))/nhalo
    print*,''
    print*,'mean uv correlation =',sum(theta(5,:nhalo))/nhalo
    print*,'mean vr correlation =',sum(theta(6,:nhalo))/nhalo
    print*,'mean ur correlation =',sum(theta(7,:nhalo))/nhalo
    print*,''
    print*,'mean ve correlation =',sum(theta(9:11,:nhalo),2)/nhalo
    print*,'mean ff correlation =',sum(theta(8,:nhalo))/nhalo

    print*,'torque_c/torque_u=',torque_u/torque_c

    close(11);close(12);close(13);close(14);close(16)

    isort_mass(:nhalo)=[(i0,i0=1,nhalo)]
    call indexedsort(nhalo,halonp,isort_mass)
    do i0=0,12
      theta(i0,:)=theta(i0,isort_mass(nhalo:1:-1))
    enddo
    open(16,file=output_name('halo_correlation'),status='replace',access='stream')
    write(16) theta
    close(16)
    open(16,file=output_name('halo_sort_index'),status='replace',access='stream')
    write(16) isort_mass(nhalo:1:-1)
    close(16)
  enddo
  sync all

  contains
  function vfactor(a)
    implicit none
    real, parameter :: n_p = -(1./4.)+(5./4.)*sqrt(1-24.*omega_nu/omega_m/25.) !~1-3f/5
    real :: a
    real :: H,km,lm
    real :: vfactor
    lm=omega_l/omega_m
    km=(1-omega_m-omega_l)/omega_m
    H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
    vfactor=a**2*H
    vfactor=vfactor*n_p!(1.-3.*(omega_nu/omega_m)/5.)
  endfunction vfactor
  function Dgrow(a)
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real a
    real Dgrow
    real g,ga,hsq,oma,ola
    hsq=om/a**3+(1-om-ol)/a**2+ol
    oma=om/(a**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
    Dgrow=a*ga/g
  end function Dgrow
end

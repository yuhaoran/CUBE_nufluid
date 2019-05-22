program forcepower
  use parameters
  use pencil_fft
  use powerspectrum
  implicit none
  save

  integer i
  real scale_fac(1000),xi(10,nbin)[*],xi0(10,222)
  real growth(2,1000)
  character(10) :: str_a

  call geometry
  if (head) then
    print*,''
    print*,'coarce force power on resolution'
    print*,'ng =',ng
    print*,'ng*nn =',ng*nn
    print*,''
    print*,'checkpoint at:'
    open(16,file=output_dir()//'zFc'//output_suffix(),status='old')
    do i=1,1000
      read(16,end=71,fmt='(f10.5)') scale_fac(i)
      write(str_a,'(f7.3)') scale_fac(i)
      print*,int(i,2),scale_fac(i),str_a
    enddo
    71 n_checkpoint=i-1
    close(16)
    print*,''
  endif

  call create_penfft_plan
  sync all
  open(11,file=opath//'image1/5.000_power_Fx_1.bin',status='old',access='stream')
  read(11) xi0
  close(11)
  do cur_checkpoint=1,n_checkpoint
    write(str_a,'(f7.3)') scale_fac(cur_checkpoint)
    if (head) print*, 'Start analyzing scale factor a =',str_a
    print*,opath//'image1/'//trim(adjustl(str_a))//'_Fc'//output_suffix()
    open(11,file=opath//'image1/'//trim(adjustl(str_a))//'_Fc'//output_suffix(),status='old',access='stream')
    read(11) r3
    close(11)

    call cross_power(xi,r3,r3)
    if (head) then
      open(15,file=opath//'image1/'//trim(adjustl(str_a))//'_powerFcx'//output_suffix(),status='replace',access='stream')
      write(15) xi
      close(15)
    endif
    print*,'wrote',opath//'image1/'//trim(adjustl(str_a))//'_powerFcx'//output_suffix()

    growth(1,cur_checkpoint)=scale_fac(cur_checkpoint)
    growth(2,cur_checkpoint)=sqrt(sum(xi(3,:4))/sum(xi0(3,:4)))
    print*, growth(:,cur_checkpoint)
  enddo
  call destroy_penfft_plan
  sync all


  growth(2,:)=growth(2,:)/1.0408
  growth(1,n_checkpoint+1)=1;
  growth(2,n_checkpoint+1)=4.851;

  open(15,file=opath//'image1/growth_z5'//output_suffix(),status='replace',access='stream')
  write(15) growth(:,1)*0
  write(15) growth(:,:n_checkpoint+1)
  close(15)
  print*,'wrote ',opath//'image1/growth_z5'//output_suffix()
  if (head) print*,'forcepower done'


endprogram

!! module for power spectrum analysis
! uses pencil_fft however work for single image only
! cx1,cx2 can be memory optimized
! auto-power can be memory optimized
! check nexp frequently
#define linear_kbin
!#define pl2d

module powerspectrum
use pencil_fft

#ifdef linear_kbin
  integer(8),parameter :: nbin=nint(nyquest*sqrt(3.))
#else
  integer(8),parameter :: nbin=floor(4*log(nyquest*sqrt(3.)/0.95)/log(2.))
#endif
#ifdef pl2d
  integer kp,kl
  integer nmode(nint(nyquest*sqrt(2.))+1,nyquest+1)
  real pow2d(nint(nyquest*sqrt(2.))+1,nyquest+1)
  real pow2drsd(nint(nyquest*sqrt(2.))+1,nyquest+1)
#endif
complex cx1(ng*nn/2+1,ng,npen),cx2(ng*nn/2+1,ng,npen)


contains

subroutine cross_power(xi,cube1,cube2)
  implicit none

  integer(8) i,j,k,ig,jg,kg,ibin
  real kr,kx(3),sincx,sincy,sincz,sinc,rbin

  real cube1(ng,ng,ng),cube2(ng,ng,ng)
  real xi(10,nbin)[*]
  real amp11,amp12,amp22
  complex cx1(ng*nn/2+1,ng,npen),cx2(ng*nn/2+1,ng,npen)

  real,parameter :: nexp=4.0 ! CIC kernel

  xi=0

  r3=cube1
  call pencil_fft_forward
  cx1=cxyz/ng_global/ng_global/ng_global

  r3=cube2
  call pencil_fft_forward
  cx2=cxyz/ng_global/ng_global/ng_global

  xi=0
  sync all
#ifdef pl2d
  print*, 'size of pow2d',nint(nyquest*sqrt(2.))+1,nyquest+1
  pow2d=0
  nmode=0
#endif

  do k=1,npen
  do j=1,ng
  do i=1,nyquest+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*ng+j
    ig=i
    kx=mod((/ig,jg,kg/)+nyquest-1,ng_global)-nyquest
    if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
    if ((ig==1.or.ig==ng*nn/2+1) .and. jg>ng*nn/2+1) cycle
    if ((ig==1.or.ig==ng*nn/2+1) .and. (jg==1.or.jg==ng*nn/2+1) .and. kg>ng*nn/2+1) cycle
    kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
    sincx=merge(1.0,sin(pi*kx(1)/ng_global)/(pi*kx(1)/ng_global),kx(1)==0.0)
    sincy=merge(1.0,sin(pi*kx(2)/ng_global)/(pi*kx(2)/ng_global),kx(2)==0.0)
    sincz=merge(1.0,sin(pi*kx(3)/ng_global)/(pi*kx(3)/ng_global),kx(3)==0.0)
    sinc=sincx*sincy*sincz
#   ifdef linear_kbin
      ibin=nint(kr)
#   else
      rbin=4.0/log(2.)*log(kr/0.95)
      ibin=merge(ceiling(rbin),floor(rbin),rbin<1)
#   endif
    xi(1,ibin)=xi(1,ibin)+1 ! number count
    xi(2,ibin)=xi(2,ibin)+kr ! k count
    amp11=real(cx1(i,j,k)*conjg(cx1(i,j,k)))/(sinc**4.0)*4*pi*kr**3
    amp22=real(cx2(i,j,k)*conjg(cx2(i,j,k)))/(sinc**4.0)*4*pi*kr**3
    amp12=real(cx1(i,j,k)*conjg(cx2(i,j,k)))/(sinc**4.0)*4*pi*kr**3
#ifdef pl2d
    kp=nint(sqrt(kx(1)**2+kx(3)**2))+1
    kl=abs(kx(2))+1
    nmode(kp,kl)=nmode(kp,kl)+1
    pow2d(kp,kl)=pow2d(kp,kl)+amp11
    pow2drsd(kp,kl)=pow2drsd(kp,kl)+amp22
#endif
    xi(3,ibin)=xi(3,ibin)+amp11 ! auto power 1
    xi(4,ibin)=xi(4,ibin)+amp22 ! auto power 2
    xi(5,ibin)=xi(5,ibin)+amp12 ! cross power
    xi(6,ibin)=xi(6,ibin)+1/sinc**2.0 ! kernel 1
    xi(7,ibin)=xi(7,ibin)+1/sinc**4.0 ! kernel 2

  enddo
  enddo
  enddo
  sync all
#ifdef pl2d
  nmode=max(1,nmode)
  pow2d=pow2d/nmode
  pow2drsd=pow2drsd/nmode
  open(55,file=output_name('pow2d'),status='replace',access='stream')
  write(55) pow2d
  write(55) pow2drsd
  close(55)
#endif

  ! co_sum
  if (head) then
    do i=2,nn**3
      xi=xi+xi(:,:)[i]
    enddo
  endif
  sync all

  ! broadcast
  xi=xi(:,:)[1]
  sync all

  ! divide and normalize
  xi(2,:)=xi(2,:)/xi(1,:)*(2*pi)/box ! k_phy
  xi(3,:)=xi(3,:)/xi(1,:) ! Delta_LL
  xi(4,:)=xi(4,:)/xi(1,:) ! Delta_RR
  xi(5,:)=xi(5,:)/xi(1,:) ! Delta_LR ! cross power
  xi(6,:)=xi(6,:)/xi(1,:) ! kernel
  xi(7,:)=xi(7,:)/xi(1,:) ! kernel
  xi(8,:)=xi(5,:)/sqrt(xi(3,:)*xi(4,:)) ! r
  xi(9,:)=sqrt(xi(4,:)/xi(3,:)) ! b
  xi(10,:)=xi(8,:)**4/xi(9,:)**2 * xi(4,:) ! P_RR*r^4/b^2 reco power
  sync all
endsubroutine

endmodule

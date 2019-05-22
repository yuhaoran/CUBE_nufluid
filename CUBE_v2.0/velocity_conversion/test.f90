use parameters
integer,parameter :: nvbin=2**(8*izipv)
integer,parameter :: nf_global=nf*nn
real vsim2phys

real vtest
integer(2) itest

real vdisp(506,2)
open(23,file='vdisp.bin',access='stream')
read(23) vdisp
close(23)
a=1
vsim2phys=(1/a)*1.5*box*h0*sqrt(omega_m)/real(nf_global)
! convert to simulation unit, divide by sqrt(DoF=3)
vdisp(:,2)=vdisp(:,2)*vdisp(:,1)/1.5/box/h0/sqrt(omega_m)*real(nf_global)/sqrt(3.)

print*,'nvbin',nvbin
print*,'vsim2phys',vsim2phys

print*, vdisp(10,:)
print*, vdisp(100,:)
print*, vdisp(300,:)
print*, vdisp(506,:)

print*, 'interp_vdisp test',interp_vdisp(0.991)

vtest=5.
sigma_vi=interp_vdisp(1.)

print*,'vtest',vtest
itest=nint((nvbin-1)*atan(sqrt(pi/2)/sigma_vi*vtest)/pi)
print*,'itest',itest
vtest=tan(pi*itest/(nvbin-1))/(sqrt(pi/2)/sigma_vi);
print*,'vtest',vtest

contains
  real function interp_vdisp(aa)
  implicit none
  integer ii,i1,i2
  real aa
  i1=1
  i2=506
  do while (i2-i1>1)
    ii=(i1+i2)/2
    if (aa>vdisp(ii,1)) then
      i1=ii
    else
      i2=ii
    endif
  enddo
  interp_vdisp=vdisp(i1,2)+(vdisp(i2,2)-vdisp(i1,2))*(aa-vdisp(i1,1))/(vdisp(i2,1)-vdisp(i1,1))
  endfunction

end

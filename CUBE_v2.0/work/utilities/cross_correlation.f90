program ccc
! compute cross correlation coefficient
use penfft_fine
use powerspectrum
implicit none

!integer,parameter :: ng=96

real xi(10,nbin)[*]
real delta_nbody(ng,ng,ng)
real delta_voronoi(ng,ng,ng)

print*, 'ng =',ng
call geometry
call create_penfft_fine_plan

open(11,file='../output/universe1/node0/0.000delta_nbody.dat',status='old',action='read',access='stream')
read(11) delta_nbody
close(11)

open(11,file='../output/universe1/node0/0.000delta_voronoi.dat',status='old',action='read',access='stream')
read(11) delta_voronoi
close(11)

call cross_power(xi,delta_nbody,delta_voronoi)
open(15,file='.'//opath//'node'//image2str(this_image()-1)//'/xi_NV.dat',status='replace',access='stream')
write(15) xi
close(15)

call destroy_penfft_fine_plan

end

module neutrinos
  use parameters
  implicit none
  save

#ifdef NEUTRINOS
  !parameters
  real,parameter :: image_buffer_nu=1.2
  real,parameter :: tile_buffer_nu=3.0
  integer(8),parameter :: np_image_nu=(nc*np_nc_nu)**3 ! average number of particles per image
  integer(8),parameter :: np_image_max_nu=np_image_nu*(nte*1./nt)**3*image_buffer_nu
  integer(8),parameter :: np_tile_max_nu=np_image_nu/nnt**3*(nte*1./nt)**3*tile_buffer_nu

  integer(izipx_nu) xp_nu(3,np_image_max_nu)[*], xp_new_nu(3,np_tile_max_nu)
  integer(izipv_nu) vp_nu(3,np_image_max_nu)[*], vp_new_nu(3,np_tile_max_nu)
  integer(4) rhoc_nu(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
  integer(8) cum_nu(1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt)[*]
  real(4) vfield_nu(3,1-ncb:nt+ncb,1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,nnt) ! cannot have >7 dims
  integer(izipi) pid_nu(np_image_max_nu)[*], pid_new_nu(np_tile_max_nu)

  integer(8) npcheck_nu
#endif

endmodule

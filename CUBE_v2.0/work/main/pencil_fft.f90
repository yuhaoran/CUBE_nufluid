! to be optimized - image1d, ctransfer more than one slab
module pencil_fft
  use parameters
  implicit none
  save

  integer(4),parameter :: NULL=0
  integer(8) planx,plany,planz,iplanx,iplany,iplanz

  ! fft arrays
  real        r3(ng,ng,ng),r0(ng,ng)
  complex     c3(ng/2,ng,ng)
  real        rxyz(ng_global+2  ,ng,npen)
  complex     cxyz(ng_global/2+1,ng,npen)
  complex     cyyyxz(npen,nn,nn,ng/2+1,npen)
  complex     cyyxz(ng,     nn,ng/2+1,npen)
  complex     czzzxy(npen,nn,nn,ng/2+1,npen)

  ! communication coarrays
  complex ctransfer1(ng/2,ng,nn)[*]
  complex ctransfer2(ng,ng/2+1,nn)[*]
  complex ctransfer3(npen,npen,nn,nn)[*]
  complex ctransfer4(ng/2+1,ng,nn)[*]

  ! equivalence statements
  equivalence(cyyyxz,cyyxz,r3,c3)
  equivalence(czzzxy,rxyz,cxyz)

  contains

  subroutine pencil_fft_forward
    implicit none
    save
    sync all
    call c2x
    call sfftw_execute(planx)
    call x2y
    call sfftw_execute(plany)
    call y2z
    call sfftw_execute(planz)
    call z2y
    call y2x
    sync all
  endsubroutine

  subroutine pencil_fft_backward
    implicit none
    save
    sync all
    call x2y
    call y2z
    call sfftw_execute(iplanz)
    call z2y
    call sfftw_execute(iplany)
    call y2x
    call sfftw_execute(iplanx)
    call x2c
    r3=r3/ng_global/ng_global/ng_global
    sync all
  endsubroutine

  subroutine c2x
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,npen ! loop over cells in z, extract slabs
      ctransfer1(:,:,1:nn)=c3(:,:,islab::npen) ! nn slabs of c3 copied to ctransfer1
      sync all
      do i1=1,nn ! loop over parts in x, get slabs from each y node
        ! i1=mod()
        cxyz(ng*(i1-1)/2+1:ng*i1/2,:,islab)=ctransfer1(:,:,m2)[image1d(i1,m1,m3)]
      enddo
      sync all
    enddo
  endsubroutine

  subroutine x2y
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,npen ! loop over z
      do i1=1,nn ! loop over squares in x direction
        ctransfer2(:,:,i1)=transpose(cxyz(ng/2*(i1-1)+1:ng/2*i1+1,:,islab))
      enddo
      sync all
      do i1=1,nn
        cyyxz(:,i1,:,islab)=ctransfer2(:,:,m1)[image1d(i1,m2,m3)]
      enddo
      sync all
    enddo
  endsubroutine

  subroutine y2z
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,ng/2+1 ! loop over slices in x direction
      do i2=1,nn
      do i1=1,nn
        ctransfer3(:,:,i1,i2)=transpose(cyyyxz(:,i1,i2,islab,:))
      enddo
      enddo
      sync all
      do i2=1,nn
      do i1=1,nn
        czzzxy(:,i1,i2,islab,:)=ctransfer3(:,:,m2,m3)[image1d(m1,i1,i2)]
      enddo
      enddo
      sync all
    enddo
  endsubroutine

  subroutine z2y
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,ng/2+1 ! loop over slices in x direction
      do i2=1,nn
      do i1=1,nn
        ctransfer3(:,:,i1,i2)=transpose(czzzxy(:,i1,i2,islab,:))
      enddo
      enddo
      sync all
      do i2=1,nn
      do i1=1,nn
        cyyyxz(:,i1,i2,islab,:)=ctransfer3(:,:,m2,m3)[image1d(m1,i1,i2)]
      enddo
      enddo
      sync all
    enddo
  endsubroutine

  subroutine y2x
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,npen ! loop over z
      do i1=1,nn ! loop over squares in x direction
        ctransfer4(:,:,i1)=transpose(cyyxz(:,i1,:,islab))
      enddo
      sync all
      do i1=1,nn
        cxyz(ng/2*(i1-1)+1:ng/2*i1+1,:,islab)=ctransfer4(:,:,m1)[image1d(i1,m2,m3)]
      enddo
      sync all
    enddo
  endsubroutine

  subroutine x2c
    implicit none
    save
    integer(8) i0,i1,i2,islab
    do islab=1,npen
      do i1=1,nn
        ctransfer1(:,:,i1)=cxyz(ng*(i1-1)/2+1:ng*i1/2,:,islab)
      enddo
      sync all
      do i1=1,nn
        c3(:,:,islab+(i1-1)*npen)=ctransfer1(:,:,m1)[image1d(m2,i1,m3)]
      enddo
      sync all
    enddo
  endsubroutine

  subroutine create_penfft_plan
    implicit none
    save
    include 'fftw3.f'
    call sfftw_plan_many_dft_r2c(planx,1,ng*nn,ng*npen,cxyz,NULL,1,ng*nn+2,cxyz,NULL,1,ng*nn/2+1,FFTW_MEASURE)
    call sfftw_plan_many_dft_c2r(iplanx,1,ng*nn,ng*npen,cxyz,NULL,1,ng*nn/2+1,cxyz,NULL,1,ng*nn+2,FFTW_MEASURE)
    call sfftw_plan_many_dft(plany,1,ng*nn,(ng/2+1)*npen,cyyxz,NULL,1,ng*nn,cyyxz,NULL,1,ng*nn,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(iplany,1,ng*nn,(ng/2+1)*npen,cyyxz,NULL,1,ng*nn,cyyxz,NULL,1,ng*nn,FFTW_BACKWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(planz,1,ng*nn,(ng/2+1)*npen,czzzxy,NULL,1,ng*nn,czzzxy,NULL,1,ng*nn,FFTW_FORWARD,FFTW_MEASURE)
    call sfftw_plan_many_dft(iplanz,1,ng*nn,(ng/2+1)*npen,czzzxy,NULL,1,ng*nn,czzzxy,NULL,1,ng*nn,FFTW_BACKWARD,FFTW_MEASURE)
    sync all
  endsubroutine

  subroutine destroy_penfft_plan
    implicit none
    save
    include 'fftw3.f'
    call sfftw_destroy_plan(planx)
    call sfftw_destroy_plan(iplanx)
    call sfftw_destroy_plan(plany)
    call sfftw_destroy_plan(iplany)
    call sfftw_destroy_plan(planz)
    call sfftw_destroy_plan(iplanz)
    sync all
  endsubroutine
endmodule

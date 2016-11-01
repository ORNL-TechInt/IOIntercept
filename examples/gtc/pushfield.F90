subroutine pushfield
  use global_parameters
  use field_array
  implicit none

  integer m,i,j,ij
  real(wp) dtime,b,r,theta

  if(irk==1)then

! 1st step of Runge-Kutta method
    dtime=0.5*tstep

!$omp parallel do private(m)
    do m=1,mgrid
        apara0(1,m)=apara(1,m)
        fluidne0(1,m)=fluidne(1,m)
    enddo

! 2nd step of Runge-Kutta method
  else
    dtime=tstep
  endif

! advance del_ne and a_parallel
!$omp parallel do private(i,j,ij,r,theta,b)
  do i=0,mpsi
     r=a0+deltar*real(i)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        theta=deltat(i)*real(j)+zeta1*qtinv(i)
        b=1.0/(1.0+r*cos(theta)) !toroidal
        fluidne(1,ij) = fluidne0(1,ij)-dtime*b*b*gradue(3,1,ij)
        apara(1,ij) = apara0(1,ij)+dtime*b*gradind(3,1,ij)
     enddo
  enddo

! poloidal and toroidal periodicity conditions
  call periodicity(fluidne)
  call periodicity(apara)

end subroutine pushfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine periodicity(fgrid)
  use global_parameters
  use field_array
  implicit none

  integer i,ii,jt,icount,idest,isource,isendtag,irecvtag,ierror,istatus(MPI_STATUS_SIZE)
  real(wp) sendr(mgrid),recvl(mgrid)
  real(wp),dimension(0:1,mgrid) :: fgrid

! toroidal end point, send E to right and receive from left
!$omp parallel do private(i)
  do i=1,mgrid
     sendr(i)=fgrid(1,i)
     recvl(i)=0.0
  enddo
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
      recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

! unpack end point data for k=0
  if(myrank_toroidal==0)then
!$omp parallel do private(i,ii,jt)
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        fgrid(0,ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
     enddo
  else
!$omp parallel do private(i)
     do i=1,mgrid
        fgrid(0,i)=recvl(i)
     enddo
  endif
  
! poloidal end point
!omp parallel do private(i,ii,jt)
  do i=0,mpsi
     ii=igrid(i)
     jt=mtheta(i)
     fgrid(0:1,ii)=fgrid(0:1,ii+jt)
  enddo
  
end subroutine periodicity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

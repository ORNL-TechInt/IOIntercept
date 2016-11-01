subroutine smooth(farray)

  use global_parameters
  use field_array
  implicit none

  integer i,j,k,ii,ij,ip,jt,kz,ismooth,mz,mzmax,mzbig,jpe,indp,indt,&
      indp1,indt1,meachtheta,jtp,icount,ierror,idest,isource,isendtag,l,&
      irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),phism(mgrid),&
      ptemp,pleft(mthetamax),pright(mthetamax),phitmp(mgrid),&
      wt,r,dt,wz,zdum,tdum,pi2_inv,farray(0:1,mgrid)

!$omp parallel do private(i)
  do i=1,mgrid
     phitmp(i)=farray(1,i)
  enddo

  ismooth=1
  if(nfilter==1)ismooth=0
  do ip=1,ismooth

! radial smoothing
    do l=1,2	
        do i=1,mpsi-1
          phitmp(igrid(i))=phitmp(igrid(i)+mtheta(i))
        enddo
        phism=0.0
!$omp parallel do private(i,j,ij)        
        do i=1,mpsi-1
          do j=1,mtheta(i)
              ij=igrid(i)+j
              phism(ij)=0.25*((1.0-wtp1(1,ij))*phitmp(jtp1(1,ij))+&
                  wtp1(1,ij)*phitmp(jtp1(1,ij)+1)+&
                  (1.0-wtp1(2,ij))*phitmp(jtp1(2,ij))+&
                  wtp1(2,ij)*phitmp(jtp1(2,ij)+1))-&
                  0.0625*((1.0-wtp2(1,ij))*phitmp(jtp2(1,ij))+&
                  wtp2(1,ij)*phitmp(jtp2(1,ij)+1)+&
                  (1.0-wtp2(2,ij))*phitmp(jtp2(2,ij))+&
                  wtp2(2,ij)*phitmp(jtp2(2,ij)+1)) 
          enddo
        enddo
        phitmp=0.625*phitmp+phism

! poloidal smoothing (-0.0625 0.25 0.625 0.25 -0.0625)
!$omp parallel do private(i,ii,jt,pright)        
        do i=1,mpsi-1
          ii=igrid(i)
          jt=mtheta(i)
          pright(1:jt)=phitmp(ii+1:ii+jt)
          phitmp(ii+1:ii+jt)=0.625*pright(1:jt)+&
              0.25*(cshift(pright(1:jt),-1)+cshift(pright(1:jt),1))-&
              0.0625*(cshift(pright(1:jt),-2)+cshift(pright(1:jt),2))
        enddo
    enddo

! parallel smoothing 
! send phi to right and receive from left
    sendr=phitmp
    recvl=0.0
    icount=mgrid
    idest=right_pe
    isource=left_pe
    isendtag=myrank_toroidal
    irecvtag=isource
    call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,recvl,icount,&
          mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
     
! send phi to left and receive from right
    sendl=phitmp
    recvr=0.0
    idest=left_pe
    isource=right_pe
    isendtag=myrank_toroidal
    irecvtag=isource
    call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
          mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
    
!$omp parallel do private(i,ii,jt,j,ij,ptemp,pleft,pright)
    do i=1,mpsi-1
        ii=igrid(i)
        jt=mtheta(i)
        if(myrank_toroidal==0)then !down-shift for zeta=0
          pleft(1:jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
          pright(1:jt)=recvr(ii+1:ii+jt)
        elseif(myrank_toroidal==mtoroidal-1)then !up-shift for zeta=2*pi
          pright(1:jt)=cshift(recvr(ii+1:ii+jt),itran(i))
          pleft(1:jt)=recvl(ii+1:ii+jt)
        else
          pleft(1:jt)=recvl(ii+1:ii+jt)
          pright(1:jt)=recvr(ii+1:ii+jt)
        endif
        
        do j=1,mtheta(i)
          ij=igrid(i)+j
          ptemp=phitmp(ij)
          phitmp(ij)=0.5*ptemp+0.25*(pleft(j)+pright(j))
        enddo
    enddo
  enddo

! toroidal BC: send phi to right and receive from left
  sendr=phitmp
  recvl=0.0
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
      recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
     
  if(myrank_toroidal==0)then
    do i=1,mpsi-1
        ii=igrid(i)
        jt=mtheta(i)
        phism(ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
    enddo
  else
    phism=recvl
  endif

! poloidal BC
  do i=1,mpsi-1
    phitmp(igrid(i))=phitmp(igrid(i)+mtheta(i))
    phism(igrid(i))=phism(igrid(i)+mtheta(i))
  enddo

! radial boundary
  phitmp(igrid(0):igrid(0)+mtheta(0))=0.0
  phitmp(igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0
  phism(igrid(0):igrid(0)+mtheta(0))=0.0
  phism(igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0

!$omp parallel do private(i)
  do i=1,mgrid
     farray(1,i)=phitmp(i)
     farray(0,i)=phism(i)
  enddo

end subroutine smooth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initiate radial interpolation for poloidal grids
subroutine rinterpolation
  use global_parameters
  use field_array
  implicit none

  integer i,ip,j,indp,indt,ij,jt
  real tdum,wt

! allocate memory
  allocate(wtp1(2,mgrid),wtp2(2,mgrid),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
  allocate(jtp1(2,mgrid),jtp2(2,mgrid)) !---wj

!$omp parallel do private(i,ip,j,indp,indt,ij,tdum,jt,wt)
  do i=1,mpsi-1
     do ip=1,2
        indp=min(mpsi,i+ip)
        indt=max(0,i-ip)
        do j=1,mtheta(i)
           ij=igrid(i)+j
! upward
           tdum=(real(j)*deltat(i)+zeta1*(qtinv(i)-qtinv(indp)))/deltat(indp)
           jt=floor(tdum)
           wt=tdum-real(jt)
           jt=mod(jt+mtheta(indp),mtheta(indp))
           if(ip==1)then
              wtp1(1,ij)=wt
              jtp1(1,ij)=igrid(indp)+jt
           else
              wtp2(1,ij)=wt
              jtp2(1,ij)=igrid(indp)+jt
           endif
! downward         
           tdum=(real(j)*deltat(i)+zeta1*(qtinv(i)-qtinv(indt)))/deltat(indt)
           jt=floor(tdum)
           wt=tdum-real(jt)
           jt=mod(jt+mtheta(indt),mtheta(indt))
           if(ip==1)then
              wtp1(2,ij)=wt
              jtp1(2,ij)=igrid(indt)+jt
           else
              wtp2(2,ij)=wt
              jtp2(2,ij)=igrid(indt)+jt
           endif
        enddo
     enddo
  enddo
  
end subroutine rinterpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! d_phieff/d_t for fluid-kinetic hybrid electron model
subroutine dphieff
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer i,j
  real deltime

  deltime=max(1.0e-10,tstep)
  j=irk+2*(ihybrid-1)

!$omp parallel do private(i)
  do i=1,mgrid
     phit(:,i)=(phi(:,i)-phisave(:,i,j))/deltime
     phisave(:,i,j)=phi(:,i)
  enddo

end subroutine dphieff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chargei
  use global_parameters
  use particle_array
  use field_array
  implicit none
  
  integer m,i,im,j,k,ip,jt,ii,j11,j10,j01,j00,ierror,igyro,ij,ipjt,&
       icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) gyrodum,weight,rdum,tdum,r,b,upara,cmratio,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
       adum(0:mpsi),sendl(mgrid,2),recvr(mgrid,2),dnitmp(0:1,mgrid),djitmp(0:1,mgrid)
  
  cmratio=qion/aion
  densityi=0.0
  currenti=0.0
#ifdef _OPENMP
! The following lines are OpenMP directives for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).
! When compiling with the OpenMP option (-qsmp=omp for the XLF compiler on
! AIX and -mp for the MIPSpro compiler on IRIX), the preprocessor variable
! _OPENMP is automatically defined and we use it here to add pieces of codes
! needed to avoid contentions between threads. In the following parallel
! loop each thread has its own private array "dnitmp()" which is used to
! prevent the thread from writing into the shared array "densityi()"
!
!$omp parallel private(dnitmp,djitmp)
  dnitmp=0.   ! Set array elements to zero
  djitmp=0.   ! Set array elements to zero
#endif
  
  if(magnetic==0)then
! scatter ion density for electrostatic simulation
!$omp do private(m,igyro,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij)
     do m=1,mi
        weight=zion(5,m)
     
        wz1=weight*wzion(m)
        wz0=weight-wz1     
        do igyro=1,ngyroi
           wp1=wpion(igyro,m)
           wp0=1.0-wp1
           
           wt10=wp0*wtion0(igyro,m)
           wt00=wp0-wt10
           
           wt11=wp1*wtion1(igyro,m)
           wt01=wp1-wt11

#ifdef _OPENMP
! Use thread-private temporary array dnitmp to store the results
           ij=jtion0(igyro,m)
           dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt00
           dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt00
        
           ij=ij+1
           dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt10
           dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt10
           
           ij=jtion1(igyro,m)
           dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt01
           dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt01
           
           ij=ij+1
           dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt11
           dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt11
#else
! If no loop-level parallelism, use original algorithm (write directly
! into array "densityi()".
           ij=jtion0(igyro,m)
           densityi(0,ij) = densityi(0,ij) + wz0*wt00
           densityi(1,ij) = densityi(1,ij) + wz1*wt00
           
           ij=ij+1
           densityi(0,ij) = densityi(0,ij) + wz0*wt10
           densityi(1,ij) = densityi(1,ij) + wz1*wt10
           
           ij=jtion1(igyro,m)
           densityi(0,ij) = densityi(0,ij) + wz0*wt01
           densityi(1,ij) = densityi(1,ij) + wz1*wt01
           
           ij=ij+1
           densityi(0,ij) = densityi(0,ij) + wz0*wt11
           densityi(1,ij) = densityi(1,ij) + wz1*wt11
#endif
        enddo
     enddo

  else
! scatter ion density and current for electromagnetic simulation 
!$omp do private(m,igyro,r,b,upara,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij)
     do m=1,mi
        r=sqrt(2.0*zion(1,m))
        b=1.0/(1.0+r*cos(zion(2,m)))
        upara=zion(4,m)*b*cmratio
        weight=zion(5,m)
        
        wz1=weight*wzion(m)
        wz0=weight-wz1     
        do igyro=1,ngyroi
           wp1=wpion(igyro,m)
           wp0=1.0-wp1
           
           wt10=wp0*wtion0(igyro,m)
           wt00=wp0-wt10
           
           wt11=wp1*wtion1(igyro,m)
           wt01=wp1-wt11

#ifdef _OPENMP
! Use thread-private temporary array dnitmp to store the results
           ij=jtion0(igyro,m)
           dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt00
           dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt00
           djitmp(0,ij) = djitmp(0,ij) + wz0*wt00*upara
           djitmp(1,ij) = djitmp(1,ij) + wz1*wt00*upara
           
           ij=ij+1
           dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt10
           dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt10
           djitmp(0,ij) = djitmp(0,ij) + wz0*wt10*upara
           djitmp(1,ij) = djitmp(1,ij) + wz1*wt10*upara
           
           ij=jtion1(igyro,m)
           dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt01
           dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt01
           djitmp(0,ij) = djitmp(0,ij) + wz0*wt01*upara
           djitmp(1,ij) = djitmp(1,ij) + wz1*wt01*upara
           
           ij=ij+1
           dnitmp(0,ij) = dnitmp(0,ij) + wz0*wt11
           dnitmp(1,ij) = dnitmp(1,ij) + wz1*wt11
           djitmp(0,ij) = djitmp(0,ij) + wz0*wt11*upara
           djitmp(1,ij) = djitmp(1,ij) + wz1*wt11*upara
#else
! If no loop-level parallelism, use original algorithm (write directly
! into array "densityi()".
           ij=jtion0(igyro,m)
           densityi(0,ij) = densityi(0,ij) + wz0*wt00
           densityi(1,ij) = densityi(1,ij) + wz1*wt00
           currenti(0,ij) = currenti(0,ij) + wz0*wt00*upara
           currenti(1,ij) = currenti(1,ij) + wz1*wt00*upara
           
           ij=ij+1
           densityi(0,ij) = densityi(0,ij) + wz0*wt10
           densityi(1,ij) = densityi(1,ij) + wz1*wt10
           currenti(0,ij) = currenti(0,ij) + wz0*wt10*upara
           currenti(1,ij) = currenti(1,ij) + wz1*wt10*upara
           
           ij=jtion1(igyro,m)
           densityi(0,ij) = densityi(0,ij) + wz0*wt01
           densityi(1,ij) = densityi(1,ij) + wz1*wt01
           currenti(0,ij) = currenti(0,ij) + wz0*wt01*upara
           currenti(1,ij) = currenti(1,ij) + wz1*wt01*upara
           
           ij=ij+1
           densityi(0,ij) = densityi(0,ij) + wz0*wt11
           densityi(1,ij) = densityi(1,ij) + wz1*wt11
           currenti(0,ij) = currenti(0,ij) + wz0*wt11*upara
           currenti(1,ij) = currenti(1,ij) + wz1*wt11*upara
#endif
        enddo
     enddo
  endif

#ifdef _OPENMP
! For OpenMP, we now write the results accumulated in each thread-private
! array dnitmp() back into the shared array densityi(). The loop is enclosed
! in a critical section so that one thread at a time updates densityi(). 
!$omp critical
  do ij=1,mgrid
     densityi(0:1,ij) = densityi(0:1,ij) + dnitmp(0:1,ij)
     currenti(0:1,ij) = currenti(0:1,ij) + djitmp(0:1,ij)
  enddo
!$omp end critical
#endif
!$omp end parallel

! If we have a particle decomposition on the toroidal domains, do a reduce
! operation to add up all the contributions to charge density on the grid
  if(npartdom>1)then
!$omp parallel do private(ij)
     do ij=1,mgrid
        dnitmp(0:1,ij)=densityi(0:1,ij)
        densityi(0:1,ij)=0.
        djitmp(0:1,ij)=currenti(0:1,ij)
        currenti(0:1,ij)=0.
     enddo
     call MPI_ALLREDUCE(dnitmp,densityi,(mgrid*2),mpi_Rsize,MPI_SUM,partd_comm,ierror)
     call MPI_ALLREDUCE(djitmp,currenti,(mgrid*2),mpi_Rsize,MPI_SUM,partd_comm,ierror)
  endif
  
! poloidal end cell, discard ghost cell j=0
  do i=0,mpsi
     densityi(:,igrid(i)+mtheta(i))=densityi(:,igrid(i)+mtheta(i))+densityi(:,igrid(i))
     currenti(:,igrid(i)+mtheta(i))=currenti(:,igrid(i)+mtheta(i))+currenti(:,igrid(i))
  enddo

! toroidal end cell
  sendl(:,1)=densityi(0,:)
  sendl(:,2)=currenti(0,:)
  recvr=0.0
  icount=2*mgrid
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource

! send densityi to left and receive from right
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
       recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
     
  if(myrank_toroidal == mtoroidal-1)then
! B.C. at zeta=2*pi is shifted
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        densityi(1,ii+1:ii+jt)=densityi(1,ii+1:ii+jt)+cshift(recvr(ii+1:ii+jt,1),itran(i))
        currenti(1,ii+1:ii+jt)=currenti(1,ii+1:ii+jt)+cshift(recvr(ii+1:ii+jt,2),itran(i))
     enddo
  else
! B.C. at zeta<2*pi is continuous
     densityi(1,:)=densityi(1,:)+recvr(:,1)
     currenti(1,:)=currenti(1,:)+recvr(:,2)
  endif
  
! zero out charge in radial boundary cell
  do i=0,nbound-1
     densityi(:,igrid(i):igrid(i)+mtheta(i))=&
          densityi(:,igrid(i):igrid(i)+mtheta(i))*real(i)/real(nbound)
     densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=&
          densityi(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*real(i)/real(nbound)
     currenti(:,igrid(i):igrid(i)+mtheta(i))=&
          currenti(:,igrid(i):igrid(i)+mtheta(i))*real(i)/real(nbound)
     currenti(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))=&
          currenti(:,igrid(mpsi-i):igrid(mpsi-i)+mtheta(mpsi-i))*real(i)/real(nbound)
  enddo
  
! flux surface average and normalization  
  gyrodum=1.0/real(ngyroi)
  zonali=0.0
  zonic=0.0
!$omp parallel do private(i,j,ij)
  do i=0,mpsi
     do j=1,mtheta(i)
	ij=igrid(i)+j
        zonali(i)=zonali(i)+gyrodum*densityi(1,ij)
        densityi(1,ij)=gyrodum*densityi(1,ij)*markeri(ij)
        zonic(i)=zonic(i)+gyrodum*currenti(1,ij)
        currenti(1,ij)=gyrodum*currenti(1,ij)*markeri(ij)
     enddo
  enddo
  
! global sum of phi00, broadcast to every toroidal PE
  call MPI_ALLREDUCE(zonali,adum,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  zonali=adum*pmarki
  call MPI_ALLREDUCE(zonic,adum,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  zonic=adum*pmarki

! densityi subtracted (0,0) mode
!$omp parallel do private(i,j,ij)
  do i=0,mpsi
     do j=1,mtheta(i)
        ij=igrid(i)+j
        densityi(1,ij)=densityi(1,ij)-zonali(i)
        currenti(1,ij)=currenti(1,ij)-zonic(i)
     enddo
! poloidal BC condition
     densityi(1,igrid(i))=densityi(1,igrid(i)+mtheta(i))
     currenti(1,igrid(i))=currenti(1,igrid(i)+mtheta(i))
  enddo
  
! enforce charge conservation for zonal flow mode
  rdum=0.0
  tdum=0.0
  do i=1,mpsi-1
     r=a0+deltar*real(i)
     rdum=rdum+r
     tdum=tdum+r*zonali(i)
  enddo
  tdum=tdum/rdum
  do i=1,mpsi-1
     zonali(i)=zonali(i)-tdum
  enddo
  
end subroutine chargei
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! find ion location for interpolation in both gathering (chargei.F90) and scattering (pushi.F90)
subroutine locatei
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer m,igyro,ip,jt,ipjt,ii,im,j00,j01
  real(wp) psitmp,thetatmp,zetatmp,rhoi,r,wz1,rdum,wp1,tflr,tdum,&
       delr,delt(0:mpsi),delz,rhoratio,pi2_inv

  delr=1.0/deltar
  delt=2.0*pi/deltat
  delz=1.0/deltaz
  rhoratio=sqrt(aion)/(abs(qion)*rho0)
  if(ngyroi==1)rhoratio=0.0
  pi2_inv=0.5/pi

!$omp parallel do private(m,igyro,psitmp,thetatmp,zetatmp,rhoi,r,ip,jt,ipjt,&
!$omp& wz1,rdum,ii,wp1,tflr,im,tdum,j00,j01)
  do m=1,mi
     psitmp=zion(1,m)
     thetatmp=zion(2,m)
     zetatmp=zion(3,m)
     rhoi=zion(6,m)*rhoratio
     
     r=sqrt(2.0*psitmp)
     ip=max(0,min(mpsi,int((r-a0)*delr+0.5)))
     jt=max(0,min(mtheta(ip),int(thetatmp*pi2_inv*delt(ip)+0.5)))
     ipjt=igrid(ip)+jt
     wz1=(zetatmp-zeta0)*delz
     wzion(m)=wz1
     
     do igyro=1,ngyroi
        rdum=delr*max(0.0_wp,min(a1-a0,r+rhoi*pgyro(igyro,ipjt)-a0))
        ii=max(0,min(mpsi-1,int(rdum)))
        wp1=rdum-real(ii)
        wpion(igyro,m)=wp1

! particle position in theta
        tflr=thetatmp+rhoi*tgyro(igyro,ipjt)

! inner flux surface
        im=ii
        tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0
        tdum=(tdum-aint(tdum))*delt(im)
        j00=max(0,min(mtheta(im)-1,int(tdum)))
        jtion0(igyro,m)=igrid(im)+j00
        wtion0(igyro,m)=tdum-real(j00)

! outer flux surface
        im=ii+1
        tdum=pi2_inv*(tflr-zetatmp*qtinv(im))+10.0
        tdum=(tdum-aint(tdum))*delt(im)
        j01=max(0,min(mtheta(im)-1,int(tdum)))
        jtion1(igyro,m)=igrid(im)+j01
        wtion1(igyro,m)=tdum-real(j01)
     enddo
  enddo
  
end subroutine locatei

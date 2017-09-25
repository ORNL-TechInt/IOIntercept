subroutine diagnosis
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer,parameter :: nfield=3,mfdiag=4,mfdata1d=2
  integer i,j,k,ii,ip,icount,ierror,mtgrid,nf
  real(wp) fieldmode(2,modes,nfield),fieldtime(mfdiag,nfield),partdata(mpdiag,nspecies),&
       field00(0:mpsi,nfield),fieldrms(0:mpsi,nfield),rmstmp(0:mpsi,nfield)
  real(wp),dimension(:,:,:),allocatable :: fieldgrid

! open output files history.out and data1d.out
  if(mype==0 .and. istep==ndiag)call opendiag(mpdiag,nfield,mfdiag,mpdata1d,mfdata1d)
  call MPI_BCAST(mstepall,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
!!! particle diagnosis: main ion, electron and/or EP, imputities, ...
  partdata=0.0
  if(mype==0) then !---wj
  do ip=1,nspecies
     if(ip==1)then
        partdata(:,ip)=diagion/diagion(10)
     elseif(ip==2)then
        partdata(:,ip)=diagelectron/diagion(10)
     endif
  enddo
  endif !---wj

!!! field diagnosis: phi, a_para, fuild_ne, ...
  ii=igrid(iflux)
  mtgrid=mtheta(iflux)
  fieldtime=0.0 
  fieldmode=0.0
  allocate(fieldgrid(0:1,0:mtgrid,nfield),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
  
  do nf=1,nfield
     if(nf==1)then
        fieldgrid(:,:,nf)=phi(:,ii:ii+mtgrid)
!time history of field quantity at theta=zeta=0 & i=iflux
        fieldtime(1,nf)=phi(0,ii)/rho0 
        fieldtime(2,nf)=phip00(iflux)/rho0 
        fieldtime(3,nf)=sqrt(sum(phip00(0:mpsi)**2)/real(mpsi+1))/rho0
!$omp parallel do private(i)
        do i=0,mpsi
           field00(i,nf)=phip00(i)
           fieldrms(i,nf)=sum(phi(0,igrid(i):igrid(i)+mtheta(i)-1)**2)
        enddo

     elseif(nf==2)then
        fieldgrid(:,:,nf)=apara(:,ii:ii+mtgrid)
        fieldtime(1,nf)=apara(0,ii)
!$omp parallel do private(i)
        do i=0,mpsi
           field00(i,nf)=apara00(i)
           fieldrms(i,nf)=sum(apara(0,igrid(i):igrid(i)+mtheta(i)-1)**2)
        enddo

     elseif(nf==3)then
        fieldgrid(:,:,nf)=fluidne(:,ii:ii+mtgrid)
        fieldtime(1,nf)=fluidne(0,ii)
!$omp parallel do private(i)
        do i=0,mpsi
           field00(i,nf)=fluidne00(i)
           fieldrms(i,nf)=sum(fluidne(0,igrid(i):igrid(i)+mtheta(i)-1)**2)
        enddo
     endif
     if ( sum(fieldrms(:,nf)) .lt. 0. ) then
       print *, 'Error0: sum(fieldrms(:,nf)) ', sum(fieldrms(:,nf)), ' nf ', nf, ' nfield ', nfield, ' on pe ', mype
     endif
  enddo

!Volume-averaged RMS
  icount=nfield*(mpsi+1)
  rmstmp = 0. !---wj
  call MPI_REDUCE(fieldrms,rmstmp,icount,mpi_Rsize,MPI_SUM,0,toroidal_COMM,ierror)
  fieldrms=rmstmp
  do nf=1,nfield
     if ( sum(fieldrms(:,nf)) .lt. 0. ) then
       print *, 'Error: sum(fieldrms(:,nf)) ', sum(fieldrms(:,nf)), ' nf ', nf, ' nfield ', nfield, ' on pe ', mype
     endif
     if ( real(mtoroidal*sum(mtheta)) .le. 0. ) then
       print *, 'Error: real(mtoroidal*sum(mtheta)) ', real(mtoroidal*sum(mtheta)), ' on pe ', mype
     endif
     fieldtime(4,nf)=sqrt(sum(fieldrms(:,nf))/real(mtoroidal*sum(mtheta)))
  enddo

  call spectrum(nfield,mtgrid,fieldgrid,fieldmode)

  deallocate(fieldgrid)

  if(mype==0)then
! write particle and field data to history file
     do i=1,nspecies
        do j=1,mpdiag
           write(iodiag,102)partdata(j,i)
        enddo
     enddo
     do i=1,nfield
        do j=1,mfdiag
           write(iodiag,102)fieldtime(j,i)
        enddo
     enddo
     do i=1,nfield
        do j=1,modes
           write(iodiag,102)fieldmode(1,j,i),fieldmode(2,j,i)
        enddo
     enddo
!     call FLUSH(iodiag)

! write radial-time data to data1d.out
     write(iodata1d,102)data1di
     if(nspecies>1)write(iodata1d,102)data1de
     write(iodata1d,102)field00
     write(iodata1d,102)fieldrms
!     call FLUSH(iodata1d)

! monitor: rms of phi, apara, fluidne, heatflux of ion & electron
     write(gtcout,"(I5,6e15.6)")istep+mstepall,(fieldtime(4,i),i=1,nfield),&
          (partdata(9,i),i=1,nspecies)
!     call FLUSH(gtcout)
  endif

101 format(i6)
102 format(e12.6)

end subroutine diagnosis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! calculate (n,m) mode amplitude for diagnostics
subroutine spectrum(nfield,mtgrid,fieldgrid,fieldmode)
  use global_parameters
  use field_array
  implicit none

  integer nfield,mtgrid
  real(wp) fieldgrid(0:1,0:mtgrid,nfield),fieldmode(2,modes,nfield)
  integer i,j,k,jt,kz,mzeach,jpe,indp,indt,indp1,indt1,mteach,icount,ierror,nf
  real(wp) wt,dt,wz,zdum,tdum,pi2_inv,xz(mtdiag),fieldflux(mtdiag/mtoroidal,mtdiag,nfield),&
       eachzeta((mtdiag/mtoroidal)*(mtdiag/mtoroidal)*nfield),&
       allzeta((mtdiag/mtoroidal)*mtdiag*nfield)
  complex(wp) yzeta(nfield,mtdiag/mtoroidal*modes),ytheta(nfield,mtdiag*modes),&
       yz(mtdiag/2+1),ye(mtdiag) 

! mesh in magnetic coordinates
  mzeach=mtdiag/mtoroidal
  mteach=mzeach
  dt=2.0*pi/real(mtdiag)
  pi2_inv=0.5/pi

! Interpolate on a flux surface from fieldline coordinates to magnetic coordinates. 
! Use mtdiag for both poloidal and toroidal grid points.
!$omp parallel do private(kz,nf,j,zdum,tdum,jt,wt,wz)           
  do kz=1,mzeach
     wz=real(kz)/real(mzeach)
     zdum=zeta0+wz*deltaz
     do nf=1,nfield
        do j=1,mtdiag
           tdum=pi2_inv*(dt*real(j)-zdum*qtinv(iflux))+10.0
           tdum=(tdum-aint(tdum))*real(mtgrid)
           jt=max(0,min(mtgrid-1,int(tdum)))
           wt=tdum-real(jt)
           
           fieldflux(kz,j,nf)=wz*((1.0-wt)*fieldgrid(1,jt,nf)+wt*fieldgrid(1,jt+1,nf))+&
                (1.0-wz)*((1.0-wt)*fieldgrid(0,jt,nf)+wt*fieldgrid(0,jt+1,nf))
        enddo
     enddo
  enddo
 
! transpose 2-d matrix from (mtoroidal,mzeach) to (1,mtoroidal*mzeach)
  allzeta=0.0
  icount=mteach*mzeach*nfield  
  do jpe=0,mtoroidal-1

!$omp parallel do private(j,nf,k,jt,indt,indp1,indp)
     do j=1,mteach
        jt=jpe*mteach+j
        indt=(j-1)*mzeach
        do nf=1,nfield
           indp1=indt+(nf-1)*mteach*mzeach
           do k=1,mzeach
              indp=indp1+k
              eachzeta(indp)=fieldflux(k,jt,nf)
           enddo
        enddo
     enddo
     
     call MPI_GATHER(eachzeta,icount,mpi_Rsize,allzeta,icount,mpi_Rsize,jpe,toroidal_comm,ierror)
  enddo
    
! transform to k space
  yz=0.0
  do j=1,mteach
     indt1=(j-1)*mzeach
     do nf=1,nfield
        indt=indt1+(nf-1)*mteach*mzeach

!$omp parallel do private(kz,k,indp)      
        do kz=0,mtoroidal-1
           do k=1,mzeach
              indp=kz*icount+indt+k
              xz(kz*mzeach+k)=allzeta(indp)
           enddo
        enddo
      
        call fftr1d(1,mtdiag,1.0,xz,yz,2)

! record toroidal mode amplitude for diagnostic
        do kz=1,modes
           yzeta(nf,j+mteach*(kz-1))=yz(nmodes(kz)+1)
        enddo

     enddo
  enddo

! gather toroidal mode amplitude for calculation of poloidal mode amplitudes
  ytheta=0.0
  icount=mteach*modes*nfield
  call MPI_GATHER(yzeta,icount,mpi_Csize,ytheta,icount,mpi_Csize,0,toroidal_comm,ierror)

  icount=mteach*modes
  if(myrank_toroidal == 0)then
     do nf=1,nfield
        do kz=1,modes

!$omp parallel do private(i,j)     
           do i=0,mtoroidal-1
              do j=1,mteach
                 ye(j+i*mteach)=ytheta(nf,j+(kz-1)*mteach+i*icount)
              enddo
           enddo

           call fftc1d(1,mtdiag,1.0,ye)
           
           fieldmode(1,kz,nf)=real(ye(mtdiag-mmodes(kz)+1))/(real(mtdiag)*rho0)**2
           fieldmode(2,kz,nf)=aimag(ye(mtdiag-mmodes(kz)+1))/(real(mtdiag)*rho0)**2
        enddo
     enddo
  endif

end subroutine spectrum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! filter out some toroidal modes in linear simulation
subroutine filter(farray)

  use global_parameters
  use field_array
  implicit none
  
  integer i,j,k,ii,ij,ip,jt,kz,jpe,indp,indt,indp1,indt1,meachtheta,mteach,mzeach,jtp,&
       icount,ierror,idest,isource,isendtag,l,irecvtag,istatus(MPI_STATUS_SIZE)
  real(wp) wt,r,dt,wz,zdum,tdum,pi2_inv,ptemp,sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),&
       phism(mgrid),pleft(mthetamax),pright(mthetamax),phitmp(mgrid),ffilter(mtdiag/2+1),&
       xz(mtdiag),farray(0:1,mgrid),phiflux(mtdiag/mtoroidal,mtdiag,mpsi),&
       eachzeta((mtdiag/mtoroidal)*(mtdiag/mtoroidal)*mpsi),&
       allzeta((mtdiag/mtoroidal)*mtdiag*mpsi)
  complex(wp) yz(mtdiag/2+1)
     
! mesh in magnetic coordinates
  mzeach=mtdiag/mtoroidal
  mteach=mzeach
  dt=2.0*pi/real(mtdiag)
  pi2_inv=0.5/pi
  
  ffilter=0.0
  ffilter(nmodes+1)=1.0
  allzeta=0.0
! Interpolate on a flux surface from fieldline coordinates to magnetic
! coordinates. Use mtdiag for both poloidal and toroidal grid points.
!$omp parallel do private(kz,i,j,wz,ii,zdum,tdum,jt,wt)           
  do kz=1,mzeach
     wz=real(kz)/real(mzeach)
     zdum=zeta0+wz*deltaz
     do i=1,mpsi
        ii=igrid(i)
        do j=1,mtdiag
           tdum=pi2_inv*(dt*real(j)-zdum*qtinv(i))+10.0
           tdum=(tdum-aint(tdum))*real(mtheta(i))
           jt=max(0,min(mtheta(i)-1,int(tdum)))
           wt=tdum-real(jt)
           
           phiflux(kz,j,i)=((1.0-wt)*farray(1,ii+jt)+wt*farray(1,ii+jt+1))*&
                wz+(1.0-wz)*((1.0-wt)*farray(0,ii+jt)+wt*farray(0,ii+jt+1))
        enddo
     enddo
  enddo
  
! transpose 2-d matrix from (mtoroidal,mzeach) to (1,mtoroidal*mzeach)
  allzeta=0.0
  icount=mteach*mzeach*mpsi  
  do jpe=0,mtoroidal-1

!$omp parallel do private(j,i,k,jt,indt,indp1,indp)
     do j=1,mteach
        jt=jpe*mteach+j
        indt=(j-1)*mzeach
        do i=1,mpsi
           indp1=indt+(i-1)*mteach*mzeach
           do k=1,mzeach
              indp=indp1+k
              eachzeta(indp)=phiflux(k,jt,i)
           enddo
        enddo
     enddo
     
     call MPI_GATHER(eachzeta,icount,mpi_Rsize,allzeta,icount,&
          mpi_Rsize,jpe,toroidal_comm,ierror)
  enddo
        
! transform to k space
  yz=0.0
!$omp parallel do private(j,i,kz,k,indt1,indt,indp,xz,yz)        
  do j=1,mteach
     indt1=(j-1)*mzeach
     do i=1,mpsi
        indt=indt1+(i-1)*mteach*mzeach
        do kz=0,mtoroidal-1
           do k=1,mzeach
              indp=kz*icount+indt+k
              xz(kz*mzeach+k)=allzeta(indp)
           enddo
        enddo
        
        call fftr1d(1,mtdiag,1.0,xz,yz,2)
        
! linear run only keep a few modes
        yz=ffilter*yz

! transform back to real space
        call fftr1d(-1,mtdiag,1.0,xz,yz,2)

! transpose back to (mtoroidal,mzeach)
        do kz=0,mtoroidal-1
           do k=1,mzeach
              indp=kz*icount+indt+k
              allzeta(indp)=xz(kz*mzeach+k)
           enddo
        enddo
     enddo
  enddo

  do jpe=0,mtoroidal-1
     call MPI_SCATTER(allzeta,icount,mpi_Rsize,eachzeta,&
          icount,mpi_Rsize,jpe,toroidal_comm,ierror)

!$omp parallel do private(j,i,k,jt,indt,indp1,indp)
     do j=1,mteach
        jt=jpe*mteach+j
        indt=(j-1)*mzeach
        do i=1,mpsi
           indp1=indt+(i-1)*mteach*mzeach
           do k=1,mzeach
              indp=indp1+k
              phiflux(k,jt,i)=eachzeta(indp)
           enddo
        enddo
     enddo
  enddo
        
! interpolate field from magnetic coordinates to fieldline coordinates
!$omp parallel do private(i,j,ii,tdum,jt,wt,jtp)           
  do i=1,mpsi
     ii=igrid(i)              
     do j=1,mtheta(i)
        tdum=pi2_inv*(deltat(i)*real(j)+zeta1*qtinv(i))+10.0
        tdum=(tdum-aint(tdum))*real(mtdiag)
        jt=max(0,min(mtdiag-1,int(tdum)))
        wt=tdum-real(jt)
        jtp=jt+1
        if(jt==0)jt=mtdiag
           
        farray(1,ii+j)=wt*phiflux(mzeach,jtp,i)+(1.0-wt)*phiflux(mzeach,jt,i)
     enddo
  enddo

! toroidal BC: send phi to right and receive from left
  sendr=farray(1,:)
  recvl=0.0
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,recvl,&
       icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
  
  if(myrank_toroidal==0)then
!$omp parallel do private(i,ii,jt)
     do i=1,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        farray(0,ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
     enddo
  else
!$omp parallel do private(i)
     do i=1,mgrid
        farray(0,i)=recvl(i)
     enddo
  endif

! poloidal BC
!$omp parallel do private(i)    
  do i=1,mpsi
     farray(:,igrid(i))=farray(:,igrid(i)+mtheta(i))
  enddo

end subroutine filter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine opendiag(mpdiag,nfield,mfdiag,mpdata1d,mfdata1d)
  use global_parameters
  use field_array,only:modes
  implicit none

  integer mpdiag,nfield,mfdiag,mpdata1d,mfdata1d
  integer i,j,iotmp,ndstep,ndata
  real(wp) tdum
  character(len=100) cdum

  iodiag=123
  iotmp=456
  iodata1d=789

! # of time steps and # of data added to history file
  ndstep=mstep/ndiag
  ndata=(nspecies*mpdiag+nfield*(2*modes+mfdiag))
  if(irun==0)then
     mstepall=0
! open time history file in new run
     open(iodiag,file='history.out',status='replace')
     write(iodiag,101)ndstep,nspecies,mpdiag,nfield,modes,mfdiag
     write(iodiag,102)tstep*ndiag

! open output file for radial-time data in new run
     open(iodata1d,file='data1d.out',status='replace')            
     write(iodata1d,101)ndstep,mpsi+1,nspecies,mpdata1d,nfield,mfdata1d
  else

! find history file from previous run
     irest=irest-1
     if(mod(irest,2)==0)then
        cdum="restart_dir1/history_restart.out"
     else
        cdum="restart_dir2/history_restart.out"
     endif
     write(*,*)'read in ',trim(cdum)

     open(iotmp,file=trim(cdum),status='old')
     rewind(iotmp)

! open time history file for restart run
     open(iodiag,file='history.out',status='replace')

! # of time steps
     read(iotmp,101)mstepall
     write(iodiag,101)mstepall+ndstep
     do i=1,5
        read(iotmp,101)j
        write(iodiag,101)j
     enddo
     
! copy restart history data file
     do i=0,ndata*mstepall
        read(iotmp,102)tdum
        write(iodiag,102)tdum
     enddo
     close(iotmp)

!!copy restart data1d.out file
     open(iodata1d,file='data1d.out',status='replace')
     if(mod(irest,2)==0)then
        cdum="restart_dir1/data1d_restart.out"
     else
        cdum="restart_dir2/data1d_restart.out"
     endif
     write(*,*)'read in ',trim(cdum)

     irest=irest+1
     
     open(iotmp,file=trim(cdum),status='old')
     rewind(iotmp)
     read(iotmp,101)j
     write(iodata1d,101)mstepall+ndstep
     do i=1,5
        read(iotmp,101)j
        write(iodata1d,101)j
     enddo

     ndata=(mpsi+1)*(nspecies*mpdata1d+nfield*mfdata1d)
     do i=1,ndata*mstepall
        read(iotmp,102)tdum
        write(iodata1d,102)tdum
     enddo
     close(iotmp)

!!calculate current total time step        
     mstepall=mstepall*ndiag
     
  endif

101 format(i6)
102 format(e12.6)

end subroutine opendiag

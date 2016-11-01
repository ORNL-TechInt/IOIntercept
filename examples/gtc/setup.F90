Subroutine setup

  use global_parameters
  use particle_array
  use field_array
  use particle_tracking
  implicit none

  logical :: dir_exist
  integer i,j,k,ierror,ij,mtheta0,mtheta1,ip,jt,indp,indt,mtest,micell,mecell
  real(wp) r0,b0,etemperature,tdum,r,q,sint,dtheta_dx,rhoi,b,zdum,&
      edensity,delr,delt,rmax,rmin,wt
  namelist /key_parameters/ numberpe,mi,mgrid,mtheta0,mtheta1,mtdiag,delr,delt,&
      ulength,utime,rho0,betae,nparam
  integer*8 eightbyteint

! Read the input file that contains the run parameters
  call read_input_params(micell,mecell,r0,b0,etemperature,edensity)

! Set up the particle decomposition within each toroidal domain
  call set_particle_decomp

! read MHD equilibrium data
  CALL EQDATA

! first run plots equilibrium profile 
  if(mype==0)CALL EQPLOT

! equilibrium unit: length (unit=cm) and time (unit=second) unit
  ulength=r0
  utime=1.0_wp/(9580._wp*b0) ! time unit = inverse gyrofrequency of proton 
! rho0 rho_s=sqrt(T_e/m_i) in equilibrium unit, v_s=sqrt(T_e/m_i)
  rho0=102.0_wp*sqrt(aion*etemperature)/(abs(qion)*b0)/ulength
  tstep=tstep/rho0
!electron beta
  betae=4.03e-11_wp*edensity*etemperature/(b0*b0)
  
! basic ion-ion collision time, Braginskii definition
  if(tauii>0.0)then
     tauii=24.0_wp-log(sqrt(edensity)/etemperature)
     tauii=2.09e7_wp*(etemperature)**1.5_wp/(edensity*tauii*2.31_wp)*sqrt(2.0_wp)/utime
     tauii=0.532_wp*tauii
  endif

! allocate array for field quantites
  CALL fieldinitial

! gyro-averaging for sqrt(mu)=rho0 on grid of magnetic coordinates
  call gyroinitial

! initiate radial interpolation on poloidal grids for smooth.F90 and field.F90
  call rinterpolation

! number of particle
  eightbyteint = 1
  mi=micell*eightbyteint*(mgrid-mpsi)/npartdom          !# of ions per MPI process
  me=mecell*eightbyteint*(mgrid-mpsi)/npartdom          !# of electrons per MPI process

! # of species
  nspecies=1
  if(nhybrid>0)nspecies=2
! Particles is tracked by tagging them with a number with an extra element to the particle array
  nparam=6
  if(track_particles == 1)nparam=8

! write out key parameter
  if(mype == 0) then
     delr=deltar/rho0
     delt=deltat(mpsi/2)*(a0+deltar*real(mpsi/2))/rho0
     mtheta0=mtheta(0)
     mtheta1=mtheta(mpsi)
     write(gtcout,key_parameters)
     call FLUSH(gtcout)
  endif

end subroutine setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fieldinitial
  use global_parameters
  use field_array
  implicit none

  integer mtest,ierror,i,j,ii,ij,n,m,mtgrid
  real envelope,theta,r,q,tdum

! allocate memory
  allocate (qtinv(0:mpsi),deltat(0:mpsi),phi00(0:mpsi),phip00(0:mpsi),apara00(0:mpsi),fluidne00(0:mpsi),STAT=mtest,SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
  if (mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate qtinv: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  allocate (itran(0:mpsi),igrid(0:mpsi),mtheta(0:mpsi),STAT=mtest) !---wj
  if (mtest /= 0) then !---wj
     write(0,*)mype,'*** Cannot allocate qtinv: mtest=',mtest !---wj
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror) !---wj
  endif !---wj

  phi00=0.0
  phip00=0.0
  apara00=0.0
  fluidne00=0.0

! --- Define poloidal grid ---
! grid spacing
  deltar=(a1-a0)/real(mpsi)

! grid shift associated with fieldline following coordinates
  tdum=pi*aminor/real(mthetamax)
  do i=0,mpsi
     r=a0+deltar*real(i)
     mtheta(i)=2*max(1,int(pi*r/tdum+0.5_wp))
     deltat(i)=2.0_wp*pi/real(mtheta(i))
     q=q0+q1*r/aminor+q2*r*r/(aminor*aminor)
     itran(i)=int(real(mtheta(i))/q+0.5_wp)
     qtinv(i)=real(mtheta(i))/real(itran(i)) !q value for coordinate transformation
     qtinv(i)=1.0/qtinv(i) !inverse q to avoid divide operation
     itran(i)=itran(i)-mtheta(i)*(itran(i)/mtheta(i))
  enddo
! un-comment the next two lines to use magnetic coordinate
!  qtinv=0.0
!  itran=0
! total number of grids on a poloidal plane
  mgrid=sum(mtheta+1)

! When doing filtering and diagnostics of toroidal mode, we need to switch from the 
! field-line following coordinates alpha-zeta to the magnetic coordinate in theta-zeta. 
! This requires a greater number of grid points in the zeta direction, which
! is mtdiag. Precisely, mtdiag should be mtheta/q but since mtheta changes
! from one flux surface to another, we use a formula that gives us enough
! grid points for all the flux surfaces considered.
  mthetamax=mtheta(mpsi)
  mtdiag=(mthetamax/mtoroidal)*mtoroidal

! starting # in the label of poloidal grid on each flux surface
  igrid(0)=1
  do i=1,mpsi
     igrid(i)=igrid(i-1)+mtheta(i-1)+1
  enddo

! flux surface for mode diagnosis and filtering; default: 8 modes and on mpsi/2 flux surface
  modes=8
  iflux=mpsi/2
  mtgrid=mtheta(iflux)
  allocate(nmodes(modes),mmodes(modes))
!!! nfilter=0: default toroidal modes # k*delta_x=0.05,0.1,...,0.4
  if(nfilter==0)then
     do i=1,modes
        mmodes(i)=0.5+real(i*(mtheta(iflux)/2-1))/(20.0*pi)
        nmodes(i)=0.5+real(mmodes(i))*qtinv(iflux)
     enddo
  else
! nfilter=1: toroidal modes explicitly specified
     nmodes=(/15,15,15,15,15,15,15,15/)
     mmodes=(/17,18,19,20,21,22,23,24/)
  endif
  if(mype==0)write(gtcout,"(a8,8i4)")"nmodes=",nmodes,"mmodes=",mmodes

! allocate memory
  allocate(phi(0:1,mgrid),phieff(0:1,mgrid),phiind(0:1,mgrid),apara(0:1,mgrid),apara0(0:1,mgrid),&
       fluidne(0:1,mgrid),fluidne0(0:1,mgrid),fluidue(0:1,mgrid),gradue(3,0:1,mgrid),&
       gradphi(3,0:1,mgrid),packphi(4,0:1,mgrid),gradind(3,0:1,mgrid),gradapara(3,0:1,mgrid),STAT=mtest,SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
  if (mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate qtinv: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! inital adiabative electron density fluidne=cos(m*mtheta-n*zeta) with a radal envelope
! initial vector potential apara=0. All other dynamical quantities=0
  phi=0.0
  phieff=0.0
  phiind=0.0
  apara=0.0
  fluidne=0.0
  fluidue=0.0
  gradue=0.0
  gradphi=0.0
  gradind=0.0
  gradapara=0.0
  if(magnetic==1)then
     n=5
     m=8
     do i=0,mpsi
        envelope=1.0e-10
        envelope=envelope*sin(pi*real(i)/real(mpsi))
        do j=1,mtheta(i)
           ij=igrid(i)+j
           theta=deltat(i)*real(j)+zeta1*qtinv(i)
           fluidne(1,ij)=envelope*(cos(real(m)*theta-real(n)*zeta1)+&
                cos(real(m-1)*theta-real(n)*zeta1))
!                cos(real(m+1)*theta-real(n)*zeta1)+cos(real(m-2)*theta-real(n)*zeta1))
        enddo
     enddo
  endif

end subroutine fieldinitial

!=============================================================================
  Subroutine read_input_params(micell,mecell,r0,b0,etemperature,edensity)
!=============================================================================

  use global_parameters
  use particle_array
  use particle_tracking
  implicit none

  logical file_exist
  integer ierror,micell,mecell
  real(wp),intent(INOUT) :: r0,b0,etemperature,edensity

  namelist /input_parameters/ irun,mstep,msnap,ndiag,nhybrid,nonlinear,paranl,&
      mode00,tstep,micell,mecell,mpsi,mthetamax,mtoroidal,ncycle,aminor,a0,a1,q0,q1,q2,&
      rc,rw,aion,qion,denion,temion,aelectron,qelectron,kappati,kappate,kappan,&
      r0,b0,etemperature,edensity,ngyroi,nbound,iload,&
      tauii,track_particles,ndata3d,magnetic,numereq,nfilter
!
! Since it is preferable to have only one MPI process reading the input file,
! we choose the master process to set the default run parameters and to read
! the input file. The parameters will then be broadcast to the other processes.
!

  if(mype==0) then

! Test if the input file gtc.in exists
    inquire(file='gtc.in',exist=file_exist)
    if (file_exist) then
      open(55,file='gtc.in',status='old')
      read(55,nml=input_parameters)
      close(55)
    else
      write(*,*)'Cannot find file gtc.in !!!'
      stop
    endif

! Changing the units of a0 and a1 from units of "a" to units of "R_0"
    a0=a0*aminor
    a1=a1*aminor

    mstep=max(2,mstep)
    msnap=min(msnap,mstep/ndiag)
    if(nonlinear==0)then
       paranl=0.0_wp
       mode00=0
    else
       nfilter=0 !nonlinear run will all toroidal modes
    endif
    rc=rc*(a0+a1)
    rw=1.0_wp/(rw*(a1-a0))

! write parameters to standard output
    write(gtcout,input_parameters)
  endif

! Now send the parameter values to all the other MPI processes
  call broadcast_input_params(micell,mecell,r0,b0,etemperature,edensity)
  
end Subroutine read_input_params

!=============================================================================
  Subroutine broadcast_input_params(micell,mecell,r0,b0,etemperature,edensity)
!=============================================================================

  use global_parameters
  use particle_array
  use particle_tracking
  implicit none

  integer,parameter :: n_integers=22,n_reals=24
  integer  :: integer_params(n_integers)
  real(wp) :: real_params(n_reals)
  integer ierror,micell,mecell
  real(wp),intent(INOUT) :: r0,b0,etemperature,edensity

! The master process, mype=0, holds all the input parameters. We need
! to broadcast their values to the other processes. Instead of issuing
! an expensive MPI_BCAST() for each parameter, it is better to pack
! everything in a single vector, broadcast it, and unpack it.

  if(mype==0)then
!   Pack all the integer parameters in integer_params() array
    integer_params(1)=irun
    integer_params(2)=mstep
    integer_params(3)=msnap
    integer_params(4)=ndiag
    integer_params(5)=nhybrid
    integer_params(6)=mode00
    integer_params(7)=micell
    integer_params(8)=mecell
    integer_params(9)=mpsi
    integer_params(10)=mthetamax
    integer_params(11)=mtoroidal
    integer_params(12)=npartdom
    integer_params(13)=ncycle
    integer_params(14)=ngyroi
    integer_params(15)=nbound
    integer_params(16)=iload
    integer_params(17)=track_particles
    integer_params(18)=ndata3d
    integer_params(19)=magnetic
    integer_params(20)=nonlinear
    integer_params(21)=numereq
    integer_params(22)=nfilter

!   Pack all the real parameters in real_params() array
    real_params(1)=paranl
    real_params(2)=tstep
    real_params(3)=aminor
    real_params(4)=a0
    real_params(5)=a1
    real_params(6)=q0
    real_params(7)=q1
    real_params(8)=q2
    real_params(9)=rc
    real_params(10)=rw
    real_params(11)=aion
    real_params(12)=qion
    real_params(13)=aelectron
    real_params(14)=qelectron
    real_params(15)=kappati
    real_params(16)=kappate
    real_params(17)=kappan
    real_params(18)=r0
    real_params(19)=b0
    real_params(20)=etemperature
    real_params(21)=edensity
    real_params(22)=tauii
    real_params(23)=denion
    real_params(24)=temion
  endif

! Send input parameters to all processes
  call MPI_BCAST(integer_params,n_integers,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(real_params,n_reals,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

  if(mype/=0)then
!   Unpack integer parameters
    irun=integer_params(1)
    mstep=integer_params(2)
    msnap=integer_params(3)
    ndiag=integer_params(4)
    nhybrid=integer_params(5)
    mode00=integer_params(6)
    micell=integer_params(7)
    mecell=integer_params(8)
    mpsi=integer_params(9)
    mthetamax=integer_params(10)
    mtoroidal=integer_params(11)
    npartdom=integer_params(12)
    ncycle=integer_params(13)
    ngyroi=integer_params(14)
    nbound=integer_params(15)
    iload=integer_params(16)
    track_particles=integer_params(17)
    ndata3d=integer_params(18)
    magnetic=integer_params(19)
    nonlinear=integer_params(20)
    numereq=integer_params(21)
    nfilter=integer_params(22)

!   Unpack real parameters
    paranl=real_params(1)
    tstep=real_params(2)
    aminor=real_params(3)
    a0=real_params(4)
    a1=real_params(5)
    q0=real_params(6)
    q1=real_params(7)
    q2=real_params(8)
    rc=real_params(9)
    rw=real_params(10)
    aion=real_params(11)
    qion=real_params(12)
    aelectron=real_params(13)
    qelectron=real_params(14)
    kappati=real_params(15)
    kappate=real_params(16)
    kappan=real_params(17)
    r0=real_params(18)
    b0=real_params(19)
    etemperature=real_params(20)
    edensity=real_params(21)
    tauii=real_params(22)
    denion=real_params(23)
    temion=real_params(24)
  endif

end subroutine broadcast_input_params

!=============================================================================
    Subroutine set_particle_decomp
!=============================================================================

  use global_parameters
  use field_array
  implicit none

  integer  :: i,j,k,pe_number,mtest,ierror

! ----- First we verify the consistency of mtoroidal and npartdom -------
! The number of toroidal domains (mtoroidal) times the number of particle
! "domains" (npartdom) needs to be equal to the number of processor "numberpe".
! numberpe cannot be changed since it is given on the command line.

! numberpe must be a multiple of mtoroidal
  if(mod(numberpe,mtoroidal) /= 0)then
    write(gtcout,*)'Wrong PE #; PE=',numberpe,'Mtoroidal=',mtoroidal
    stop
  endif

!number of particle decomposition.
  npartdom=numberpe/mtoroidal
  if(mype==0)then
    write(gtcout,*)'*******************************************************'
    write(gtcout,*)'  Using npartdom=',npartdom,' and mtoroidal=',mtoroidal
    write(gtcout,*)'*******************************************************'
    write(gtcout,*)
  endif

! Make sure that "mpsi", the total number of flux surfaces, is an even
! number since this quantity will be used in Fast Fourier Transforms
  mpsi=2*(mpsi/2)

! We now give each PE (task) a unique domain identified by 2 numbers: the
! particle and toroidal domain numbers.
!    particle_domain_location = rank of the particle domain holding mype
!    toroidal_domain_location = rank of the toroidal domain holding mype
! 
! On the IBM SP, the MPI tasks are distributed in an orderly fashion to each
! node unless the LoadLeveler instruction "#@ blocking = unlimited" is used.
! On Seaborg for example, the first 16 tasks (mype=0-15) will be assigned to
! the first node that has been allocated to the job, then the next 16
! (mype=16-31) will be assigned to the second node, etc. When not using the
! OpenMP, we want the particle domains to sit on the same node because
! communication is more intensive. To achieve this, successive PE numbers are
! assigned to the particle domains first.
! It is easy to achieve this ordering by simply using mype/npartdom for
! the toroidal domain and mod(mype,npartdom) for the particle domain.
!
!  pe_number=0
!  do j=0,mtoroidal-1
!     do i=0,npartdom-1
!        pe_grid(i,j)=pe_number
!        particle_domain_location(pe_number)=i
!        toroidal_domain_location(pe_number)=j
!        pe_number=pe_number+1
!     enddo
!  enddo

  !---wj - Changed for improved communication performance
  !particle_domain_location=mod(mype,npartdom) !---wj
  !toroidal_domain_location=mype/npartdom !---wj
  particle_domain_location=mype/mtoroidal !---wj
  toroidal_domain_location=mod(mype,mtoroidal) !---wj

!  write(0,*)'mype=',mype,"  particle_domain_location =",&
!            particle_domain_location,' toroidal_domain_location =',&
!            toroidal_domain_location,' pi=',pi

! Domain decomposition in toroidal direction.
  zeta0=2.0_wp*pi*real(toroidal_domain_location)/real(mtoroidal)
  zeta1=2.0_wp*pi*real(toroidal_domain_location+1)/real(mtoroidal)

!  write(0,*)mype,' in set_particle_decomp:  zeta0=',&
!            zeta0,'  zeta1=',zeta1

! grid spacing in the toroidal direction, one toroidal grip per toroidal domain
  deltaz=zeta1-zeta0

! ---- Create particle domain communicator and toroidal communicator -----
! We now need to create a new communicator which will include only the
! processes located in the same toroidal domain. The particles inside
! each toroidal domain are divided equally between "npartdom" processes.
! Each one of these processes will do a charge deposition on a copy of
! the same grid, requiring a toroidal-domain-wide reduction after the
! deposition. The new communicator will allow the reduction to be done
! only between those processes located in the same toroidal domain.
!
! We also need to create a purely toroidal communicator so that the
! particles with the same particle domain id can exchange with their
! toroidal neighbors.
!
! Form 2 subcommunicators: one that includes all the processes located in
! the same toroidal domain (partd_comm), and one that includes all the
! processes part of the same particle domain (toroidal_comm).
! Here is how to create a new communicator from an old one by using
! the MPI call "MPI_COMM_SPLIT()".
! All the processes passing the same value of "color" will be placed in
! the same communicator. The "rank_in_new_comm" value will be used to
! set the rank of that process on the communicator.
!  call MPI_COMM_SPLIT(old_comm,color,rank_in_new_comm,new_comm,ierror)

! particle domain communicator (for communications between the particle
! domains WITHIN the same toroidal domain)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,toroidal_domain_location,&
                      particle_domain_location,partd_comm,ierror)

! toroidal communicator (for communications BETWEEN toroidal domains of same
! particle domain number)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,particle_domain_location,&
                      toroidal_domain_location,toroidal_comm,ierror)

  call mpi_comm_size(partd_comm,nproc_partd,ierror)
  call mpi_comm_rank(partd_comm,myrank_partd,ierror)

  call mpi_comm_size(toroidal_comm,nproc_toroidal,ierror)
  call mpi_comm_rank(toroidal_comm,myrank_toroidal,ierror)

!  write(0,*)'mype=',mype,'  nproc_toroidal=',nproc_toroidal,&
!       ' myrank_toroidal=',myrank_toroidal,'  nproc_partd=',nproc_partd,&
!       ' myrank_partd=',myrank_partd

  if(nproc_partd/=npartdom)then
    write(0,*)'*** nproc_partd=',nproc_partd,' NOT EQUAL to npartdom=',npartdom
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  if(nproc_toroidal/=mtoroidal)then
    write(*,*)'*** nproc_toroidal=',nproc_toroidal,' NOT EQUAL to mtoroidal=',&
              mtoroidal
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! We now find the toroidal neighbors of the current toroidal domain and
! store that information in 2 easily accessible variables. This information
! is needed several times inside the code, such as when particles cross
! the domain boundaries. We will use the toroidal communicator to do these
! transfers so we don't need to worry about the value of myrank_partd.
! We have periodic boundary conditions in the toroidal direction so the
! neighbor to the left of myrank_toroidal=0 is (mtoroidal-1).

  left_pe=mod(myrank_toroidal-1+mtoroidal,mtoroidal)
  right_pe=mod(myrank_toroidal+1,mtoroidal)

end subroutine set_particle_decomp

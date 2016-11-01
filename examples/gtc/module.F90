module precision
  include 'mpif.h'
  !use mpi
  integer, parameter :: doubleprec=selected_real_kind(12),&
    singleprec=selected_real_kind(6),&
    defaultprec=kind(0.0)
#ifdef DOUBLE_PRECISION
  integer, parameter :: wp=doubleprec,mpi_Rsize=MPI_DOUBLE_PRECISION,&
                        mpi_Csize=MPI_DOUBLE_COMPLEX
#else
  integer, parameter :: wp=singleprec,mpi_Rsize=MPI_REAL,&
                        mpi_Csize=MPI_COMPLEX
#endif
  save
end module precision

module global_parameters
  use precision
  integer,parameter :: gtcout=11

! control parameters
  integer track_particles,ndata3d,mgrid,mpsi,mthetamax,mtoroidal,iodiag,iodata1d,&
       nspecies,istep,ndiag,msnap,mstep,mstepall,mode00,nbound,irun,iload,irk,idiag,&
       ncycle,mtdiag,nhybrid,ihybrid,nparam,magnetic,nonlinear,numereq,nfilter
  real(wp) paranl,a0,a1,aminor,q0,q1,q2,pi,tstep,ulength,utime,rho0,rc,rw

! MPI toroidal and particle decomposion
  integer :: mype,numberpe,npartdom,toroidal_comm,partd_comm,nproc_partd,myrank_partd,&
       nproc_toroidal,myrank_toroidal,left_pe,right_pe,&
       toroidal_domain_location,particle_domain_location
#if RESTART_TIME 
  integer,allocatable,dimension(:) :: data_size
#endif

!!XY rolling restart
  integer :: irest,FileExit
  character(len=10):: restart_dir1,restart_dir2
  save
end module global_parameters

module equilibrium
  integer lsp,lst
  real psiw,spdpsi,spdtheta,spdrg,spdtor
  real,dimension(:),allocatable :: stpp,mipp,mapp
  real,dimension(:,:),allocatable :: qpsi,gpsi,ppsi,rpsi,torpsi,tpsi,npsi,nepp,tepp,&
       tipp,zepp,ropp,erpp,spcos,spsin,rgpsi,psirg,psitor,cpsi
  real,dimension(:,:,:),allocatable :: b,x,z,g,rd,nu,dl,ha,hb
  save
end module equilibrium

module particle_array
  use precision
! particle diagnostics: # of quantities per species in history.out and data1d.out  
  integer,parameter :: mpdiag=10,mpdata1d=2

! electron
  integer me,me1,memax
  real(wp) qelectron,aelectron,kappate,kappan,betae
  real(wp),dimension(mpdiag) :: diagelectron
  real(wp),dimension(:),allocatable :: rteme,pmarke,zonale,markere,rdteme,pfluxe
  real(wp),dimension(:,:),allocatable :: zelectron,zelectron0,zelectron1,densitye,phit,data1de
  real(wp),dimension(:,:,:),allocatable :: phisave

! ion
  integer mi,mimax,ngyroi
  real(wp) qion,aion,denion,temion,kappati,tauii
  real(wp),dimension(mpdiag) :: diagion
  integer,dimension(:,:),allocatable :: jtion0,jtion1
  real(wp),dimension(:),allocatable :: wzion,rtemi,rden,pmarki,zonali,zonic,markeri,rdtemi,pfluxi
  real(wp),dimension(:,:),allocatable :: zion,zion0,wpion,wtion0,wtion1,densityi,currenti,data1di
  save
end module particle_array

module particle_tracking
  use precision
  real(wp),dimension(:,:),allocatable :: ptrackedi,ptrackede
  integer,dimension(2) :: ntrackp
  integer,dimension(:),allocatable :: np_all,offset
  integer*8 :: my_offset,np_total
  save
end module particle_tracking

module field_array
  use precision
! PIC global fieldline following mesh
  real(wp) deltar,deltaz,zeta1,zeta0
  integer,dimension(:),allocatable :: itran,igrid,mtheta
  real(wp),dimension(:),allocatable :: deltat,qtinv

! fields on mesh: phi, phieff, phind, apara, fluidne, fluidue, zonal and gradient
  real(wp),dimension(:),allocatable :: phi00,phip00,apara00,fluidne00
  real(wp),dimension(:,:),allocatable :: phi,phieff,phiind,apara,fluidne,fluidue,apara0,fluidne0
  real(wp),dimension(:,:,:),allocatable :: gradphi,gradind,gradapara,gradue,packphi

! diagnostics and filtering
  integer iflux,modes
  integer,dimension(:),allocatable :: nmodes,mmodes

! radial interpolation
  integer,dimension(:,:),allocatable :: jtp1,jtp2
  real(wp),dimension(:,:),allocatable :: wtp1,wtp2

! gyro averaging
  real(wp),dimension(:,:),allocatable :: pgyro,tgyro
  save
end module field_array

module petsc_array
  use precision
  integer :: newcomm,nproc_newcomm,myrank_newcomm
  integer,dimension(:),allocatable :: userj,users
  real(wp),dimension(:),allocatable :: usera,userb,userx
  save
end module petsc_array

module data_type
  integer, parameter :: bp_char=0,bp_short=1,&
                bp_int=2,bp_long=3,bp_longlong=4,bp_float=5,&
                bp_double=6,bp_longdouble=7,bp_pointer=8,&
                bp_string=9,bp_uchar=50,bp_ushort=51,bp_uint=52,&
                bp_ulong=53,bp_ulonglong=54
  save
end module data_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!                 Gyrokinetic Toroidal Code (GTC)                            !
!                          Version 2, 2008                                   !
!                 University of California, Irvine                           !
!                                                                            !
! Current developers:                                                        !
! UCI: Z. Lin (zhihongl@uci.edu), I. Holod, Y. Xiao, W. L. Zhang             ! 
! ORNL: S. Klasky                                                            ! 
! PPPL: S. Ethier                                                            ! 
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program GTC
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer i
  real(doubleprec) timecpu(9),t1,tchargei,tchargee
  character(len=12) cdate(4)
  tchargei=0.
  tchargee=0.

! Initialize MPI, OpenMP, Adios, PETSc, Random number generator etc
  CALL INITIAL(cdate,timecpu)

! input parameters, equilibrium & profile data, coefficients for gyroaveraging & interpolation
  CALL SETUP

! initialize particle position and velocity for all species
  CALL LOAD

! calculate ion gather-scatter coefficients for chargei and pushi
  CALL LOCATEI

! main time loop
  do istep=1,mstep
     do irk=1,2

! idiag=0: do time history diagnosis at irk=1 & istep=ndiag*N
        idiag=mod(irk+1,2)+mod(istep,ndiag)
        if(idiag==0)then
           if(ndata3d==1)call DATAOUT3D !write 3D fluid data
           if(track_particles==1)CALL PARTOUT !write particle data
           call timer(timecpu(9),timecpu(7))
        endif

! gradients of phi, phiind, fluidue, apara; and zonal flow
        CALL FIELD_GRADIENT
        call timer(timecpu(9),timecpu(6))

! push ion
        
        call packup_phi(1)
        CALL PUSHI
! push fields of apara, fluidne
        if(magnetic==1)CALL PUSHFIELD
        call timer(timecpu(9),timecpu(1))
        
! redistribute ion across PEs
        CALL SHIFTI(mimax,mi)
        call timer(timecpu(9),timecpu(2))

! ion perturbed density
        CALL LOCATEI
        t1=mpi_wtime()
        CALL CHARGEI
        tchargei=tchargei+(mpi_wtime()-t1)
        call timer(timecpu(9),timecpu(3))

! smooth ion density & current
        CALL SMOOTH(densityi)
        if(magnetic==1)CALL SMOOTH(currenti)
        call timer(timecpu(9),timecpu(5))

        if(magnetic==0)then
! solve GK Poisson equation for phi using adiabatic electron
           CALL POISSON_SOLVER("adiabatic-electron")
        else
! solver fields of phi, phiind, fluidue in EM simulation
           CALL FIELD_SOLVER
        endif
           call timer(timecpu(9),timecpu(4))

! smooth potential
        CALL SMOOTH(phi)
        if(magnetic==1)then
           CALL SMOOTH(phiind)
           CALL SMOOTH(fluidue)
        endif
        call timer(timecpu(9),timecpu(5))
! use one or a few toroidal modes in linear simulation
        if(nfilter==1)CALL FILTER(phi)
        call timer(timecpu(9),timecpu(4))

        call packup_phi(1)
        do ihybrid=1,nhybrid

! time derivative of effective potential
           CALL DPHIEFF
           call timer(timecpu(9),timecpu(5))
           call packup_phi(2)

! push electron, sub-cycling
           do i=1,ncycle*irk
! 1st RK step
              CALL PUSHE(i,1)
              call timer(timecpu(9),timecpu(1))
             
              CALL SHIFTE(memax,me)
              call timer(timecpu(9),timecpu(2))
             
! 2nd RK step
              CALL PUSHE(i,2)
              call timer(timecpu(9),timecpu(1))

              CALL SHIFTE(memax,me)
              call timer(timecpu(9),timecpu(2))
           enddo
           
! nonadiabatic electron charge density
           t1=mpi_wtime()
           CALL CHARGEE
           tchargee=tchargee+(mpi_wtime()-t1)
           call timer(timecpu(9),timecpu(3))
           
! smooth electron density
           CALL SMOOTH(densitye)
           call timer(timecpu(9),timecpu(5))

! solve GK Poisson equation using non-adiabatic electron
           CALL POISSON_SOLVER("kinetic-electron")
           call timer(timecpu(9),timecpu(4))

! smooth potential
           CALL SMOOTH(phi)
           call timer(timecpu(9),timecpu(5))
! use one or a few toroidal modes in linear simulation; filter at istep=ndiag
           if(nfilter==1)CALL FILTER(phi)
           call timer(timecpu(9),timecpu(4))
        enddo

!write diagnostics and 1D data
        if(idiag==0)CALL DIAGNOSIS
     enddo

! profile snapshots, write particle information to restart file
     if(mod(istep,mstep/msnap)==0)then
!        CALL SNAPSHOT
!        CALL RESTART_IO("write")
        call timer(timecpu(9),timecpu(7))
     endif
  enddo

!write(6,*)'Chargei time=',tchargei,'Chargee time=',tchargee
! Program finalize
 CALL FINAL(cdate,timecpu)

end program GTC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initial(cdate,timecpu)
  use global_parameters
  implicit none

  integer ierror,i,m,nsize
  real(doubleprec) timecpu(9)
  integer,dimension(:),allocatable :: mget,mput
  character(len=12) cdate(4)
#ifdef _OPENMP
  integer nthreads,omp_get_num_threads
#endif

! MPI initialize, total # of PE, and rank of PE
  call mpi_init(ierror)
  call mpi_comm_size(mpi_comm_world,numberpe,ierror)
  call mpi_comm_rank(mpi_comm_world,mype,ierror)

! initialize timer
  istep=0
  timecpu=0.0
  timecpu(8)=mpi_wtime()
  timecpu(9)=timecpu(8)

! open standard output file, record program starting time
  if(mype==0)then
     open(gtcout,file='gtc.out',status='replace')
     call date_and_time(cdate(1),cdate(2))

! check OPENMP
#ifdef _OPENMP
!$omp parallel private(nthreads)
     nthreads=omp_get_num_threads()  !Get the number of threads if using OMP
!$omp single
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Number of OpenMP threads = ',nthreads
     write(gtcout,'("===================================",/)')
!$omp end single nowait
!$omp end parallel
#else
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Run without OpenMP threads'
     write(gtcout,'("===================================",/)')
#endif
  endif

! initiate adios
#if ADIOS
  !call adios_init ("config.xml"//char(0), ierror)! PLEASE KEEP THIS CHANGE!
#endif

! PETSc initialize
#ifdef __PETSc
  call petsc_init
#endif

!!XY rolling restart, initialize restart directory and control var
  irest=0
  FileExit=0

! numerical constant
  pi=4.0_wp*atan(1.0_wp)

! **** Use the intrinsic F90 random number generator *****
! initialize f90 random number generator
!    call random_seed !comemt out this line for identitical random number 
  call random_seed(size=nsize)
  allocate(mget(nsize),mput(nsize))
  call random_seed(get=mget)
  do i=1,nsize 
     call system_clock(m)   !random initialization for collision
     if(irun==0)m=1         !same initialization
     mput(i)=111111*(mype+1)+m+mget(i)
  enddo
  call random_seed(put=mput)
  deallocate(mget,mput)
  
end subroutine initial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine final(cdate,timecpu)
  use global_parameters
  implicit none

  integer ierror,iorestart
  real(doubleprec) timecpu(9)
  character(len=12) cdate(4)
  character(len=8) ic(8)
!total cpu and wall clock time
  timecpu(9)=mpi_wtime()
  timecpu(8)=timecpu(9)-timecpu(8)
  ic(1)='pusher'
  ic(2)='shift'
  ic(3)='charge'
  ic(4)='poisson'
  ic(5)='smooth'
  ic(6)='field'
  ic(7)='diag'
  ic(8)='total'
  if(mype==0)then
     write(gtcout,*)'CPU TIME USAGE (in SEC):'
     write(gtcout,*)ic
     write(gtcout,'(8(1pe8.1),/)')timecpu(1:8)

! Restart file info
     FileExit=1
     iorestart=345
     open(iorestart,file="FileExit.out",status="replace")
     write(iorestart,"(A9,i1)")"FileExit=",FileExit
     write(iorestart,"(A9,i5)")"irest   =",irest
     if(mod(irest+1,2)==0)then
        write(iorestart,"(A12)")"restart_dir1"
     else
        write(iorestart,"(A12)")"restart_dir2"
     endif
     close(iorestart)
     
! record program end time
     call date_and_time(cdate(3),cdate(4))
     write(gtcout,*)'Program starts at DATE=', cdate(1), 'TIME=', cdate(2)
     write(gtcout,*)'Program ends at   DATE=', cdate(3), 'TIME=', cdate(4)

! close standard output file
     if(gtcout /= 6 .and. gtcout /= 0)close(gtcout)
  endif

! PETSc finalize
#ifdef __PETSc
  call petsc_final
#endif

#if ADIOS
  !call adios_finalize (mype,ierror)
#endif

! MPI finalize
  call mpi_finalize(ierror)

end subroutine final
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=========================================
subroutine timer(t0,timecpu)
!=========================================
  use precision
  use global_parameters
  implicit none
  real(doubleprec) t0,t1,dt,timecpu

! Get cpu usage time since the begine of run and subtract value from the previous call
!  call cpu_time(t1)
  t1=mpi_wtime()
  dt=t1-t0
  t0=t1
  timecpu=timecpu+dt

! Get wall clock time and subtract value from the previous call
!  t1wc=MPI_WTIME()
!  dtwc=t1wc-t0wc
!  t0wc=t1wc
!  timewc=timewc+dtwc

end subroutine timer

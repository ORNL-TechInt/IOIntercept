subroutine shifti(mpmax,mp)
! this subroutine is identical to shifte,except for changing "electron" to "ion"
  use global_parameters
  use particle_array
  use field_array
  implicit none
  
  integer mpmax,mp,i,m,msendleft(2),msendright(2),mrecvleft(2),mrecvright(2),mtop,m0,msend,&
       msendtmp,mrecvtmp,idest,isource,isendtag,irecvtag,nzphase,kzparticle(mpmax),iright(mpmax),&
       ileft(mpmax),isendcount,irecvcount,istatus(MPI_STATUS_SIZE),ierror,iteration,lasth
  real(wp),dimension(:,:),allocatable :: recvleft,recvright,sendleft,sendright
  real(wp) zetaright,zetaleft,pi_inv
#ifdef _OPENMP
  integer msleft(32,0:15),msright(32,0:15)
  integer nthreads,gnthreads,iam,delm,mbeg,mend,omp_get_num_threads,omp_get_thread_num
#endif

  nzphase=2*nparam ! nzion=14 if track_particles=1, =12 otherwise
  pi_inv=1.0/pi
  m0=1
  iteration=0
  
100 iteration=iteration+1
!  if(iteration>mtoroidal)then
!     write(*,*)'endless particle sorting loop at PE=',mype
!     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
!  endif
! change - wj
  if(iteration>mtoroidal.and.m0<=mp)then
     write(*,*)'shifti endless particle sorting loop at PE=',mype
     write(*,*)'At PE=',mype," m0 and mp equal ",m0,mp,"zelectron(3,m0)=",zelectron(3,m0)
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  msend=0
  msendright=0
  msendleft=0

  if(m0 <= mp)then
!$omp parallel do private(m)
     do m=1,mp
        kzparticle(m)=0
     enddo

#ifdef _OPENMP
! This section of code (down to #else) is included by the preprocessor when
! the compilation is performed with OpenMP support. We must then use a few
! temporary arrays and add some work distribution code.
! First, we initialize the shared (and temporary) arrays msleft and msright
! to zero.
     msleft=0
     msright=0

! Then we start the parallel region with
!$omp parallel private(nthreads,iam,delm,i,mbeg,mend,m,zetaright,zetaleft) &
!$omp& shared(gnthreads,msleft,msright)
     nthreads=omp_get_num_threads()    !Get the number of threads ready to work
     iam=omp_get_thread_num()          !Get my thread number (position)
     delm=(mp-m0+1)/nthreads       !Calculate the number of steps per thread
     i=mod((mp-m0+1),nthreads)
!$omp single                 !Put nthread in global variable for later use.
     gnthreads=nthreads      !nthread is the same for all threads so only one
!$omp end single nowait      !of them needs to copy the value in gnthreads

! We now distribute the work between the threads. The loop over the particles
! is distributed equally (as much as possible) between them.
     mbeg=m0+min(iam,i)*(delm+1)+max(0,(iam-i))*delm
     mend=mbeg+delm+(min((iam+1),i)/(iam+1))-1

! label particle to be moved
     do m=mbeg,mend
        zetaright=min(2.0*pi,zion(3,m))-zeta1
        zetaleft=zion(3,m)-zeta0
        
        if(zetaright*zetaleft > 0)then
           zetaright=zetaright*0.5*pi_inv
           zetaright=zetaright-real(floor(zetaright))
           msright(3,iam)=msright(3,iam)+1
           kzparticle(mbeg+msright(3,iam)-1)=m
           
           if(zetaright < 0.5)then
! particle to move right
              msright(1,iam)=msright(1,iam)+1
              iright(mbeg+msright(1,iam)-1)=m

! particle to move left
           else
              msleft(1,iam)=msleft(1,iam)+1
              ileft(mbeg+msleft(1,iam)-1)=m
           endif
        endif
     enddo
! End of the OpenMP parallel region
!$omp end parallel

! Now that we are out of the parallel region we need to gather and rearrange
! the results of the multi-thread calculation. We need to end up with the
! same arrays as for the sequential (single-threaded) calculation.
     do m=0,gnthreads-1
        delm=(mp-m0+1)/gnthreads
        i=mod((mp-m0+1),gnthreads)
        mbeg=m0+min(m,i)*(delm+1)+max(0,(m-i))*delm
        if( msleft(2,m) /= 0 )msendleft(2)=msendleft(1)+msleft(2,m)
        do i=1,msleft(1,m)
           ileft(msendleft(1)+i)=ileft(mbeg+i-1)
        enddo
        msendleft(1)=msendleft(1)+msleft(1,m)
        if( msright(2,m) /= 0 )msendright(2)=msendright(1)+msright(2,m)
        do i=1,msright(1,m)
           iright(msendright(1)+i)=iright(mbeg+i-1)
        enddo
        msendright(1)=msendright(1)+msright(1,m)
        do i=1,msright(3,m)
           kzparticle(msend+i)=kzparticle(mbeg+i-1)
        enddo
        msend=msend+msright(3,m)
     enddo

#else
!  This section of code replaces the section above when the compilation does
!  NOT include the OpenMP support option. Temporary arrays msleft and msright
!  are not needed as well as the extra code for thread work distribution.

     do m=m0,mp
        zetaright=min(2.0*pi,zion(3,m))-zeta1
        zetaleft=zion(3,m)-zeta0
        
        if( zetaright*zetaleft > 0 )then
           zetaright=zetaright*0.5*pi_inv
           zetaright=zetaright-real(floor(zetaright))
           msend=msend+1
           kzparticle(msend)=m
           
           if( zetaright < 0.5 )then
! # of particle to move right
              msendright(1)=msendright(1)+1
              iright(msendright(1))=m

! # of particle to move left
           else
              msendleft(1)=msendleft(1)+1
              ileft(msendleft(1))=m
           endif
        endif
     enddo
     
#endif

  endif

  if (msend /= (msendleft(1)+msendright(1))) then
     write(*,*)'mype=',mype,'  msend NOT equal to msendleft+msendright'
     msend=msendleft(1)+msendright(1)
  endif

  if(iteration>1)then

! test: at least 1 particle needs to be shifted.
     mrecvtmp=0
     msendtmp=0
     if(msend>0)msendtmp=1

     call MPI_ALLREDUCE(msendtmp,mrecvtmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)

! no particle to be shifted, return
     if ( mrecvtmp == 0 ) then
!        write(0,*)istep,irk,mype,mp,m0,iteration
        return
     endif
  endif

! an extra space to prevent zero size when msendright(1)=msendleft(1)=0
  allocate(sendright(nzphase,max(1,msendright(1))),sendleft(nzphase,max(1,msendleft(1))),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj

! pack particle to move right
!$omp parallel do private(m)
  do m=1,msendright(1)
     sendright(1:nparam,m)=zion(1:nparam,iright(m))
     sendright(nparam+1:nzphase,m)=zion0(1:nparam,iright(m))
  enddo

! pack particle to move left
!$omp parallel do private(m)
  do m=1,msendleft(1)    
     sendleft(1:nparam,m)=zion(1:nparam,ileft(m))
     sendleft(nparam+1:nzphase,m)=zion0(1:nparam,ileft(m))
  enddo

  mtop=mp
! # of particles remain on local PE
  mp=mp-msendleft(1)-msendright(1)
! fill the hole
  lasth=msend
  do i=1,msend
     m=kzparticle(i)
     if (m > mp) exit  !Break out of the DO loop if m > mp
     do while(mtop == kzparticle(lasth))
        mtop=mtop-1
        lasth=lasth-1
     enddo
     zion(1:nparam,m)=zion(1:nparam,mtop)
     zion0(1:nparam,m)=zion0(1:nparam,mtop)
     mtop=mtop-1
     if (mtop == mp) exit  !Break out of the DO loop if mtop=mp
  enddo

! send # of particle to move right
  mrecvleft=0
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(msendright,2,MPI_INTEGER,idest,isendtag,&
       mrecvleft,2,MPI_INTEGER,isource,irecvtag,toroidal_comm,istatus,ierror)
  
  allocate(recvleft(nzphase,max(1,mrecvleft(1))),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
 
! send particle to right and receive from left
  recvleft=0.0
  isendcount=max(1,msendright(1))*nzphase
  irecvcount=max(1,mrecvleft(1))*nzphase
  call MPI_SENDRECV(sendright,isendcount,mpi_Rsize,idest,isendtag,recvleft,&
       irecvcount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
  
! send # of particle to move left
  mrecvright=0
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(msendleft,2,MPI_INTEGER,idest,isendtag,&
       mrecvright,2,MPI_INTEGER,isource,irecvtag,toroidal_comm,istatus,ierror)
  
  allocate(recvright(nzphase,max(1,mrecvright(1))),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj

! send particle to left and receive from right
  recvright=0.0
  isendcount=max(1,msendleft(1))*nzphase
  irecvcount=max(1,mrecvright(1))*nzphase
  call MPI_SENDRECV(sendleft,isendcount,mpi_Rsize,idest,isendtag,recvright,&
       irecvcount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
  
! need extra particle array
  if(mp+mrecvleft(1)+mrecvright(1) > mpmax)then
     write(*,*)"need bigger particle array",mype,mpmax,mp+mrecvleft(1)+mrecvright(1)
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! unpack particle, particle moved from left
!$omp parallel do private(m)
  do m=1,mrecvleft(1)
     zion(1:nparam,m+mp)=recvleft(1:nparam,m)
     zion0(1:nparam,m+mp)=recvleft(nparam+1:nzphase,m)
  enddo

! particle moved from right
!$omp parallel do private(m)
  do m=1,mrecvright(1)
     zion(1:nparam,m+mp+mrecvleft(1))=recvright(1:nparam,m)
     zion0(1:nparam,m+mp+mrecvleft(1))=recvright(nparam+1:nzphase,m)
  enddo
  
  mp=mp+mrecvleft(1)+mrecvright(1)
  
  deallocate(sendleft,sendright,recvleft,recvright)
  m0=mp-mrecvright(1)-mrecvleft(1)+1
  goto 100
  
end subroutine shifti

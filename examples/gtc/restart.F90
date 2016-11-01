#define RRESTART_CHECK
subroutine restart_io(iop)
  use global_parameters
  use particle_array
  use field_array
  use data_type
  use particle_tracking
  implicit none
  integer ierror, mrequest
  integer i,j,k,m,subsize,startidx,endidx,ndata
  integer mquantity,mflx,n_mode,mstepfinal,noutputs
  real(wp) dum
  character(*),intent(in)::iop
  character(len=50) :: restart_fname
  character(len=18) cdum
  integer :: save_restart_files
  INTEGER(KIND=MPI_OFFSET_KIND) mype_filesize, sum_filesize

#if ADIOS
integer*8 adios_handle, adios_groupsize, adios_totalsize, adios_err, adios_buf_size
    real*8 :: start_time, end_time, total_time,gbs
    integer :: sz
    real*8 :: lsz

  real(wp),dimension(:),allocatable::zion0_read,zelectron0_read
! THIS IS JUST FOR TESTNG
real(wp),dimension(:,:),allocatable::zion_read
integer :: err

#endif

  if(iop/="read" .and. iop/="write")then
     write(*,*)'*** subroutine restart_io (iop <> "read" or "write")',iop
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
     return
  endif
  
#if ADIOS
  !!xy rolling restart

!  write(restart_fname,'(a,i5.5,"_",i5.5,".bp")')"restart_dir/restart_",myrank_toroidal, irest 
  
  if(mod(irest,2)==0)then
     write(restart_fname,'(a,i5.5,".bp")')"restart_dir1/restart_",myrank_toroidal
  else
     write(restart_fname,'(a,i5.5,".bp")')"restart_dir2/restart_",myrank_toroidal
  endif
  
  ! setup the element path for this node
  !   call MPI_BARRIER(MPI_COMM_WORLD,ierror)
       start_time = MPI_WTIME()
  restart_fname=trim(restart_fname)//char(0)
  ! we include these new adios routines SAK

  if (iop=="write") then
     call adios_open (adios_handle, 'restart'//char(0), restart_fname, 'w'//char(0),adios_err)
#if RESTART_TIME
    call MPI_BARRIER(MPI_COMM_WORLD,err)
    start_time = MPI_WTIME()

#endif
#include "gwrite_restart.fh"
     call adios_close (adios_handle,adios_err)
#if RESTART_TIME
!    call MPI_BARRIER(MPI_COMM_WORLD,err)
    end_time = MPI_WTIME()
    total_time = end_time - start_time
    sz = adios_totalsize 
    write(6,200) mype,istep,total_time,(1.0*sz)/1024.0/1024.0
    call MPI_BARRIER(MPI_COMM_WORLD,err)
    call MPI_Gather(sz,1,MPI_INTEGER,data_size,1,MPI_INTEGER,0,mpi_comm_world,ierror)
    if (mype==0) then
    lsz = 0.0
    do i =1,nproc_toroidal
      lsz = lsz + data_size(i)*1.0
    end do
    lsz = lsz/1024.0/1024.0/1024.0
    gbs = (1.0*lsz) /total_time 
      write(*,199) total_time,lsz,gbs
    end if
199 format('restart time=',f10.2,f12.1,f12.2)
200 format(i4,i4,f12.2,f12.2)
	!    if(mype==0)write(*,*)'restart_time=',total_time
   
#endif
	!------------------------------------------------------------------------------
!                   SCOTT: THIS WILL CHECK THE WRITE, AND READ IT
!                          WE WILL VERIFY THAT WE READ IN WHAT WE WRITE!!!
#ifdef RESTART_CHECK
    !if (mype.eq.0) then
! please change the gread_restart.fh file and change the zion to zion_read
!    if(mype.eq.0) print *,nparam,mimax,'=nparam,mimax'
     allocate(zion_read(nparam,mimax),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
     allocate(zion0_read(mimax),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
     if(nhybrid>0)then
        allocate(zelectron0_read(memax),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
     end if
     if (mype.eq.0) write(6,*) 'opening up the restarts'
     call adios_open (adios_handle, 'restart'//char(0), restart_fname, 'r'//char(0), adios_err)
     write(6,*) restart_fname,'for',mype
     write(6,*) mype,mpsi,mgrid,nparam,memax,nhybrid
#include "gread_restart.fh"
     call adios_close (adios_handle,adios_err)

     do i=1,mi
        if (abs(zion0_read(i)-zion0(6,i))>1.0e-10) then
           write(10,*)'PE=', mype, ' i=',i, ' zion0_read=', zion0_read(i), ' zion0=',zion0(6,i)
           call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
        endif
        do j=1,nparam
            if (abs(zion_read(j,i)-zion(j,i))>1.0e-10) then
               write(10,*)'PE=', mype, ' i=',i,' j=',j,' zion_read=', zion_read(j,i), 'zion=',zion(j,i)
               call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
            endif
         enddo  
     enddo  


     if(allocated(zelectron0_read))then
        zelectron0(6,:)=zelectron0_read
        deallocate(zelectron0_read)
     endif
     zion0(6,:)=zion0_read
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
     if (mype.eq.0) write(6,*) '--------------------------------------------'
#endif

!     end if
!
!
!                      END OF CHECKING THE WRITE!!!
!------------------------------------------------------------------------------

     call restart_hist
!!***********************************************

  else
     allocate(zion0_read(mimax),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
     if(nhybrid>0)then
        allocate(zelectron0_read(memax),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
     end if
     call adios_open (adios_handle, 'restart'//char(0), restart_fname, 'r'//char(0), adios_err)
#include "gread_restart.fh"
     call adios_close (adios_handle,adios_err)
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)


     if(allocated(zelectron0_read))then
        zelectron0(6,:)=zelectron0_read
        deallocate(zelectron0_read)
     endif
     zion0(6,:)=zion0_read
     deallocate(zion0_read)
     if (mype.eq.0) write(*,*)'read in ',trim(restart_fname),mype
     irest=irest+1
  endif  ! end of read
#else
!---wj removed 2013-04-16, reinserted 2014-04-14
  call restart_native(iop)
#endif

end subroutine restart_io
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!the unformatted bin file
subroutine restart_native(iop)
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer ierror
  character(*),intent(in)::iop
  character(len=80) :: restart_fname
  character(len=19) cdum
  real(doubleprec) t1, t2 !---wj - changes to time i/o

  if(mype < 10)then
     write(cdum,'("DATA_RESTART.00000",i1)')mype
  elseif(mype < 100)then
     write(cdum,'("DATA_RESTART.0000",i2)')mype
  elseif(mype < 1000)then
     write(cdum,'("DATA_RESTART.000",i3)')mype
  elseif(mype < 10000)then
     write(cdum,'("DATA_RESTART.00",i4)')mype
  elseif(mype < 100000)then
     write(cdum,'("DATA_RESTART.0",i5)')mype
  else
     write(cdum,'("DATA_RESTART.",i6)')mype
  endif

!!XY rolling restart: save two copies of restart data
  if(mod(irest,2)==0)then
     restart_fname="restart_dir1/"//trim(cdum)
  else
     restart_fname="restart_dir2/"//trim(cdum)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)     

! record particle information for future restart run
  if(iop=="write")then
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)     
     t1 = mpi_wtime()
     open(222,file=trim(restart_fname),status='replace',form='unformatted')
     write(222)mi,me,rdtemi,pfluxi,phi,phip00,phi00,zonali,rtemi,rteme
     write(222)zion(1:nparam,1:mi),zion0(6,1:mi)
     if(nhybrid>0)write(222)rdteme,pfluxe,zonale,phisave,zelectron(1:nparam,1:me),zelectron0(6,1:me)
     if(magnetic>0)write(222)apara,fluidne,phiind,fluidue
     close(222)
     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
     t2 = mpi_wtime()
     if ( mype==0 ) write(gtcout,*) 'Time required to write restart files (sec): ', t2-t1
  if(mype==0)write(*,*)'write to ',trim(restart_fname)
  call restart_hist

!!***********************************************
! read particle information to restart previous run
  else
     if(mype==0)write(*,*)'read in ',trim(restart_fname)
     open(333,file=trim(restart_fname),status='old',form='unformatted')
     read(333)mi,me,rdtemi,pfluxi,phi,phip00,phi00,zonali,rtemi,rteme
     read(333)zion(1:nparam,1:mi),zion0(6,1:mi)
     if(nhybrid>0)read(333)rdteme,pfluxe,zonale,phisave,zelectron(1:nparam,1:me),zelectron0(6,1:me)
     if(magnetic>0)read(333)apara,fluidne,phiind,fluidue
     close(333)     
     irest=irest+1 !!XY rolling rstart
     if(mype==0)write(*,*)restart_fname,'read over'
     
  endif

end subroutine restart_native

subroutine restart_hist
  use global_parameters
  implicit none
  integer i,ii(5),mstepfinal,ndstep,ndata
  real dum
  character(len=80) :: restart_fname
!! write out two copies of history.out and data1d.out for restart from a crash
     if(mype==0 .and. istep<=mstep)then
        if(mod(irest,2)==0)then
           restart_fname="restart_dir1/"
        else
           restart_fname="restart_dir2/"
        endif        
#if RESTART_TIME
#else
! history.out
        open(333,file=trim(restart_fname)//'history_restart.out',status='replace')
        rewind(iodiag)
        read(iodiag,101)mstepfinal
        ndstep=mstepfinal-mstep/ndiag+istep/ndiag
        write(333,101)ndstep
        do i=1,5
           read(iodiag,101)ii(i)
          write(333,101)ii(i)
        enddo

        ndata=ii(1)*ii(2)+ii(3)*(2*ii(4)+ii(5))
        do i=0,ndata*ndstep
           read(iodiag,102)dum
          write(333,102)dum
        enddo
        close(333)
        write(*,*)'write to ',trim(restart_fname)//'history_restart.out'
! data1d.out        
        open(444,file=trim(restart_fname)//'data1d_restart.out',status='replace')
        rewind(iodata1d)
        read(iodata1d,101)mstepfinal
        ndstep=mstepfinal-mstep/ndiag+istep/ndiag
        write(444,101)ndstep
        do i=1,5
           read(iodata1d,101)ii(i)
           write(444,101)ii(i)
       enddo


        ndata=ii(1)*(ii(2)*ii(3)+ii(4)*ii(5))
        do i=1,ndata*ndstep
           read(iodata1d,102)dum
           write(444,102)dum
        enddo
        close(444)
        write(*,*)'write to ',trim(restart_fname)//'data1d_restart.out'
! save current location of restart files
        open(555,file="FileExit.out",status="replace")
        write(555,"(A9,i1)")"FileExit=",FileExit
        write(555,"(A9,i5)")"irest   =",irest+1
        if(mod(irest,2)==0)then
           write(555,103)"restart_dir1"
        else
           write(555,103)"restart_dir2"
        endif
103     format(A12)
        close(555)        

        if(istep==mstep)then
           close(iodiag)
           close(iodata1d)
        endif
#endif

     endif
     irest=irest+1 !!XY rolling rstart
101  format(i6)
102  format(e12.6)

end subroutine restart_hist

subroutine tag
  use global_parameters
  use particle_array
  use particle_tracking
  implicit none

! Each MPI process has its own "ptracked" array to accumulate the tracked
! particles that may or may not reside in its subdomain.
! The vector "ntrackp" contains the number of tracked particles currently residing on a MPI process

  ntrackp=0
  allocate(ptrackedi(nparam,max(mimax,1))) !!array to store ion tracked particles
  ptrackedi=0.0
  if(nhybrid>0)then
     allocate(ptrackede(nparam,max(memax,1)))  !!array to store tracked electrons
     ptrackede=0.0
  endif

! If irun==0 and if particle tracking is "on", tag each particle with a unique number
  if(irun==0)call tag_particles

end subroutine tag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine partout

#if ADIOS
     call locate_tracked_particles
     call write_tracked_particles
#endif

end subroutine partout

!========================================
!! originally written by s. either, modified by Y. Xiao April 2008
subroutine tag_particles
  
  !========================================
  
  use global_parameters
  use particle_array
  use particle_tracking
  implicit none
  integer :: i,m,np,npp
  integer :: ierror
  real(wp) :: r,ainside,aoutside,thickratio
  ! We keep track of the particles by tagging them with a unique number.
  ! We add an extra element to the particle array, which holds the
  ! particle tag, i.e. just a number.
  ! The input parameter "nptrack" is the total number of particles that
  ! we track.
  
  ! We tag each particle with a unique number that will be carried along
  ! with it as it moves from one processor/domain to another. To facilitate
  ! the search for the tracked particles, we give the non-tracked particles
  ! a negative value.
  ! CAVEAT: We are storing the tag in a floating point array, which means that
  !         for a large number of particles, the truncation error associated
  !         with the precision may make some particles have the same tag.
  !         This is particularly true in single precision where the maximum
  !         number of particles that can be distinguished with tags differing
  !         by unity is 2**24-1 = 16,777,215 (16.7 million). The most important
  !         is to distinguish the tracked particles between them. The other
  !         particles just need to have a negative tag. In order to minimize
  !         the effect of the precision truncation, we retag the tracked
  !         particles with positive numbers starting at 1.
  
  ! first allocate the np_all array...
  
  !  allocate(np_all(numberpe),offset(numberpe))
  do m=1,mi
     zion(7,m)=-real(m,wp) !!!-real(m+mype*mi)
     zion(8,m)=real(mype+1,wp)
     zion0(7,m)=-real(m,wp)
     zion0(8,m)=real(mype+1,wp)
  enddo
  
  if (nhybrid>0) then
     do m=1,me
        zelectron(7,m)=-real(m,wp) !!!-real(m+mype*mi)
        zelectron(8,m)=real(mype+1,wp)
        zelectron0(7,m)=-real(m,wp)
        zelectron0(8,m)=real(mype+1,wp)
     enddo
     
  endif
  
  thickratio=0.95
  ainside=0.5*(1.0-thickratio)*a1+0.5*(1.0+thickratio)*a0
  ainside=ainside*ainside*0.5
  
  aoutside=0.5*(1.0+thickratio)*a1+0.5*(1.0-thickratio)*a0
  aoutside=aoutside*aoutside*0.5
  
  if (nhybrid>0) then
     np=0
     do m=1,me
        
        if(zelectron(1,m)>ainside .and. zelectron(1,m)<aoutside) then
           np=np+1
           zelectron(7,m)=real(np,wp)
           zelectron0(7,m)=real(np,wp)
        endif
     enddo
     ! now communicate this out to all processors
     call MPI_Gather(np,1,MPI_INTEGER,offset,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
     if (mype==0) then
        np_total =  0
        do i=1,numberpe
           np_total = np_total + offset(i)
        end do
     end if
  else
     np=0
     do m=1,mi
        
        if(zion(1,m)>ainside .and. zion(1,m)<aoutside) then
           np=np+1
           zion(7,m)=real(np,wp)
           zion(7,m)=real(np,wp)
        endif
     enddo
     ! now communicate this out to all processors
     call MPI_Gather(np,1,MPI_INTEGER,offset,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
     if (mype==0) then
         np_total =  0
         do i=1,numberpe
            np_total = np_total + offset(i)
         end do
      end if
      
   endif
   
   if(mype==0)write(*,*)'np=',np,'ainside=',ainside,'aoutside=',aoutside
   
   call MPI_BCAST(np_total,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
   
   
   !!  if(mype==0)write(*,*)'zelectron(1,1:20)=',zelectron(1,1:20),'zelectron(1:8,1)=',zelectron(1:8,1)
   !! we only condiser those particles originally satify a certain condition
   ! namely within an annulous ain<r<aout
   ! We divide the number of tracked particles equally between all processors
   ! as much as possible. If nptrack is not a multiple of numberpe, we add
   ! 1 particle to some of the the processors until we get to nptrack. We also
   ! start at mi/2 so that the lower numbers fall around r=0.5 (see subroutine 
   ! load for the initial particle distribution in r).
   !  acentre=0.5*(a0+a1)
   !  awidth=0.1*(a1-a0)*0.5
   !  ainside= acentre-awidth
   !  aoutside=acentre+awidth
   
   ! np=0
   ! do m=1,mi
   !    r=sqrt(2.0*zion(1,m))
   !    if(r>ainside .and. r<aoutside) then
   !       np=np+1
   !       zion(7,m)=real(np,wp)
   !       zion0(7,m)=real(np,wp)
   !!    endif
   !! enddo
   
   !! if (nhybrid>0) then
   !!      np=0
   !!      do m=1,me
   !!         r=sqrt(2.0*zelectron(1,m))
   !!         if(r>ainside .and. r<aoutside) then
   !!            np=np+1
   !!            zelectron(7,m)=real(np,wp)
   !!            zelectron0(7,m)=real(np,wp)
   !!         endif
   !!      enddo
   !! endif
   
   !  if(mype==0)write(78,*)'np=',np,'  npp=',npp,'  zion(7,mi/2)=',zion(7,mi/2),&
   !                        '  zion(7,1)=',zion(7,1)
   
   !  if(mype==0)then
   !  ! On the master process (mype=0), we pick "nptrack" particles that will
   !  ! be followed at every time step. To facilitate the search of those
   !  ! particles among all the others in the zion arrays of each processor,
   !  ! we give them a positive number from 1 to nptrack. The particles are
   !  ! picked around r=0.5, which is mi/2.
   !    do m=(mi-nptrack)/2,(mi+nptrack)/2-1
   !       zion(7,m)=-zion(7,m)
   !       zion0(7,m)=-zion(7,m)
   !    enddo
   !  endif
   
 end subroutine tag_particles

!========================================
 
 subroutine locate_tracked_particles

!========================================
   
   use global_parameters
   use particle_array
   use particle_tracking
   implicit none
   integer :: i,m,npp
   
   ! Check if tracked particles are located on this processor. Particles that
   ! are tracked at every time step have a positive zion(7,m) value. All the
   ! others have a negative value.
   ntrackp=0
   ptrackedi=0.0
   
   npp=0
   !!  do m=1,mi
   !!    !! if(zion(7,m)>0.0)then
   !!       npp=npp+1
   !!       ptrackedi(1:nparam,npp)=zion(1:nparam,m)
   !!    !! endif
   !!  enddo
  ntrackp(1)=npp+1
  npp=0
  
  if (nhybrid>0) then
     ptrackede=0.0
     do m=1,me
        if(zelectron(7,m)>0.0)then
           npp=npp+1
           ptrackede(1:nparam,npp)=zelectron(1:nparam,m)
        endif
     enddo
     ntrackp(2)=npp
  else
     npp=0
     ptrackedi=0.0
     do m=1,mi
        if(zion(7,m)>0.0)then
           npp=npp+1
           ptrackedi(1:nparam,npp)=zion(1:nparam,m)
        endif
     enddo
     ntrackp(1)=npp
     
  endif
  !!  if(mype==0)write(*,*)'me=',me,'np=',npp,'zelectron(1,1:6)=',zelectron(7,1:6),&
  !!                        'zelectron(1:8,1)= ',zelectron(1:8,1)
  
end subroutine locate_tracked_particles

!========================================

subroutine write_tracked_particles
  
  !========================================
  
  use global_parameters
  use particle_tracking
  implicit none
  
  integer :: i,j,ntpart(0:numberpe-1)
  character(len=10) :: cdum
#if ADIOS
  character(len=50),SAVE:: fname
  character(len=50)::dirstr
  integer*8 :: group_id, buf_id,comm, ierror
  integer*8 :: adios_handle, adios_groupsize, adios_totalsize,adios_err
  integer ntrack_sort, g_iteration
#endif
  
  
#if ADIOS
  comm=MPI_COMM_WORLD
  write(fname,'("trackp_dir/TRACKP_",i5.5,".bp")')mstepall+istep  
  fname=trim(fname)//char(0)
  if (nhybrid>0) then
     ntrack_sort = ntrackp(2)
     call sort_ptrackede(ptrackede(:,1:ntrackp(2)), nparam, ntrack_sort, numberpe, mstepall, istep, comm, adios_err)     
     ntrackp(2) = ntrack_sort
     ! communicate the number of particles through the MPI_Allgather
     call MPI_Gather(ntrackp(2),1,MPI_INTEGER,np_all,1,MPI_INTEGER,0,comm,ierror)
     if (mype==0) then
        offset(1) = 0
        my_offset = 0
        do i=2,numberpe
           offset(i) =np_all(i-1) +  offset(i-1)
        end do
     end if
     call MPI_Scatter(offset,1,MPI_INTEGER,my_offset,1,MPI_INTEGER,0,comm,ierror)
     ! now my_offset is the offset in the big particle array....
     call adios_open (adios_handle, "particles"//char(0), fname, "w"//char(0),adios_err)
     !call adios_set_path (adios_handle,dirstr//char(0),adios_err)
#include "gwrite_particles.fh"     
     call adios_close (adios_handle,adios_err)     
     if(mype==0)write(*,*)"tracking filename ",fname     
  else !ions
     ntrack_sort = ntrackp(1)
     call sort_ptrackede(ptrackedi(:,1:ntrackp(1)), nparam, ntrack_sort, numberpe, mstepall, istep, comm, adios_err)     
     ntrackp(1) = ntrack_sort
     ! communicate the number of particles through the MPI_Allgather
     call MPI_Gather(ntrackp(1),1,MPI_INTEGER,np_all,1,MPI_INTEGER,0,comm,ierror)
     if (mype==0) then
        offset(1) = 0
        my_offset = 0
        do i=2,numberpe
           offset(i) =np_all(i-1) +  offset(i-1)
        end do
     end if
     call MPI_Scatter(offset,1,MPI_INTEGER,my_offset,1,MPI_INTEGER,0,comm,ierror)
     ! now my_offset is the offset in the big particle array....
     call adios_open (adios_handle, "particles"//char(0), fname, "w"//char(0),adios_err)
     !call adios_set_path (adios_handle,dirstr//char(0),adios_err)
#include "gwrite_particles.fh"     
     call adios_close (adios_handle,adios_err)     
     if(mype==0)write(*,*)"tracking filename ",fname
  end if
#else
  
     write(cdum,'("trackp_dir/TRACKP.",i5.5)')mype
     open(57,file=cdum,status='unknown',position='append')
     write(57,*)istep
     write(57,*)ntrackp(1:nspecies)

     do j=1,ntrackp(1)
        write(57,*)ptrackedi(1:nparam,j)
     enddo
     if(nhybrid>0)then
        do j=1,ntrackp(2)
           write(57,*)ptrackede(1:nparam,j)
        enddo
     endif

     close(57)
#endif

end subroutine write_tracked_particles

!========================================


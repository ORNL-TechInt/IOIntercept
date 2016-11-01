!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine petsc_init
      use global_parameters,only:mype
      use particle_decomp,only:partd_comm
      use petsc_array,only: newcomm,nproc_newcomm,myrank_newcomm

      implicit none 

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscviewer.h"

      PetscErrorCode ierr
      integer ierror

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

      call MPI_Comm_dup(partd_comm,newcomm,ierror)

      call mpi_comm_size(newcomm,nproc_newcomm,ierror)

      call mpi_comm_rank(newcomm,myrank_newcomm,ierror)

      if(ierror /= 0) then
         write(*,*)mype,'main_petsc: newcomm error',ierror
         call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
      endif
      
      return
      end subroutine petsc_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine petsc_final
      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscis.h"

      Vec       x,b,x_recv
      Mat       Apetsc1
      KSP       ksp
      PC        pc
      IS        is_from,is_to
      VecScatter     gather

      common/petsc_common/Apetsc1,x,b,ksp,pc,x_recv,is_from,is_to,gather

      PetscErrorCode ierr

      call VecDestroy(x,ierr)
      call VecDestroy(b,ierr)
      call VecDestroy(x_recv,ierr)
      call VecScatterDestroy(gather,ierr)
      call ISDestroy(is_from,ierr)
      call ISDestroy(is_to,ierr)
      call MatDestroy(Apetsc1,ierr)
      call KSPDestroy(ksp,ierr)

      call PetscFinalize(ierr)

      return
      end subroutine petsc_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine petsc_create(mindex)
      use global_parameters,only: mgrid
      use petsc_array, only:newcomm,usera,userj,users
      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"

      Vec       x,b,x_recv
      Mat       Apetsc1
      KSP       ksp
      PC        pc
      IS        is_from,is_to
      VecScatter     gather

      common/petsc_common/Apetsc1,x,b,ksp,pc,x_recv,is_from,is_to,gather

      integer mindex
      double precision tol,value(mindex)
      PetscErrorCode ierr
      PetscInt n,i,j,ij,ilow,ihigh,nz
      PetscTruth flg
      integer col(mindex)

      n=mgrid

      call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',n,flg,ierr)
      call VecCreate(newcomm,x,ierr)
      call VecSetSizes(x,PETSC_DECIDE,n,ierr)
      call VecSetFromOptions(x,ierr)
      call VecDuplicate(x,b,ierr)

!*
      call VecCreate(PETSC_COMM_SELF,x_recv,ierr)
      call VecSetSizes(x_recv,n,n,ierr)
      call VecSetType(x_recv,VECSEQ,ierr)
      call VecSetFromOptions(x_recv,ierr)

      call ISCreateStride(newcomm,n,0,1,is_from,ierr)
      call ISCreateStride(PETSC_COMM_SELF,n,0,1,is_to,ierr)

      call VecScatterCreate(x,is_from,x_recv,is_to,gather,ierr)
!*

!     call MatCreate(newcomm,Apetsc1,ierr)
!     CHKERRQ(ierr)
!     call MatSetSizes(Apetsc1,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
!     CHKERRQ(ierr)
!     call MatSetType(Apetsc1,MATMPIAIJ,ierr)
!     CHKERRQ(ierr)
!*optimized matrix creation
      call MatCreateMPIAIJ(newcomm,PETSC_DECIDE,PETSC_DECIDE,n,n,
     .mindex,PETSC_NULL_INTEGER,mindex,PETSC_NULL_INTEGER,Apetsc1,ierr)
      call MatSetFromOptions(Apetsc1,ierr)
      call MatGetOwnershipRange(Apetsc1,ilow,ihigh,ierr)

      ij=users(ilow-1)
      do i=ilow,ihigh-1
        nz=users(i)-users(i-1)
        do j=1,nz
          ij=ij+1
          col(j) = userj(ij)
          value(j) = usera(ij)
        enddo
        call MatSetValues(Apetsc1,1,i,nz,col,value,INSERT_VALUES,ierr)
      enddo

      call MatAssemblyBegin(Apetsc1,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(Apetsc1,MAT_FINAL_ASSEMBLY,ierr)
!     call MatTranspose(Apetsc1,Bpetsc1,ierr)
!     call MatAXPY(Bpetsc1,1.d0,Apetsc1,DIFFERENT_NONZERO_PATTERN,ierr)
!     call MatPartitioningSetAdjacency(mpart, Bpetsc1, ierr )
!     call MatCopy(Apetsc1,Bpetsc1,DIFFERENT_NONZERO_PATTERN,ierr)


      call KSPCreate(newcomm,ksp,ierr)
      call KSPSetOperators(ksp,Apetsc1,Apetsc1,
     .                     DIFFERENT_NONZERO_PATTERN,ierr)
!     call KSPGetPC(ksp,pc,ierr)
!     call PCSetType(pc,PCJACOBI,ierr)
      tol = 1.d-4
!     tol = 1.d-20
      call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_DOUBLE_PRECISION,     &
     &     PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)
      call KSPSetFromOptions(ksp,ierr)

      return
      end subroutine petsc_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine petsc_solver
      use global_parameters,only:mgrid,istep,irk,mype
      use petsc_array,only:userx,userb
      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscdef.h"
#include "include/finclude/petscis.h"
 
      Vec       x,b,x_recv
      Mat       Apetsc1
      KSP       ksp
      PC        pc
      IS        is_from,is_to
      VecScatter     gather

      common/petsc_common/Apetsc1,x,b,ksp,pc,x_recv,is_from,is_to,gather

      double precision temp
      PetscErrorCode ierr
      PetscInt i,n,ilow,ihigh,isize,ilo,ihi
      PetscInt isi
      PetscOffset idx
      PetscScalar vec_x(1)

      n=mgrid

      call VecGetLocalSize(b,isize,ierr)
      call VecGetOwnershipRange(b,ilow,ihigh,ierr)
      do i=ilow,ihigh-1
        temp=userb(i)
        call VecSetValues(b,1,i,temp,INSERT_VALUES,ierr)
      enddo
      call VecAssemblyBegin(b,ierr)
      call VecAssemblyEnd(b,ierr)

!     if(istep==1 .and. irk==1 .and. mype==0)then
!       write(*,*)'checking vector b:'
!       call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)
!       CHKERRQ(ierr)
!     endif
!     if(istep==1 .and. irk==1 .and. mype==0)then
!       write(*,*)'checking matrix a:'
!       call MatView(Apetsc1,PETSC_VIEWER_STDOUT_SELF,ierr)
!       CHKERRQ(ierr)
!     endif

      call KSPSetInitialGuessNonzero(ksp,Petsc_True,ierr)
      call KSPSolve(ksp,b,x,ierr)

!* 
      call VecScatterBegin(gather,x,x_recv,INSERT_VALUES,
     .  SCATTER_FORWARD,ierr)
      call VecScatterEnd(gather,x,x_recv,INSERT_VALUES,
     .  SCATTER_FORWARD,ierr)
!*
      call VecGetOwnershipRange(x_recv,ilo,ihi,ierr)
      call VecGetLocalSize(x_recv,isi,ierr)
      call VecGetArray(x_recv,vec_x,idx,ierr)
      do i=ilo,ihi-1
        userx(i)=vec_x(idx+i-ilo+1)
      enddo
      call VecRestoreArray(x,vec_x,idx,ierr)

      return
      end subroutine petsc_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

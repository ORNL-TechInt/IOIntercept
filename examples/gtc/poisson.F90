subroutine poisson_solver(scenario)
  use precision
  use global_parameters
  use particle_array
  use field_array
  use petsc_array,only: userb,userx
  implicit none

  integer :: i,it,ij,j,k,n,iteration,mring,mindex,mtest,ierr
  integer,dimension(:),allocatable :: nindex
  integer,dimension(:,:),allocatable :: indexp
  real(wp),dimension(:,:),allocatable :: ring
  real(wp) gamma,tmp,prms,perr(mgrid)
  real(wp) ptilde(mgrid),phitmp(mgrid),dentmp(mgrid)
  character(*),intent(in) :: scenario

  save nindex,indexp,ring

  interface
     subroutine field_solver_initial(mring,mindex,nindex,indexp,ring)
       use precision
       use global_parameters,only: mgrid
       implicit none
       
       integer :: mring,mindex
       integer,dimension(mgrid) :: nindex
       integer,dimension(mindex,mgrid) :: indexp
       real(wp),dimension(mindex,mgrid) :: ring
     end subroutine field_solver_initial
  end interface

! number of gyro-ring
  mring=2

! number of summation: maximum is 32*mring+1
  mindex=32*mring+1

! gamma=0.75: max. resolution for k=0.577
  gamma=0.75
  iteration=5

  tmp=1.0/(2.0-gamma) !Assuming T_i=T_e

! initialize poisson solver
  if(istep==1 .and. irk==1 .and. scenario=='adiabatic-electron')then
    allocate(indexp(mindex,mgrid),ring(mindex,mgrid),nindex(mgrid),STAT=mtest)
    if (mtest /= 0) then
        write(0,*)mype,'*** Cannot allocate indexp: mtest=',mtest
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    endif

!*initialize
    call poisson_solver_initial(mring,mindex,nindex,indexp,ring)
  endif

  if(scenario=='adiabatic-electron')then
!$omp parallel do private(i)
     do i=1,mgrid
        dentmp(i)=qion*densityi(1,i)
     enddo
  else
!$omp parallel do private(i)
     do i=1,mgrid
        dentmp(i)=qion*densityi(1,i)+qelectron*densitye(1,i)
     enddo
  endif

#ifdef __PETSc

!$omp parallel do private(i)
  do i=1,mgrid
     userb(i-1)=dentmp(i)
  enddo
  
  call petsc_solver

!$omp parallel do private(i)
  do i=1,mgrid
     phi(1,i)=userx(i-1)
  enddo

! radial boundary
  phi(1,igrid(0):igrid(0)+mtheta(0))=0.0
  phi(1,igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0

#else

! first iteration, first guess of phi. (1+T_i/T_e) phi - phi_title = n_i
!$omp parallel do private(i)
  do i=1,mgrid
     phitmp(i)=dentmp(i)*tmp
  enddo

  do it=2,iteration
!$omp parallel do private(i,j)
     do i=1,mgrid
        ptilde(i)=0.0
        do j=1,nindex(i)
           ptilde(i)=ptilde(i)+ring(j,i)*phitmp(indexp(j,i))
        enddo
     enddo

!$omp parallel do private(i)
     do i=1,mgrid
        perr(i)=ptilde(i)-gamma*phitmp(i)
        phitmp(i)=(dentmp(i)+perr(i))*tmp
     enddo
!        prms=sum(perr*perr)/sum(phitmp*phitmp)
!        if(mype==0)write(*,*)istep,prms

! radial boundary
     phitmp(igrid(0):igrid(0)+mtheta(0))=0.0
     phitmp(igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0
  enddo

! store final results
!$omp parallel do private(i)
  do i=1,mgrid
     phi(1,i)=phitmp(i)
  enddo

#endif

! in equilibrium unit
!$omp parallel do private(i)
  do i=0,mpsi
     phi(1,igrid(i)+1:igrid(i)+mtheta(i))=phi(1,igrid(i)+1:igrid(i)+&
          mtheta(i))*rtemi(i)*(qion*rho0)**2/aion
! poloidal BC
     phi(1,igrid(i))=phi(1,igrid(i)+mtheta(i))
  enddo
end subroutine poisson_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poisson_solver_initial(mring,mindex,nindex,indexp,ring)
  use precision
  use global_parameters
  use field_array
  use petsc_array,only: userj,users,usera,userb,userx
  implicit none

  integer mring,mindex,ierr,nindex0max
  integer,dimension(mgrid) :: nindex
  integer,dimension(mindex,mgrid) :: indexp
  real(wp),dimension(mindex,mgrid) :: ring        
  integer i,ii,ij,ij0,ipjt,j,jt,j1,j0,jm,k,kr,kp,nt,np,minn,maxn
  real(wp) vring(3),fring(3),rg,tgrid,ddelr,ddelt,wght,r,rr,t,rdum,wr,&
       wt1,wt0,tdum,zdum,b,pi2_inv,delr,delt(0:mpsi)
  
#ifdef __PETSc
  interface
     subroutine matrix(mindex,nindex,indexp,ring,nindex0max)
       use precision
       use global_parameters,only: mgrid
       implicit none

       integer :: mindex,nindex0max
       integer,dimension(mgrid) :: nindex
       integer,dimension(mindex,mgrid) :: indexp
       real(wp),dimension(mindex,mgrid) :: ring
     end subroutine matrix
     
     subroutine petsc_create(mindex)
       implicit none
       integer mindex
     end subroutine petsc_create
  end interface
#endif

  if(mring==1)then
! one ring, velocity in unit of rho0
     vring(1)=sqrt(2.0)
     fring(1)=1.0
     
  elseif(mring==2)then
! two rings good for up to k_perp=1.5
     vring(1)=0.9129713024553
     vring(2)=2.233935334042
     fring(1)=0.7193896325719
     fring(2)=0.2806103674281
     
  else      
! three rings: exact(<0.8%) for up to k_perp=1.5
     vring(1)=0.388479356715
     vring(2)=1.414213562373
     vring(3)=2.647840808818
     fring(1)=0.3043424333839
     fring(2)=0.5833550690524
     fring(3)=0.1123024975637
  endif

  pi2_inv=0.5/pi
  delr=1.0/deltar
  delt=2.0*pi/deltat
  nindex=0
  ring=0.0
  indexp=1
  do i=0,mpsi
     do j=1,mtheta(i)
        ij0=igrid(i)+j

! 1st point is the original grid point
        nindex(ij0)=1
        indexp(1,ij0)=ij0
        ring(1,ij0)=0.25

! position of grid points
        rg=a0+deltar*real(i)
        tgrid=deltat(i)*real(j)+zeta1*qtinv(i)
        tgrid=tgrid*pi2_inv
        tgrid=2.0*pi*(tgrid-aint(tgrid))
        jt=max(0,min(mtheta(i),int(pi2_inv*delt(i)*tgrid+0.5)))
        
! B-field
        b=1.0/(1.0+rg*cos(tgrid))
        ipjt=igrid(i)+jt
        
        do kr=1,mring

! FLR from grid point and weight of 8-point for each ring. Counter-clockwise: 1,5,2,6,3,7,4,8
           do kp=1,8
              if(kp<5)then
                 ddelr=pgyro(kp,ipjt)
                 ddelt=tgyro(kp,ipjt)
                 wght=0.0625*fring(kr)
                 
              elseif(kp==5)then
                 ddelr=0.5*(pgyro(1,ipjt)+pgyro(2,ipjt))
                 ddelt=0.5*(tgyro(1,ipjt)+tgyro(2,ipjt))
                 wght=0.125*fring(kr)
              elseif(kp==6)then
                 ddelr=0.5*(pgyro(2,ipjt)+pgyro(3,ipjt))
                 ddelt=0.5*(tgyro(2,ipjt)+tgyro(3,ipjt))
              elseif(kp==7)then
                 ddelr=0.5*(pgyro(3,ipjt)+pgyro(4,ipjt))
                 ddelt=0.5*(tgyro(3,ipjt)+tgyro(4,ipjt))
              elseif(kp==8)then
                 ddelr=0.5*(pgyro(4,ipjt)+pgyro(1,ipjt))
                 ddelt=0.5*(tgyro(4,ipjt)+tgyro(1,ipjt))
              endif
              
! position for each point with rho_i=2.0*vring
              r=rg+ddelr*2.0*vring(kr)*sqrt(0.5/b)
              t=tgrid+ddelt*2.0*vring(kr)*sqrt(0.5/b)

! linear interpolation
              rdum=delr*max(0.0_wp,min(a1-a0,r-a0))
              ii=max(0,min(mpsi-1,int(rdum)))
              wr=rdum-real(ii)
              if(wr>0.95)wr=1.0
              if(wr<0.05)wr=0.0

! outer flux surface
              tdum=t-zeta1*qtinv(ii+1)
              tdum=tdum*pi2_inv+10.0
              tdum=delt(ii+1)*(tdum-aint(tdum))
              j1=max(0,min(mtheta(ii+1)-1,int(tdum)))
              wt1=tdum-real(j1)
              if(wt1>0.95)wt1=1.0
              if(wt1<0.05)wt1=0.0

! inner flux surface
              tdum=t-zeta1*qtinv(ii)
              tdum=tdum*pi2_inv+10.0
              tdum=delt(ii)*(tdum-aint(tdum))
              j0=max(0,min(mtheta(ii)-1,int(tdum)))
              wt0=tdum-real(j0)
              if(wt0>0.95)wt0=1.0
              if(wt0<0.05)wt0=0.0
              
! index and weight of each point
              do np=1,4
                 if(np==1)then
                    ij=igrid(ii+1)+j1+1
                    rr=wght*wr*wt1
                 elseif(np==2)then
                    if(j1==0)j1=mtheta(ii+1)
                    ij=igrid(ii+1)+j1
                    rr=wght*wr*(1.0-wt1)
                 elseif(np==3)then
                    ij=igrid(ii)+j0+1
                    rr=wght*(1.0-wr)*wt0
                 else
                    if(j0==0)j0=mtheta(ii)
                    ij=igrid(ii)+j0
                    rr=wght*(1.0-wr)*(1.0-wt0)
                 endif

! insignificant point replaced by the original grid point
                 if(rr<0.001)then
                    ring(1,ij0)=ring(1,ij0)+rr
                    goto 100
                 endif

                 do nt=1,nindex(ij0)
! redundant point
                    if(ij==indexp(nt,ij0))then
                       ring(nt,ij0)=ring(nt,ij0)+rr
                       goto 100
                    endif
                 enddo
! new point
                 nindex(ij0)=nindex(ij0)+1
                 nt=nindex(ij0)
                 indexp(nt,ij0)=ij
                 ring(  nt,ij0)=rr
                
100             continue
              enddo  !end of 4-point interpolation loop
           enddo     !end of 8-point-per-ring loop
        enddo        !end of ring loop
     enddo           !end of poloidal loop
  enddo              !end of radial loop     
  
! check array size
  if(maxval(nindex)>mindex)then
     write(gtcout,*)'Poisson error',mype,maxval(nindex),' > ',mindex
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
  endif
  
  rdum=0.0
  tdum=0.0
  zdum=1.0
  do i=1,mgrid
     rdum=sum(ring(1:nindex(i),i))
     tdum=max(tdum,rdum)
     zdum=min(zdum,rdum)
  enddo
  if(mype==0)then
     write(gtcout,*)'poisson solver=',maxval(nindex),minval(nindex),&
          tdum,zdum,mgrid,sum(nindex)
  end if

#ifdef __PETSc
  call matrix(mindex,nindex,indexp,ring,nindex0max)

  call petsc_create(nindex0max)
#endif
end subroutine poisson_solver_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gyroinitial
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer i,j,ij,ngyromax
  real(wp) r,tdum,q,b,dtheta_dx,rhos

! for N-point gyro-averaging, N=4 for now. Anothe four point could be added for N=8
  ngyromax=ngyroi
  allocate(pgyro(ngyromax,mgrid),tgyro(ngyromax,mgrid),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj

! guiding center approximation; no gyro-averaging
  if(ngyromax==1)then
     pgyro=0.0
     tgyro=0.0
  else
! 4-point gyro-averaging on grids; starting from outermost point, counter-clockwise: 1,2,3,4.
! Calculate coeffecients here for rho_local=rho0*sqrt(2B_0/B), _0 represent quantities on axis.
! In subroutines poisson and locatei: rho_particle=rho_local*sqrt(m_i/m_proton)*sqrt(mu)/q_i.
! dtheta/delta_x=1/(r*(1+r*cos(theta))), delta_x=poloidal length increase.
!$omp parallel do private(i,j,r,ij,tdum,q,b,dtheta_dx,rhos)
     do i=0,mpsi
        r=a0+deltar*real(i)
        do j=0,mtheta(i)
           ij=igrid(i)+j
           tdum=deltat(i)*real(j)
           q=q0+q1*r/aminor+q2*r*r/(aminor*aminor)
           b=1.0/(1.0+r*cos(tdum))
           dtheta_dx=1.0/r
! first two points perpendicular to field line on poloidal surface            
           rhos=sqrt(2.0/b)*rho0
           pgyro(1,ij)=rhos
           pgyro(3,ij)=0.0-pgyro(1,ij)
! non-orthorgonality between psi and theta: tgyro=-rhos*dtheta_dx*r*sin(tdum)
           tgyro(1,ij)=0.0
           tgyro(3,ij)=0.0-tgyro(1,ij)

! the other two points tangential to field line
           tgyro(2,ij)=rhos*dtheta_dx
           tgyro(4,ij)=0.0-tgyro(2,ij)
           pgyro(2,ij)=rhos*0.5*rhos/r
           pgyro(4,ij)=pgyro(2,ij)

! add another 4-point for 8-point averaging; Starting from 1st quadrant, counter-clockwise: 5,6,7,8
           if(ngyromax==8)then
              pgyro(5,ij)=0.707*pgyro(1,ij)
              tgyro(5,ij)=0.707*tgyro(2,ij)
              pgyro(6,ij)=0.707*pgyro(3,ij)
              tgyro(6,ij)=0.707*tgyro(2,ij)
              pgyro(7,ij)=0.707*pgyro(3,ij)
              tgyro(7,ij)=0.707*tgyro(4,ij)
              pgyro(8,ij)=0.707*pgyro(1,ij)
              tgyro(8,ij)=0.707*tgyro(4,ij)
           endif
        enddo
     enddo
  endif

end subroutine gyroinitial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix(mindex,nindex,indexp,ring,nindex0max)
!*Using nindex,indexp,ring generate matrices users, usera, and userj
  use precision
  use global_parameters
  use field_array
  use petsc_array,only: usera,userx,userb,users,userj
  implicit none

  integer :: mindex,nnz,windexp,msort,nindex0max
  integer,dimension(mgrid) :: nindex,nindex0
  integer,dimension(mindex,mgrid) :: indexp,indexp0
  real(wp),dimension(mindex,mgrid) :: ring,ring0
  
  integer :: i,j,k,ij0,ierror,mtest
  real(wp) :: wring,diagonal

!*local copy of variables
  nindex0=nindex
  indexp0=indexp
  ring0=ring

!*poloidal BC
  do i=0,mpsi
    nindex0(igrid(i))=nindex0(igrid(i)+mtheta(i))
    indexp0(:,igrid(i))=indexp0(:,igrid(i)+mtheta(i))
    ring0(:,igrid(i))=ring0(:,igrid(i)+mtheta(i))
    indexp0(1,igrid(i))=indexp0(1,igrid(i))-mtheta(i)  !*"j=0" refer to self
  enddo

  diagonal=2.0-real(magnetic)
  do i=1,mgrid
    ring0(1,i)=diagonal-ring0(1,i)
    do j=2,nindex0(i)
      ring0(j,i)=0.0-ring0(j,i)
    enddo
  enddo

!*radial conditions.
  do i=0,mpsi,mpsi
    do j=0,mtheta(i)
      ij0=igrid(i)+j
      nindex0(ij0)=1
      indexp0(1,ij0)=ij0
      ring0(1,ij0)=1.0
    enddo
  enddo

!*sort indexp and ring
  do i=1,mgrid
    do msort=nindex0(i)-1,1,-1
      do j=1,msort
        if(indexp0(j,i)>indexp0(j+1,i)) then
          windexp=indexp0(j,i)
          wring=ring0(j,i)
          indexp0(j,i)=indexp0(j+1,i)
          ring0(j,i)=ring0(j+1,i)
          indexp0(j+1,i)=windexp
          ring0(j+1,i)=wring
        endif
      enddo
    enddo
  enddo

!*Count nonzero and allocate
  nindex0max=0
  nnz=0
  do i=1,mgrid
    if(nindex0(i).gt.nindex0max) nindex0max=nindex0(i)
    do j=1,nindex0(i)
      nnz=nnz+1
    enddo
  enddo

  allocate(usera(nnz),userj(nnz),userb(0:mgrid-1),users(-1:mgrid-1),&
    userx(0:mgrid-1),stat=mtest)
  if (mtest /= 0) then
    write(0,*)mype,'matrix: Cannot allocate userX'
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

!*Make matrix usera,userj
  nnz=0
  users(-1)=0
  do i=1,mgrid
    users(i-1)=users(i-2)+nindex0(i)
    do j=1,nindex0(i)
      nnz=nnz+1
      usera(nnz)=ring0(j,i)
      userj(nnz)=indexp0(j,i)-1
    enddo
  enddo
end subroutine matrix

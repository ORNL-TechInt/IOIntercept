! Load or restart all particle species
subroutine load
  use global_parameters
  use particle_tracking

! initialize particle position and velocity: ions first, then electron
  CALL LOADI
  if(nhybrid>0)CALL LOADE
  if(irun>0)CALL RESTART_LOAD


!!tag particles for diagnosis and initialize ptrack arrays
  if(track_particles==1) then 
     allocate(np_all(numberpe),offset(numberpe))
     CALL TAG
  end if

#if RESTART_TIME
    if (mype==0) write(6,*) 'timing the restarts',micell
    allocate(data_size(numberpe))
#endif
end subroutine load
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! load the first ion species
subroutine loadi
  use global_parameters
  use particle_array
  use field_array
  use particle_tracking
  implicit none

  integer i,j,ij,m,ierr,iorestart,mtest,ierror
  real(wp) :: c0=2.515517,c1=0.802853,c2=0.010328,d1=1.432788,d2=0.189269,d3=0.001308
  real(wp) w_initial,energy,momentum,r,rmi,pi2_inv,delr,ainv,&
       vthi,cost,b,upara,eperp,cmratio,z4tmp,en1,en2,zdum,tdum,rmax,rmin,umax
  integer*8 eightbyteint
  
! allocate memory
  mimax=mi+100*ceiling(sqrt(real(mi))) !ions array upper bound

  allocate(zion(nparam,mimax),zion0(nparam,mimax),pmarki(0:mpsi),markeri(mgrid),&
       densityi(0:1,mgrid),currenti(0:1,mgrid),rtemi(0:mpsi),pfluxi(0:mpsi),rteme(0:mpsi),&
       rdtemi(0:mpsi),zonali(0:mpsi),zonic(0:mpsi),rden(0:mpsi),data1di(0:mpsi,mpdata1d),&
       wzion(mimax),wpion(ngyroi,mimax),&
       wtion0(ngyroi,mimax),wtion1(ngyroi,mimax),STAT=mtest,SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
  if(mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate ion: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  data1di = 0. !---wj

  allocate(jtion0(ngyroi,mimax),jtion1(ngyroi,mimax),STAT=mtest) !---wj
  if(mtest /= 0) then !---wj
     write(0,*)mype,'*** Cannot allocate ion: mtest=',mtest !---wj
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror) !---wj
  endif !---wj

  rtemi=1.0
  rteme=1.0
  rden=1.0
  pfluxi=0.0
  rdtemi=0.0
  zonali=0.0
  zonic=0.0
  markeri=0.0
  densityi=0.0
  currenti=0.0

  eightbyteint = 1

! number of marker per grid, Jacobian=(1.0+r*cos(theta+r*sin(theta)))*(1.0+r*cos(theta))
  pmarki=0.0
!$omp parallel do private(i,j,r,ij,zdum,tdum,rmax,rmin)
  do i=0,mpsi
     r=a0+deltar*real(i)
     do j=1,mtheta(i)
        ij=igrid(i)+j
	zdum=zeta1
        tdum=real(j)*deltat(i)+zdum*qtinv(i)
        markeri(ij)=(1.0+r*cos(tdum))**2
        pmarki(i)=pmarki(i)+markeri(ij)
     enddo
     rmax=min(a1,r+0.5*deltar)
     rmin=max(a0,r-0.5*deltar)
     tdum=real(mi*eightbyteint*npartdom)*(rmax*rmax-rmin*rmin)/(a1*a1-a0*a0)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        markeri(ij)=tdum*markeri(ij)/pmarki(i)
        markeri(ij)=1.0/markeri(ij) !to avoid divide operation
     enddo
     pmarki(i)=1.0/(real(mtoroidal)*tdum)
     markeri(igrid(i))=markeri(igrid(i)+mtheta(i))
  enddo

  if(irun>1)return

  rmi=1.0/real(mi*eightbyteint*npartdom)
  pi2_inv=0.5/pi
  delr=1.0/deltar
  ainv=1.0/aminor
  w_initial=1.0e-3
  if(magnetic==1)w_initial=0.0
  umax=4.0

! radial: uniformly distributed in r^2, later transform to psi
!$omp parallel do private(m)
  do m=1,mi
!!   zion(1,m)=sqrt(a0*a0+(real(m)-0.5)*(a1*a1-a0*a0)*rmi)
     zion(1,m)=sqrt(a0*a0+(real((m-1)*eightbyteint*npartdom+myrank_partd+1)-0.5)*(a1*a1-a0*a0)*rmi)
  enddo

! Set zion(2:6,1:mi) to uniformly distributed random values between 0 and 1
  call random_number(zion(2,1:mi))
  call random_number(zion(3,1:mi))
  call random_number(zion(4,1:mi))
  call random_number(zion(5,1:mi))
  call random_number(zion(6,1:mi))

! poloidal: uniform in alpha=theta_0+r*sin(alpha_0), theta_0=theta+r*sin(theta)
!$omp parallel do private(m)
  do m=1,mi                
     zion(2,m)=2.0*pi*(zion(2,m)-0.5)
     zion0(2,m)=zion(2,m) !zion0(2,:) for temporary storage
  enddo
  do i=1,10
!$omp parallel do private(m)
     do m=1,mi
        zion(2,m)=zion0(2,m)-2.0*zion(1,m)*sin(zion(2,m))
     enddo
  enddo
!$omp parallel do private(m)
  do m=1,mi                
     zion(2,m)=zion(2,m)*pi2_inv+10.0 !period of 1
     zion(2,m)=2.0*pi*(zion(2,m)-aint(zion(2,m))) ![0,2*pi)
  enddo

! Maxwellian distribution in v_para, <v_para^2>=1.0, use zion0(4,:) as temporary storage
!$omp parallel do private(m,z4tmp)
  do m=1,mi                
     z4tmp=zion(4,m)
     zion(4,m)=zion(4,m)-0.5
     zion0(4,m)=sign(1.0_wp,zion(4,m))
     zion(4,m)=sqrt(max(1.0e-20_wp,log(1.0_wp/max(1.0e-20_wp,zion(4,m)**2))))
     zion(4,m)=zion(4,m)-(c0+c1*zion(4,m)+c2*zion(4,m)**2)/&
          (1.0+d1*zion(4,m)+d2*zion(4,m)**2+d3*zion(4,m)**3)
     if(zion(4,m)>umax)zion(4,m)=z4tmp
  enddo

!$omp parallel do private(m)
  do m=1,mi

! toroidal:  uniformly distributed in zeta
     zion(3,m)=zeta0+(zeta1-zeta0)*zion(3,m)

     zion(4,m)=zion0(4,m)*min(umax,zion(4,m))

! initial random weight
     zion(5,m)=2.0*w_initial*(zion(5,m)-0.5)*(1.0+cos(zion(2,m)))

! Maxwellian distribution in v_perp, <v_perp^2>=1.0
     zion(6,m)=max(1.0e-20_wp,min(umax*umax,-log(max(1.0e-20_wp,zion(6,m)))))
  enddo

! transform zion(1,:) to psi, zion(4,:) to rho_para, zion(6,:) to sqrt(mu) 
  vthi=sqrt(temion)*rho0*abs(qion)/aion
!$omp parallel do private(m)
  do m=1,mi
     zion0(1,m)=1.0/(1.0+zion(1,m)*cos(zion(2,m))) !B-field
     zion(1,m)=0.5*zion(1,m)*zion(1,m)
     zion(4,m)=vthi*zion(4,m)*aion/(qion*zion0(1,m))
     zion(6,m)=sqrt(aion*vthi*vthi*zion(6,m)/zion0(1,m))
  enddo

  cmratio=qion/aion !!XY

  if(iload==0)then
! uniform loading
!$omp parallel do private(m)
     do m=1,mi
        zion0(6,m)=1.0
     enddo

  else if(iload==1)then
! analytic temp profile
!$omp parallel do private(m,r,cost,b,upara,energy,en1,en2)
     do m=1,mi
!!XY initialize zion0(6,m)
        r=sqrt(2.0*zion(1,m))
        cost=cos(zion(2,m))
        b=1.0_wp/(1.0_wp+r*cost)
        upara=zion(4,m)*b*cmratio
        energy=0.5_wp*aion*upara*upara+zion(6,m)*zion(6,m)*b
        energy=energy/(aion*vthi*vthi)
        en1=exp((r-rc)*kappati) 
        en2=sqrt(en1*en1*en1)

        zion0(6,m)=en2*exp(-energy*(en1-1.0_wp))! for non-uniform temperature
!!!*exp((rc-r)*kappan) ! for non-uniform density

     enddo
  else
! true nonuniform for temperature profile, density profile by weight; need work here
!$omp parallel do private(m,r,i)
     do m=1,mi
        r=sqrt(2.0*zion(1,m))
        i=max(0,min(mpsi,int((r-a0)*delr+0.5)))
        zion(4,m)=zion(4,m)*sqrt(rtemi(i))
        zion(6,m)=zion(6,m)*sqrt(rtemi(i))
        zion0(6,m)=max(0.1_wp,min(10.0_wp,rden(i)))
     enddo
  endif
  
end subroutine loadi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine loade
  use global_parameters
  use particle_array
  use field_array
  use particle_tracking
  implicit none

  integer m,mtest,ierror
  real vthe,tsqrt,r,cost,b,upara,energy,eperp,cmratio

! load electron on top of first species of ion if mi=me
  vthe=sqrt(aelectron/aion)*(qion/qelectron) !assuming T_i=T_e
  tsqrt=1.0 !assuming T_i=T_e
  cmratio=qion/aion !! cmratio=qion/aion

  memax=me+100*ceiling(sqrt(real(me))) !electrons array upper bound
  allocate(zelectron(nparam,memax),zelectron0(nparam,memax),&
       markere(mgrid),densitye(0:1,mgrid),zelectron1(nparam,memax),&
       phisave(0:1,mgrid,2*nhybrid),phit(0:1,mgrid),pfluxe(0:mpsi),&
       rdteme(0:mpsi),pmarke(0:mpsi),zonale(0:mpsi),&
       data1de(0:mpsi,mpdata1d),STAT=mtest,SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
 
  if(mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate electron: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  phisave = 0. !---wj

  rdteme=0.0
  pfluxe=0.0
  zonale=0.0

  markere=markeri*real(mi)/real(me)
  pmarke=pmarki*real(mi)/real(me)
 
  if(irun>1)return

!! load only trapped electonr
  me=0
  do m=1,mi
     r=sqrt(2.0*zion(1,m))
     cost=cos(zion(2,m))
     b=1.0/(1.0+r*cost)
     upara=zion(4,m)*b*cmratio
     energy=0.5*aion*upara*upara+zion(6,m)*zion(6,m)*b
     eperp=zion(6,m)*zion(6,m)/(1.0-r)
     if(eperp>energy)then     ! keep trapped electrons
        me=me+1
        zelectron(1,me)=zion(1,m)
        zelectron(2,me)=zion(2,m)
        zelectron(3,me)=zion(3,m)
        zelectron(4,me)=zion(4,m)*vthe
        zelectron(5,me)=zion(5,m)
        zelectron(6,me)=zion(6,m)*tsqrt
        zelectron0(6,me)=zion0(6,m)
     endif
  enddo
  
end subroutine loade
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! restart from previous runs
subroutine restart_load
  use global_parameters
  implicit none

  integer iorestart,ierr
  character(len=9)cdum1

  if(mype==0)then
!!get restart directory info
     iorestart=345
     open(iorestart,file="FileExit.out",status="old")
     read(iorestart,"(A9,i1)")cdum1,FileExit
     read(iorestart,"(A9,i5)")cdum1,irest
     close(iorestart)         
     irest=irest-1
  endif
  call mpi_bcast(irest,1,mpi_integer,0,mpi_comm_world,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call restart_io("read")
  return

end subroutine restart_load

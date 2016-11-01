! Solve fields
subroutine field_solver
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer i,j,ij
  real(wp) r,theta,b,gperp2,g2_inv

! faked k_perp for benchmarking petsc
  gperp2=0.1 !k_perp rho0 value
  gperp2=(gperp2/rho0)**2 !k_perp^2 in basic unit
  g2_inv=1.0/gperp2

! Reduced MHD: phieff=0, k_perp rho_s <<1; GK Poisson Eq. becomes Laplace's Eq.
!$omp parallel do private(i,j,ij)
  do i=0,mpsi
     do j=1,mtheta(i)
        ij=igrid(i)+j
        
        phieff(1,ij)=fluidne(1,ij)*rho0*rho0
!        phieff(1,ij)=0.0

! Electrostatic potential; Could ignore ion contriution in reduced MHD.
        phi(1,ij)=(denion*densityi(1,ij)-fluidne(1,ij))*g2_inv
!        phi(1,ij)=-g2_inv*fluidne(1,ij)

! Inductive potential
        phiind(1,ij)=phieff(1,ij)-phi(1,ij)

! Electron parallel flow from Ampere's law
        fluidue(1,ij)=-gperp2*apara(1,ij)*rho0*rho0/betae+denion*currenti(1,ij)
     enddo
  enddo

! fluidue=delta_ue*B_0(R_0)/B_0(R), phi=[ B_0(R)/B_0(R_0)]^2 (delta_ne)
!$omp parallel do private(i,j,ij,r,theta,b)
  do i=0,mpsi
     r=a0+deltar*real(i)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        theta=deltat(i)*real(j)+zeta1*qtinv(i)
        b=1.0/(1.0+r*cos(theta)) !toroidal
        fluidue(1,ij) = fluidue(1,ij)/b
        phiind(1,ij) = phiind(1,ij)*b*b
     enddo
  enddo

end subroutine field_solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine field_gradient
  use global_parameters
  use field_array
  implicit none

  integer i
  real(wp) r

! field gradients
  if(magnetic==1)then
     call gradient(phiind,gradind)
     call gradient(fluidue,gradue)
     call gradient(apara,gradapara)
  endif

! grad phi_electrostatic
  call gradient(phi,gradphi)

! add (0,0) mode to d phi/d psi
  if(mode00==1)then
! solve zonal flow
     call zonal

!$omp parallel do private(i,r)
     do i=1,mpsi-1
        r=a0+deltar*real(i)
        gradphi(1,0:1,igrid(i):igrid(i)+mtheta(i))=gradphi(1,0:1,igrid(i):igrid(i)+mtheta(i))+&
             phip00(i)/r
     enddo
  endif
  
end subroutine field_gradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zonal
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer i,ismooth
  real den00(0:mpsi),r

! phi00=r*E_r, E_r(a0)=0. Trapezoid rule
  if(nhybrid==0)phip00=qion*zonali
  if(nhybrid>0)phip00=qion*zonali+qelectron*zonale
! (-0.0625 0.25 0.625 0.25 -0.0625) radial smoothing of (0,0) mode density phip00
  do ismooth=1,2
     den00(0)=phip00(0)
     den00(mpsi)=phip00(mpsi)
     den00(1)=phip00(3)
     den00(mpsi-1)=phip00(mpsi-3)
     den00(2:mpsi-2)=phip00(0:mpsi-4)+phip00(4:mpsi)
     den00(1:mpsi-1)=0.625*phip00(1:mpsi-1)+0.25*(phip00(0:mpsi-2)+phip00(2:mpsi))-&
          0.0625*den00(1:mpsi-1)
     phip00=den00
  enddo

  den00=phip00
  phip00=0.0

  do i=1,mpsi
     r=a0+deltar*real(i)
     phip00(i)=phip00(i-1)+0.5*deltar*((r-deltar)*den00(i-1)+r*den00(i))
  enddo

! subtract net momentum
!     phip00=phip00-sum(phip00)/real(mpsi+1)

! d phi/dr, in equilibrium unit
  do i=0,mpsi
     r=a0+deltar*real(i)
     phip00(i)=-phip00(i)/r
  enddo

! add FLR contribution using Pade approximation: b*<phi>=(1+b)*<n>
  phi00=den00*rho0*rho0 ! potential in equilibrium unit
  do i=1,mpsi-1
     phip00(i)=phip00(i)+0.5*(phi00(i+1)-phi00(i-1))/deltar
  enddo
  
! (0,0) mode potential store in phi00
  phi00=0.0
  do i=1,mpsi
     phi00(i)=phi00(i-1)+0.5*deltar*(phip00(i-1)+phip00(i))
  enddo
  if(mode00==0)phip00=0.0

end subroutine zonal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gradient(scalar,vector)
  use global_parameters
  use field_array
  implicit none

  integer i,ii,ij,j,k,icount,idest,isource,isendtag,irecvtag,ierror,ip,jt,&
      istatus(MPI_STATUS_SIZE)
  real(wp) diffr,difft(0:mpsi),diffz,r,drdp,pleft(mthetamax),pright(mthetamax),&
      sendl(mgrid),recvr(mgrid),sendr(mgrid),recvl(mgrid),sendrs(3,mgrid),&
      recvls(3,mgrid),q,delq
  real(wp),dimension(0:1,mgrid),intent(in) :: scalar
  real(wp),dimension(3,0:1,mgrid),intent(out) :: vector

! finite difference for e-field in equilibrium unit
  diffr=0.5_wp/deltar
  difft=0.5_wp/deltat
  diffz=0.5_wp/deltaz
!$omp parallel do private(i,j,k)
  do i=1,mgrid
     do j=0,1
        do k=1,3
          vector(k,j,i)=0.0
        enddo
     enddo
  enddo

! d_scalar/d_psi
!$omp parallel do private(i,j,r,drdp,ij)
  do i=1,mpsi-1
     r=a0+deltar*real(i)
     drdp=1.0/r
     do j=1,mtheta(i)
        ij=igrid(i)+j
        
        vector(1,1,ij)=drdp*diffr*((1.0-wtp1(1,ij))*scalar(1,jtp1(1,ij))+&
            wtp1(1,ij)*scalar(1,jtp1(1,ij)+1)-&
            ((1.0-wtp1(2,ij))*scalar(1,jtp1(2,ij))+wtp1(2,ij)*scalar(1,jtp1(2,ij)+1)))
     enddo
   enddo

! d_scalar/d_theta
!$omp parallel do private(i,j,ij,jt)
  do i=1,mpsi-1
     do j=1,mtheta(i)
        ij=igrid(i)+j
        jt=j+1-mtheta(i)*(j/mtheta(i))
        vector(2,1,ij)=difft(i)*(scalar(1,igrid(i)+jt)-scalar(1,igrid(i)+j-1))
     enddo
  enddo
  
! send scalar to right and receive from left
  sendr=scalar(1,:)
  recvl=0.0
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
      recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
  
! send scalar to left and receive from right
  sendl=scalar(1,:)
  recvr=0.0
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
      recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

! unpack scalar_boundary and calculate E_zeta at boundaries toroidal
!$omp parallel do private(i,j,ii,jt,ij,pleft,pright)
  do i=1,mpsi-1
     ii=igrid(i)
     jt=mtheta(i)
     if(myrank_toroidal==0)then !down-shift for zeta=0
        pleft(1:jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
        pright(1:jt)=recvr(ii+1:ii+jt)
     elseif(myrank_toroidal==mtoroidal-1)then !up-shift for zeta=2*pi
        pright(1:jt)=cshift(recvr(ii+1:ii+jt),itran(i))
        pleft(1:jt)=recvl(ii+1:ii+jt)
     else
        pleft(1:jt)=recvl(ii+1:ii+jt)
        pright(1:jt)=recvr(ii+1:ii+jt)
     endif

! d_scalar/d_zeta
     do j=1,mtheta(i)
        ij=igrid(i)+j
        vector(3,1,ij)=(pright(j)-pleft(j))*diffz
     enddo
  enddo

! adjust the difference between safety factor q and qtinv for fieldline coordinate
!omp parallel do private(i,j,r,q,delq,ij)
  do i=1,mpsi-1
     r=a0+deltar*real(i)
     q=q0+q1*r/aminor+q2*r*r/(aminor*aminor)
     delq=(1.0/q-qtinv(i))

     do j=1,mtheta(i)
        ij=igrid(i)+j
        vector(3,:,ij)=vector(3,:,ij)+delq*vector(2,:,ij)
     enddo
  enddo

  call periodicity(vector(1,:,:))
  call periodicity(vector(2,:,:))
  call periodicity(vector(3,:,:))

end subroutine gradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

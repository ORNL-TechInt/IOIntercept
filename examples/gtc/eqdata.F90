subroutine eqdata
  use global_parameters
  use equilibrium
  implicit none

  integer i,icount,ierror

! # of grids for spline, lsp and lst should be same as those in spdata.dat and profile.dat
  lsp=81
  lst=81

! equilibrium and profile are represented using 2D-spline on (psi, theta) or 1D-spline on (psi) 
! poloidal flux=[0,psiw], spdpsi=pw/(lsp-1); theta=[0,2pi), spdtheta=2pi/lst
  spdtheta=2.0*pi/real(lst)

! allocate memonry for equilibrium spline array
  allocate (stpp(lsp),mipp(lsp),mapp(lsp),nepp(3,lsp),tepp(3,lsp),tipp(3,lsp),zepp(3,lsp),&
       ropp(3,lsp),erpp(3,lsp),qpsi(3,lsp),gpsi(3,lsp),ppsi(3,lsp),rpsi(3,lsp),torpsi(3,lsp),&
       b(9,lsp,lst),x(9,lsp,lst),z(9,lsp,lst),g(9,lsp,lst),rd(9,lsp,lst),nu(9,lsp,lst),&
       dl(9,lsp,lst),spcos(3,lst),spsin(3,lst),rgpsi(3,lsp),cpsi(3,lsp),psirg(3,lsp),psitor(3,lsp),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj

  stpp = 0. !---wj
  mipp = 0. !---wj
  mapp = 0. !---wj
  nepp = 0. !---wj
  tepp = 0. !---wj
  tipp = 0. !---wj
  zepp = 0. !---wj
  ropp = 0. !---wj
  erpp = 0. !---wj
  qpsi = 0. !---wj
  gpsi = 0. !---wj
  ppsi = 0. !---wj
  rpsi = 0. !---wj
  torpsi = 0. !---wj

  spcos = 0. !---wj
  spsin = 0. !---wj
  rgpsi = 0. !---wj
  cpsi = 0. !---wj
  psirg = 0. !---wj
  psitor = 0. !---wj

  g = 0.0 !---wj
  rd = 0.0 !---wj
  nu = 0.0 !---wj
  dl = 0.0 !---wj

  if(mype==0)then

     if(numereq==1)then
! use EFIT & TRANSP data
! EFIT spdata used in GTC: 2D spline b,x,z,rd; 1D spline qpsi,gpsi,torpsi
        call spdata

! TRANSP prodata used in GTC: tipp,tepp,nepp,zeff,ropp,erpp
        call prodata

     else

! Construct analytic equilibrium and profile assuming high aspect-ratio, concentric cross-section
        call analyticeq

     endif

! End of equilibrium and profile construction
! Now, define a radial coordinate for radial grid: d(rpsi)/d(rgpsi)=sqrt(T_i)
! Using rgpsi, grid cell size for gather-scatter operations is propotioal to local ion gyroradius
     rgpsi(1,1)=0.0
     do i=2,lsp
        rgpsi(1,i)=rgpsi(1,i-1)+(rpsi(1,i)-rpsi(1,i-1))/sqrt(tipp(1,i))
     enddo
     call construct_spline1(lsp,spdpsi,rgpsi)
!     write(gtcout,*)'rpsi=',rpsi(1,:),'rgpsi=',rgpsi(1,:)

! inverse spline fit: psitor,psirg
! spline fit cell size     
     spdtor=torpsi(1,lsp)/real(lsp-1)
     spdrg=rgpsi(1,lsp)/real(lsp-1)
     call invert_spline0(1,lsp,spdpsi,spdtor,torpsi,psitor)
     call invert_spline1(lsp,spdpsi,spdrg,rgpsi,psirg)
  endif

! broadcast eqquilibrium and proile data to all MPI processes
! radial spline fit of plasma profiles
  icount=3*lsp
  call MPI_BCAST(nepp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(tepp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(tipp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(zepp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(ropp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(erpp,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

! radial spline fit of geometry
  call MPI_BCAST(qpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(gpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(cpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(rpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(rgpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(torpsi,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(psirg,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(psitor,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

! poloidal spline fit of sin and cos functions
  icount=3*lst
  call MPI_BCAST(spcos,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(spsin,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

! radial and poloildal spline fit of B field and (X,Z) coordinates
  icount=9*lsp*lst
  call MPI_BCAST(b,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(x,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(z,icount,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

end subroutine eqdata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine spdata
  use global_parameters
  use equilibrium
  implicit none

  integer i,j,isp,lemax,lrmax,krip,nripa,nripb,nrip,iformat,lspdum,lstdum
  real ped,rmaj,d0,brip,wrip,xrip
  character(len=11) cdum

  iformat=2 !1: for equal arc and Boozer; 2 : for Boozer only

! read equilibrium spline parameters
  isp=110
  open(isp,file='spdata.dat',status='old')
  read(isp,'(a11)') cdum
  read(isp,102) lspdum,lstdum,lemax,lrmax
  if(iformat==1)read(isp,101) psiw,ped,rmaj !psiw: poloidal flux at wall
  if(iformat==2)read(isp,101) psiw,ped
  write(gtcout,*) "eqdata=",cdum,iformat,lspdum,lstdum,lemax,lrmax,psiw,ped
  if(lspdum /= lsp .or. lstdum /= lst)write(gtcout,*)'WARNING: WRONG spline',lspdum,lstdum,lsp,lst

  spdpsi=psiw/real(lsp-1)

! read spline array
  if(iformat==1)then
! for equal arc length and Boozer coordinates
     do i=1,lsp
!        read(isp,101)b(1,i,:) !B-field, b(1,0,0)=1
!        read(isp,101)b(2,i,:)
!        read(isp,101)b(3,i,:)
!        read(isp,101)b(4,i,:)
!        read(isp,101)b(5,i,:)
!        read(isp,101)b(6,i,:)
!        read(isp,101)b(7,i,:)
!        read(isp,101)b(8,i,:)
!        read(isp,101)b(9,i,:)
!        read(isp,101)x(1,i,:) !X-coordinate, x(0,0)=1
!        read(isp,101)x(2,i,:)
!        read(isp,101)x(3,i,:)
!        read(isp,101)x(4,i,:)
!        read(isp,101)x(5,i,:)
!        read(isp,101)x(6,i,:)
!        read(isp,101)x(7,i,:)
!        read(isp,101)x(8,i,:)
!        read(isp,101)x(9,i,:)
!        read(isp,101)z(1,i,:) !Z-coordinate
!        read(isp,101)z(2,i,:)
!        read(isp,101)z(3,i,:)
!        read(isp,101)z(4,i,:)
!        read(isp,101)z(5,i,:)
!        read(isp,101)z(6,i,:)
!        read(isp,101)z(7,i,:)
!        read(isp,101)z(8,i,:)
!        read(isp,101)z(9,i,:)
!        read(isp,101)g(1,i,:) !Gicobian; use (gq+I)/b^2 in GTC
!        read(isp,101)g(2,i,:)
!        read(isp,101)g(3,i,:)
!        read(isp,101)g(4,i,:)
!        read(isp,101)g(5,i,:)
!        read(isp,101)g(6,i,:)
!        read(isp,101)g(7,i,:)
!        read(isp,101)g(8,i,:)
!        read(isp,101)g(9,i,:)
!        read(isp,101)rd(1,i,:) !I
!        read(isp,101)rd(2,i,:)
!        read(isp,101)rd(3,i,:)
!        read(isp,101)rd(4,i,:)
!        read(isp,101)rd(5,i,:)
!        read(isp,101)rd(6,i,:)
!        read(isp,101)rd(7,i,:)
!        read(isp,101)rd(8,i,:)
!        read(isp,101)rd(9,i,:)
!        read(isp,101)nu(1,i,:) !zeta-phi
!        read(isp,101)nu(2,i,:)
!        read(isp,101)nu(3,i,:)
!        read(isp,101)nu(4,i,:)
!        read(isp,101)nu(5,i,:)
!        read(isp,101)nu(6,i,:)
!        read(isp,101)nu(7,i,:)
!        read(isp,101)nu(8,i,:)
!        read(isp,101)nu(9,i,:)
!        read(isp,101)dl(1,i,:) !del
!        read(isp,101)dl(2,i,:)
!        read(isp,101)dl(3,i,:)
!        read(isp,101)dl(4,i,:)
!        read(isp,101)dl(5,i,:)
!        read(isp,101)dl(6,i,:)
!        read(isp,101)dl(7,i,:)
!        read(isp,101)dl(8,i,:)
!        read(isp,101)dl(9,i,:)
!        read(isp,101)qpsi(:,i) !q
!        read(isp,101)gpsi(:,i) !g-current
!        read(isp,101)ppsi(:,i) !pressure; not used in GTC
!        read(isp,101)rpsi(:,i) !minor radius; use (x(:,0)-x(0,0)) in GTC
!        read(isp,101)torpsi(:,i) !toroidal flux; not used in GTC
     enddo
     read(isp,102) krip,nripa,nripb
!        do i=1,lsp
!           read(isp,101)ha(1,i,:),ha(2,i,:),ha(3,i,:),ha(4,i,:),ha(5,i,:),ha(6,i,:),ha(7,i,:),&
!            ha(8,i,:),ha(9,i,:)
!           read(isp,101)hb(1,i,:),hb(2,i,:),hb(3,i,:),hb(4,i,:),hb(5,i,:),hb(6,i,:),hb(7,i,:),&
!            hb(8,i,:),hb(9,i,:)
!        enddo
     write(gtcout,*)"ripple=",krip,nripa,nripb,rmaj
  else
! for Boozer coordinates only
     do i=1,lsp
        read(isp,101)b(1,i,:) !B-field
        read(isp,101)b(2,i,:)
        read(isp,101)b(3,i,:)
        read(isp,101)b(4,i,:)
        read(isp,101)b(5,i,:)
        read(isp,101)b(6,i,:)
        read(isp,101)b(7,i,:)
        read(isp,101)b(8,i,:)
        read(isp,101)b(9,i,:)
        read(isp,101)x(1,i,:) !X-coordinate
        read(isp,101)x(2,i,:)
        read(isp,101)x(3,i,:)
        read(isp,101)x(4,i,:)
        read(isp,101)x(5,i,:)
        read(isp,101)x(6,i,:)
        read(isp,101)x(7,i,:)
        read(isp,101)x(8,i,:)
        read(isp,101)x(9,i,:)
        read(isp,101)z(1,i,:) !Z-coordinate
        read(isp,101)z(2,i,:)
        read(isp,101)z(3,i,:)
        read(isp,101)z(4,i,:)
        read(isp,101)z(5,i,:)
        read(isp,101)z(6,i,:)
        read(isp,101)z(7,i,:)
        read(isp,101)z(8,i,:)
        read(isp,101)z(9,i,:)
        read(isp,101)g(1,i,:) !Gicobian
        read(isp,101)g(2,i,:)
        read(isp,101)g(3,i,:)
        read(isp,101)g(4,i,:)
        read(isp,101)g(5,i,:)
        read(isp,101)g(6,i,:)
        read(isp,101)g(7,i,:)
        read(isp,101)g(8,i,:)
        read(isp,101)g(9,i,:)
        read(isp,101)qpsi(:,i) !q
        read(isp,101)gpsi(:,i) !g-current
        read(isp,101)rd(1:3,i,1) !i-current
        read(isp,101)ppsi(:,i) !pressure
        read(isp,101)rpsi(:,i) !minor radius
        read(isp,101)torpsi(:,i) !toroidal flux
        rd(1,i,:)=rd(1,i,1)
        rd(2,i,:)=rd(2,i,1)
        rd(3,i,:)=rd(3,i,1)
     enddo
     rd(4:9,:,:)=0.0
     nu=0.0
     dl=0.0
     read(isp,102) krip,nrip
     read(isp,101) rmaj,d0,brip
     read(isp,101) wrip,xrip
!        do i=1,lsp
!           read(isp,101)ha(1,i,:),ha(2,i,:),ha(3,i,:),ha(4,i,:),ha(5,i,:),ha(6,i,:),ha(7,i,:),&
!            ha(8,i,:),ha(9,i,:)
!           read(isp,101)hb(1,i,:),hb(2,i,:),hb(3,i,:),hb(4,i,:),hb(5,i,:),hb(6,i,:),hb(7,i,:),&
!            hb(8,i,:),hb(9,i,:)
!        enddo
     write(gtcout,*)"ripple=",krip,nrip,rmaj,d0,brip,wrip,xrip
  endif
  
  close(isp)
  write(gtcout,*)"X(0,0)=",x(1,1,1),"B(0,0)=",b(1,1,1)

101 format(1p4e18.10)
102 format(6i4)
  
end subroutine spdata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine prodata
  use global_parameters
  use equilibrium
  implicit none
  
  integer, parameter :: lpp=11,npp=22 !lpp=# of quantities; npp=2+# of radial grids of profile data
  integer i,j,ioprof
  real ppdata(lpp,npp),px,dpx
  character(len=11) ppname(lpp)

! read equilibrium spline parameters
  ioprof=110

! read NSTX plasma profile data
  open(ioprof,file='profile.dat',status='old')
  read(ioprof,*)ppname 
  write(gtcout,*)"NSTX plasma profile data"
  write(gtcout,*)ppname
  read(ioprof,*)ppdata(:,2:npp-1) 
  close(ioprof)
  
! extend TRANSP data to magnetic axis
  ppdata(1,1)=0.0
  ppdata(2,1)=0.0
  ppdata(3,1)=0.0
  ppdata(4,1)=ppdata(4,2)
  ppdata(5,1)=ppdata(4,2)
  ppdata(6,1)=ppdata(6,2)
  ppdata(7,1)=ppdata(7,2)
  ppdata(8,1)=ppdata(8,2)
  ppdata(9,1)=ppdata(9,2)
  ppdata(10,1)=ppdata(10,2)
  ppdata(11,1)=0.0

! extend TRANSP data to wall
  ppdata(1,npp)=ppdata(1,21)+0.5*(ppdata(1,21)-ppdata(1,20))
  ppdata(2,npp)=ppdata(2,21)+0.5*(ppdata(2,21)-ppdata(2,20))
  ppdata(3,npp)=ppdata(3,21)+0.5*(ppdata(3,21)-ppdata(3,20))
  ppdata(4,npp)=ppdata(4,21)+0.5*(ppdata(4,21)-ppdata(4,20))
  ppdata(5,npp)=ppdata(3,npp)+ppdata(4,npp)
  ppdata(6,npp)=ppdata(6,21)+0.5*(ppdata(6,21)-ppdata(6,20))
  ppdata(7,npp)=ppdata(7,21)+0.5*(ppdata(7,21)-ppdata(7,20))
  ppdata(8,npp)=ppdata(8,21)+0.5*(ppdata(8,21)-ppdata(8,20))
  ppdata(9,npp)=ppdata(9,21)+0.5*(ppdata(9,21)-ppdata(9,20))
  ppdata(10,npp)=ppdata(10,21)+0.5*(ppdata(10,21)-ppdata(10,20))
  ppdata(11,npp)=ppdata(11,21)+0.5*(ppdata(11,21)-ppdata(11,20))     

! spline grid data using linear interpolation on poloidal flux
  do i=2,lsp-1
     px=real(i-1)*ppdata(1,npp)/real(lsp-1)
     j=0
1    continue
     j=j+1
     if(px>ppdata(1,j)) goto 1
     dpx=(ppdata(1,j)-px)/(ppdata(1,j)-ppdata(1,j-1))
     stpp(i)=(1.0-dpx)*ppdata(2,j)+dpx*ppdata(2,j-1)
     mipp(i)=(1.0-dpx)*ppdata(3,j)+dpx*ppdata(3,j-1)
     mapp(i)=(1.0-dpx)*ppdata(5,j)+dpx*ppdata(5,j-1)
     tepp(1,i)=(1.0-dpx)*ppdata(6,j)+dpx*ppdata(6,j-1)
     nepp(1,i)=(1.0-dpx)*ppdata(7,j)+dpx*ppdata(7,j-1)
     tipp(1,i)=(1.0-dpx)*ppdata(8,j)+dpx*ppdata(8,j-1)
     zepp(1,i)=(1.0-dpx)*ppdata(9,j)+dpx*ppdata(9,j-1)
     ropp(1,i)=(1.0-dpx)*ppdata(10,j)+dpx*ppdata(10,j-1)
     erpp(1,i)=(1.0-dpx)*ppdata(11,j)+dpx*ppdata(11,j-1)
  enddo

! magnetix axis
  stpp(1)=ppdata(2,1)
  mipp(1)=ppdata(3,1)
  mapp(1)=ppdata(5,1)
  tepp(1,1)=ppdata(6,1)
  nepp(1,1)=ppdata(7,1)
  tipp(1,1)=ppdata(8,1)
  zepp(1,1)=ppdata(9,1)
  ropp(1,1)=ppdata(10,1)
  erpp(1,1)=ppdata(11,1)

! wall
  stpp(lsp)=ppdata(2,npp)
  mipp(lsp)=ppdata(3,npp)
  mapp(lsp)=ppdata(5,npp)
  tepp(1,lsp)=ppdata(6,npp)
  nepp(1,lsp)=ppdata(7,npp)
  tipp(1,lsp)=ppdata(8,npp)
  zepp(1,lsp)=ppdata(9,npp)
  ropp(1,lsp)=ppdata(10,npp)
  erpp(1,lsp)=ppdata(11,npp)  
! end of NSTX profile data

! spline fit plasma profile data
  call construct_spline0(0,lsp,spdpsi,tepp)
  call construct_spline0(0,lsp,spdpsi,nepp)
  call construct_spline0(0,lsp,spdpsi,tipp)
  call construct_spline0(0,lsp,spdpsi,zepp)
  call construct_spline0(0,lsp,spdpsi,ropp)
  call construct_spline0(0,lsp,spdpsi,erpp)

101 format(1p4e18.10)
102 format(6i4)

end subroutine prodata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Construct analytic equilibrium and profile assuming high aspect-ratio, concentric cross-section
subroutine analyticeq
  use global_parameters
  use equilibrium
  implicit none

  integer i,j
  real dum,spdum,spdum_inv

! spline fit for cos and sin functions in domain [0,2pi]: spcos,spsin
  spdum=2.0*pi/real(lst-1)
  do i=1,lst
     dum=spdum*real(i-1)
     spcos(1,i)=cos(dum)
  enddo
  call construct_spline0(0,lsp,spdum,spcos)

! sin function obtained from cos function using sin(theta+pi/2)=cos(theta)
! requires lst=4*n+1, n=positive integer
  if(mod(lst-1,4)>0)write(gtcout,*)"WARNING: WRONG LST #. lst=", lst
  do i=1,lst/4
     j=i+3*lst/4
     spsin(1:3,i)=spcos(1:3,j)
  enddo
  do i=lst/4+1,lst
     j=i-lst/4
     spsin(1:3,i)=spcos(1:3,j)
  enddo
  
! 1D spline data functions: qpsi,gpsi,torpsi,tipp,tepp,nepp,zeff,ropp,erpp
! Use Boozer coordinates; Asuming a parabolic q profile q=q_0+q_2(r/a)^2
! Poloidal flux at wall and spline grids on poloidal flux
  spdum_inv=0.5*aminor*aminor/q2
  psiw=spdum_inv*log(q2/q0+1.0)
  spdpsi=psiw/real(lsp-1)
        
! toroidal flux torpsi=0.5*r*r, and q factor
  spdum=2.0*q2/(aminor*aminor)
  do i=1,lsp
     dum=spdpsi*real(i-1)
     torpsi(1,i)=q0*spdum_inv*(exp(spdum*dum)-1.0)
     qpsi(1,i)=q0+spdum*torpsi(1,i)
     rpsi(1,i)=sqrt(2.0*torpsi(1,i))
  enddo
  
! g-current is 1 to order of epsilon
  gpsi(1,:)=1.0
        
! c-current is 0 to order of epsilon
  cpsi(1,:)=1.0
        
! density and temperature profiles
  nepp(1,:)=1.0        
  tepp(1,:)=1.0
  tipp(1,:)=1.0

! 1D spline fit: spline0 for y'=0 at x=0
  call construct_spline0(0,lsp,spdpsi,qpsi)
  call construct_spline0(0,lsp,spdpsi,gpsi)
  call construct_spline0(0,lsp,spdpsi,cpsi)
  call construct_spline0(0,lsp,spdpsi,nepp)
  call construct_spline0(0,lsp,spdpsi,tepp)
  call construct_spline0(0,lsp,spdpsi,tipp)

! 1D spline fit: spline0 for y as a linear function of x near x=0
  call construct_spline0(1,lsp,spdpsi,torpsi)
! near r=0, r=sqrt(2*q_0*psi), expand r in term of sqrt(psi) in first spline cell only
! spline2 for y as a linear function of sqrt(x) near x=0
  call construct_spline1(lsp,spdpsi,rpsi)

! 2D spline functions: b,x,z,rd
! x=1+rcos(theta)
  do i=1,lsp
     do j=1,lst
!        x(:,i,j)=rpsi(:,i)*spcos(:,j) !matrix multiplication not allowed in debug mode
        x(1,i,j)=rpsi(1,i)*spcos(1,j)
        x(2,i,j)=rpsi(2,i)*spcos(1,j)
        x(3,i,j)=rpsi(3,i)*spcos(1,j)
        x(4,i,j)=rpsi(1,i)*spcos(2,j)
        x(5,i,j)=rpsi(2,i)*spcos(2,j)
        x(6,i,j)=rpsi(3,i)*spcos(2,j)
        x(7,i,j)=rpsi(1,i)*spcos(3,j)
        x(8,i,j)=rpsi(2,i)*spcos(3,j)
        x(9,i,j)=rpsi(3,i)*spcos(3,j)
     enddo
  enddo

! z=rsin(theta)
  do i=1,lsp
     do j=1,lst
!        z(:,i,j)=rpsi(:,i)*spsin(:,j)
        z(1,i,j)=rpsi(1,i)*spsin(1,j)
        z(2,i,j)=rpsi(2,i)*spsin(1,j)
        z(3,i,j)=rpsi(3,i)*spsin(1,j)
        z(4,i,j)=rpsi(1,i)*spsin(2,j)
        z(5,i,j)=rpsi(2,i)*spsin(2,j)
        z(6,i,j)=rpsi(3,i)*spsin(2,j)
        z(7,i,j)=rpsi(1,i)*spsin(3,j)
        z(8,i,j)=rpsi(2,i)*spsin(3,j)
        z(9,i,j)=rpsi(3,i)*spsin(3,j)
     enddo
  enddo

! B=1/(1+rcos(theta))=1-rcos(theta) to order of epsilon
  b=x
  b(1,:,:)=1.0-b(1,:,:)
  x(1,:,:)=x(1,:,:)+1.0
        
! I-current is zero to order of epsilon
  rd=0.0

end subroutine analyticeq
!!!!!!!!!!!!!

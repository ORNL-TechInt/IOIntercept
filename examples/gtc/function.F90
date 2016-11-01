!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! function of plasma profile
! electron density
real function edensity(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline0

  edensity=spline0(pdum,lsp,spdpsi,nepp)
end function edensity

! d_densitye/d_psi
real function edensity_dp(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: dspline0

  edensity_dp=dspline0(pdum,lsp,spdpsi,nepp)
end function edensity_dp

! electron temperature
real function etemperature(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline0

  etemperature=spline0(pdum,lsp,spdpsi,tepp)
end function etemperature

! d_etemperature/d_psi
real function etemperature_dp(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: dspline0

  etemperature_dp=dspline0(pdum,lsp,spdpsi,tepp)
end function etemperature_dp

! ion temperature
real function itemperature(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline0

  itemperature=spline0(pdum,lsp,spdpsi,tipp)
end function itemperature

! d_itemperature/d_psi
real function itemperature_dp(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: dspline0

  itemperature_dp=dspline0(pdum,lsp,spdpsi,tipp)
end function itemperature_dp

! Z_eff
real function zeff(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline0

  zeff=spline0(pdum,lsp,spdpsi,zepp)
end function zeff

! toroidal rotation
real function trotation(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline0

  trotation=spline0(pdum,lsp,spdpsi,ropp)
end function trotation

! spline cos function
real function splcos(pdum,delx)
  use equilibrium  
  implicit none
  real pdum,delx
  real,external :: spline0

  splcos=spline0(pdum,lsp,delx,spcos)
end function splcos

! spline sin function
real function splsin(pdum,delx)
  use equilibrium  
  implicit none
  real pdum,delx
  real,external :: spline0

  splsin=spline0(pdum,lsp,delx,spsin)
end function splsin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function from equilibrium data
! safety factor
real function qfactor(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline0

  qfactor=spline0(pdum,lsp,spdpsi,qpsi)
end function qfactor

! g and gp
real function gcurrent(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline0

  gcurrent=spline0(pdum,lsp,spdpsi,gpsi)
end function gcurrent

real function gcurrent_dp(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: dspline0

  gcurrent_dp=dspline0(pdum,lsp,spdpsi,gpsi)
end function gcurrent_dp

! pressure
real function pressure(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline0

  pressure=spline0(pdum,lsp,spdpsi,ppsi)
end function pressure

! toroidal flux
real function torflux(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline0

  torflux=spline0(pdum,lsp,spdpsi,torpsi)
end function torflux

! minor radius
real function radius(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline1

  radius=spline1(pdum,lsp,spdpsi,rpsi)
end function radius

! radial grid
real function rgrid(pdum)
  use equilibrium  
  implicit none
  real pdum
  real,external :: spline1

  rgrid=spline1(pdum,lsp,spdpsi,rgpsi)
end function rgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! magnetic field amplitude
real function bfield(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  bfield=spline2d(0,pdum,tdum,lsp,lst,spdpsi,spdtheta,b)
end function bfield

! db_field/dpsi
real function dbdp(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dbdp=spline2d(1,pdum,tdum,lsp,lst,spdpsi,spdtheta,b)
end function dbdp

! db_field/dtheta
real function dbdt(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  dbdt=spline2d(2,pdum,tdum,lsp,lst,spdpsi,spdtheta,b)
end function dbdt

! transform from (psi,theta) to (X,Z)
real function xcoor(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  xcoor=spline2d(0,pdum,tdum,lsp,lst,spdpsi,spdtheta,x)
end function xcoor

! transform from (psi,theta) to (X,Z)
real function zcoor(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  zcoor=spline2d(0,pdum,tdum,lsp,lst,spdpsi,spdtheta,z)
end function zcoor

! Jacobian
real function jacobian(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  jacobian=spline2d(0,pdum,tdum,lsp,lst,spdpsi,spdtheta,g)
end function jacobian

! toroidal current I
real function icurrent(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  icurrent=spline2d(0,pdum,tdum,lsp,lst,spdpsi,spdtheta,rd)
end function icurrent

real function didp(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  didp=spline2d(1,pdum,tdum,lsp,lst,spdpsi,spdtheta,rd)
end function didp

! difference between magnetic angle zeta and cylindrical angle phi
real function zeta2phi(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  zeta2phi=spline2d(0,pdum,tdum,lsp,lst,spdpsi,spdtheta,nu)
end function zeta2phi

! delta in B-field contravariant representation
real function delb(pdum,tdum)
  use equilibrium
  implicit none
  real pdum,tdum
  real,external :: spline2d

  delb=spline2d(0,pdum,tdum,lsp,lst,spdpsi,spdtheta,dl)
end function delb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1D spline for radial profiles
real function spline0(pdum,lsp,delx,y)
  implicit none
  integer lsp,i
  real pdum,y(3,lsp),delx,dpx

  i=max(1,min(lsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
  spline0=y(1,i)+dpx*y(2,i)+dpx*dpx*y(3,i)

end function spline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! derivative of 1D spline function
real function dspline0(pdum,lsp,delx,y)
  implicit none
  integer lsp,i
  real pdum,y(3,lsp),delx,dpx

  i=max(1,min(lsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
  dspline0=y(2,i)+2.0*dpx*y(3,i)

end function dspline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1D spline with first point being linear function y=sqrt(x)
real function spline1(pdum,lsp,delx,y)
  implicit none
  integer lsp,i
  real pdum,y(3,lsp),delx,dpx,sdum

  i=max(1,min(lsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
! expand y(x) using sprt(x) near x=0   
  sdum=0.5-0.5*sign(1.0,i-1.5) !sdum=0 for i>1; sdum=1 for i=1
  dpx=(1.0-sdum)*dpx + sdum*sqrt(max(1.0e-10,dpx))
  spline1=y(1,i)+dpx*y(2,i)+dpx*dpx*y(3,i)

end function spline1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2D spline on poloidal plane
real function spline2d(iflag,pdum,tdum,lsp,lst,delx,dely,f)
  implicit none
  integer iflag,lsp,lst,i,j
  real pdum,tdum,f(9,lsp,lst),dpx,dp2,dtx,dt2,sdum,delx,dely

  i=max(1,min(lsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
! remove sigularity near psi=1 due to functional form of sqrt(psi)
! expand f(x,y) using sprt(x) near x=0; i.e., subroutine spline1 in 1D case   
  sdum=0.5-0.5*sign(1.0,i-1.5) !sdum=0 for i>1; sdum=1 for i=1
  dpx=(1.0-sdum)*dpx + sdum*sqrt(max(1.0e-10,dpx))
  dp2=dpx*dpx

  j=max(1,min(lst,ceiling(tdum/dely)))
  dtx=tdum-dely*real(j-1)
  dt2=dtx*dtx
  
  if(iflag==0)then !2D spline value
     spline2d=f(1,i,j)    +f(2,i,j)*dpx    +f(3,i,j)*dp2 &
          +f(4,i,j)*dtx+f(5,i,j)*dtx*dpx+f(6,i,j)*dtx*dp2 &
          +f(7,i,j)*dt2+f(8,i,j)*dt2*dpx+f(9,i,j)*dt2*dp2
  elseif(iflag==1)then !derivative with respect to psi
     spline2d=f(2,i,j)+f(5,i,j)*dtx+f(8,i,j)*dt2+2.0*dpx*(f(3,i,j)+f(6,i,j)*dtx+f(9,i,j)*dt2)
  else !derivative with respect to theta
     spline2d=f(4,i,j)+f(5,i,j)*dpx+f(6,i,j)*dp2+2.0*dtx*(f(7,i,j)+f(8,i,j)*dpx+f(9,i,j)*dp2)
  endif
end function spline2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! construct 1D spline with on x=[0,xmax], grid x_i = (i-1)*delx, delx=xmax/(nsp-1)
! in domain i, y = y(1,i) + y(2,i)*delx + y(3,i)*delx**2
subroutine construct_spline0(iflag,nsp,delx,y)
  implicit none
  integer i,nsp,iflag
  real delx,y(3,nsp)

! first point
  if(iflag==0)then
! iflag=0: first point derivative being zero y'=0
     y(2,1)=0.0
     y(3,1)=(y(1,2)-y(1,1)-delx*y(2,1))/(delx*delx)
     elseif(iflag==1)then
! iflag=1: first point being linear function y=x
        y(1,1)=0.0
        y(2,1)=y(1,2)/delx
        y(3,1)=0.0
     elseif(iflag==2)then
! iflag=2: first point being quadratic function y=x*x
        y(1,1)=0.0
        y(2,1)=0.0
        y(3,1)=y(1,2)/(delx*delx)
     endif

! other points y(2,i) by derivative at grid i; y(3,i) by value at grid (i+1)
  do i=2,nsp-1
     y(2,i)=y(2,i-1)+2.0*delx*y(3,i-1)
     y(3,i)=(y(1,i+1)-y(1,i)-delx*y(2,i))/(delx*delx)
  enddo

! last point is not used;
  y(2,nsp)=0.0
  y(3,nsp)=0.0

end subroutine construct_spline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! inversion of spline1 y(x) to x(y)
subroutine invert_spline0(iflag,nsp,delx,dely,y,x)
  implicit none
  integer i,nsp,j,iflag
  real delx,dely,y(3,nsp),x(3,nsp),ydum,y0,y1,y2

! first point
  x(1,1)=0.0

! other points
  do i=2,nsp-1
! y grid
     ydum=dely*real(i-1)
! search x grid for ydum     
     j=1
     do while (ydum>y(1,j+1))
        j=j+1
     enddo

! x(1,i)=j grid location + distance from j grid     
     y0=y(1,j)
     y1=y(2,j)
     y2=y(3,j)
     x(1,i)=delx*real(j-1)+(sqrt(y1*y1+4.0*y2*(ydum-y0))-y1)/(2.0*y2)
  enddo

! last point
  x(1,nsp)=delx*real(nsp-1)

! spline fit x function
  if(iflag==1)call construct_spline0(1,nsp,dely,x)

end subroutine invert_spline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! spline1 for cases with first point being linear function y=sqrt(x)
subroutine construct_spline1(nsp,delx,y)
  implicit none
  integer i,nsp
  real delx,y(3,nsp)

! first point
  y(1,1)=0.0
  y(2,1)=y(1,2)/sqrt(delx)
  y(3,1)=0.0

! second point
  y(2,2)=0.5*y(2,1)/sqrt(delx)
  y(3,2)=(y(1,3)-y(1,2)-delx*y(2,2))/(delx*delx)

! other points f(2,i) by derivative at grid i; f(3,i) by value at grid (i+1)
  do i=3,nsp-1
     y(2,i)=y(2,i-1)+2.0*delx*y(3,i-1)
     y(3,i)=(y(1,i+1)-y(1,i)-delx*y(2,i))/(delx*delx)
  enddo

! last point is not used;
  y(2,nsp)=0.0
  y(3,nsp)=0.0

end subroutine construct_spline1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! inversion of spline2 y(x) to x(y)
subroutine invert_spline1(nsp,delx,dely,y,x)
  implicit none
  integer i,nsp,j
  real delx,dely,y(3,nsp),x(3,nsp),ydum,y0,y1,y2

! first point
  x(1,1)=0.0

! other points
  do i=2,nsp-1
! y grid
     ydum=dely*real(i-1)
! search x grid for ydum     
     j=1
     do while (ydum>y(1,j+1))
        j=j+1
     enddo

! x(1,i)=j grid location + distance from j grid     
     if(j==1)then
! first point y=y1*sqrt(x)
        x(1,i)=(ydum/y(2,1))**2
     else
! other points as usual        
        y0=y(1,j)
        y1=y(2,j)
        y2=y(3,j)
        x(1,i)=delx*real(j-1)+(sqrt(y1*y1+4.0*y2*(ydum-y0))-y1)/(2.0*y2)
     endif
  enddo

! last point
  x(1,nsp)=delx*real(nsp-1)

! spline fit x function
! call spline with first point being quadratic
  call construct_spline0(2,nsp,dely,x)

end subroutine invert_spline1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

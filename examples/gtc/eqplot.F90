!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eqplot
  use global_parameters
  use equilibrium
  implicit none
  
  integer i,j,neq,ieq,nplot
  real dum,spdum,psi(lsp),datap(lsp),theta(lst),datat(lst),datax(lsp,lst),&
       dataz(lsp,lst),data(lsp,lst)
  real,external :: splcos,splsin

  ieq=120 !open equilibiurm plot output file
  open(ieq,file='equilibrium.out',status='replace')
101 format(i6)
102 format(e16.8)

!# of 1D radial plots
  nplot=19
  write(ieq,101)nplot,lsp

!radial axis using poloidal flux function
  do i=1,lsp
     psi(i)=spdpsi*real(i-1)
  enddo
  write(ieq,102)psi 
  
!1: sqaure-root of normalized toroidal flux function
  datap=stpp(:)
  write(ieq,102)datap 

!2: minor radius
  datap=mipp(:)
  write(ieq,102)datap 

!3: major radius
  datap=mapp(:)
  write(ieq,102)datap 

!4: Te
  datap=tepp(1,:)
  write(ieq,102)datap 

!5: ne
  datap=nepp(1,:)
  write(ieq,102)datap 

!6: ti
  datap=tipp(1,:)
  write(ieq,102)datap 

!7: zeff
  datap=zepp(1,:)
  write(ieq,102)datap 

!8: toroidal rotation
  datap=ropp(1,:)
  write(ieq,102)datap 

!9: radial electric field
  datap=erpp(1,:)
  write(ieq,102)datap 

!10: q profile
  datap=qpsi(1,:)
  write(ieq,102)datap 

!11: gcurrent profile
  datap=gpsi(1,:)
  write(ieq,102)datap 

!12: pressure profile
  datap=ppsi(1,:)
  write(ieq,102)datap 

!13: minor radius
  datap=rpsi(1,:)
  write(ieq,102)datap 

!14: toroidal flux
  datap=torpsi(1,:)
  write(ieq,102)datap

!15: error of spline cos in [0, pi/2]
  spdum=2.0*pi/real(lsp-1)
  do i=1,lsp
     dum=0.25*spdum*real(i-1)
     datap(i)=splcos(dum,spdum)-cos(dum)
  enddo
  write(ieq,102)datap
!  write(gtcout,*)'error in cos=',datap

!16: error of spline sin in [0, pi/2]
  spdum=2.0*pi/real(lsp-1)
  do i=1,lsp
     dum=0.25*spdum*real(i-1)
     datap(i)=splsin(dum,spdum)-sin(dum)
  enddo
  write(ieq,102)datap
!  write(gtcout,*)'error in sin=',datap

!17: radial grid: rgpsi
  datap=rgpsi(1,:)
  write(ieq,102)datap

!18: inverse of spline torpsi: psitor
  datap=psitor(1,:)
  write(ieq,102)datap

!19: inverse of spline rgpsi: psirg
  datap=psirg(1,:)
  write(ieq,102)datap

!# of 2D contour plots on poloidal plane
  nplot=5
  write(ieq,101)nplot,lsp,lst

! mesh points on (X,Z)
  datax=x(1,:,:)
  dataz=z(1,:,:)
  write(ieq,102)datax,dataz 

!1: b-field
  data=b(1,:,:)
  write(ieq,102)data

!2: Jacobian
  data=g(1,:,:)
  write(ieq,102)data

!3: icurrent
  data=rd(1,:,:)
  write(ieq,102)data

!4: zeta2psi
  data=nu(1,:,:)
  write(ieq,102)data

!5: delb
  data=dl(1,:,:)
  write(ieq,102)data

close(ieq)
end subroutine eqplot

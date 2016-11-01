subroutine pushi
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer i,j,k,m,ii,ierror,igyro,ij,irsmooth,icount,iir(mi)
  real(wp) cost0,sint0,cost,sint,rdum,tdum,dtime,q,r,b,g,gp,ri,rip,dbdp,dbdt,dedb,deni,&
       deltaf,upara,energy,kappa,dptdp,dptdt,dptdz,vdr,pdot,tdot,zdot,rdot,wdot,wdrift,&
       vdrenergy,delr,pi2,dmark(0:mpsi),xdum,ydum,xdot,ydot,tem,tem_inv,rfac,&
       wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,temp(0:mpsi),dtemp(0:mpsi),epara,wpara,&
       perturb,psimax,psimin,paxis,pdum,wdrive,ainv,rinv,qinv,psitmp,thetatmp,zetatmp,&
       e1,e2,e3,cmratio,cinv,wpgc(3,mi),vthi,dtem(0:mpsi),dden(0:mpsi),ddum(0:mpsi),&
       vdrtmp(0:mpsi),diagtmp(mpdiag),data1dtmp(0:mpsi,mpdata1d),&
       b1,b2,b3,b4,wbgc(4,mi),dapdp,dapdt,dapdz,vap,gyrodum,ev(4)

  delr=1.0/deltar
  pi2=2.0*pi
  psimax=0.5*a1*a1
  psimin=0.5*a0*a0
  paxis=0.5*temion*(8.0*rho0)**2
  cmratio=qion/aion
  cinv=1.0/qion
  vthi=sqrt(temion)*rho0*abs(qion)/aion
  tem_inv=1.0/(aion*vthi*vthi)
  perturb=real(nonlinear)
  vdrtmp=0.0

  if(irk==1)then
! 1st step of Runge-Kutta method
     dtime=0.5*tstep
!$omp parallel do private(m)
     do m=1,mi
        zion0(1:5,m)=zion(1:5,m)
     enddo

! 2nd step of Runge-Kutta method
  else
     dtime=tstep
     if(nonlinear==1)vdrtmp=pfluxi
  endif

  gyrodum=1.0/real(ngyroi)
! gather e_field using ngyroi-point gyro-averaging, sorting in poloidal angle
! The following line is an OpenMP directive for loop-level parallelism
! on shared memory machines (see http://www.openmp.org).

  if(magnetic==0)then
! electromagnetic fluctuations
!$omp parallel do private(m,e1,e2,e3,ev,wz1,wz0,igyro,ij,wp0,wt00,wt10,wp1,wt01,wt11)
!dir$ nextscalar
     do m=1,mi
        e1=0.0
        e2=0.0
        e3=0.0
        ev=0.0
        
        wz1=wzion(m)
        wz0=1.0-wz1

!dir$ nextscalar
        do igyro=1,ngyroi
           
           ij=jtion0(igyro,m)
           wp0=1.0-wpion(igyro,m)
           wt00=1.0-wtion0(igyro,m)
!           e1=e1+wp0*wt00*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
!           e2=e2+wp0*wt00*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
!           e3=e3+wp0*wt00*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           ev=ev+wp0*wt00*(wz0*packphi(:,0,ij)+wz1*packphi(:,1,ij))
           
           ij=ij+1
           wt10=1.0-wt00
!           e1=e1+wp0*wt10*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
!           e2=e2+wp0*wt10*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
!           e3=e3+wp0*wt10*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           ev=ev+wp0*wt10*(wz0*packphi(:,0,ij)+wz1*packphi(:,1,ij))
           
           ij=jtion1(igyro,m)
           wp1=1.0-wp0
           wt01=1.0-wtion1(igyro,m)
!           e1=e1+wp1*wt01*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
!           e2=e2+wp1*wt01*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
!           e3=e3+wp1*wt01*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           ev=ev+wp1*wt01*(wz0*packphi(:,0,ij)+wz1*packphi(:,1,ij))
           
           ij=ij+1
           wt11=1.0-wt01
!           e1=e1+wp1*wt11*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
!           e2=e2+wp1*wt11*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
!           e3=e3+wp1*wt11*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           ev=ev+wp1*wt11*(wz0*packphi(:,0,ij)+wz1*packphi(:,1,ij))
        enddo

        wpgc(1,m)=gyrodum*ev(1)
        wpgc(2,m)=gyrodum*ev(2)
        wpgc(3,m)=gyrodum*ev(3) 

        wbgc(1,m)=0.0
        wbgc(2,m)=0.0
        wbgc(3,m)=0.0
        wbgc(4,m)=0.0
     enddo

  else
! electromagnetic fields
!$omp parallel do private(m,e1,e2,e3,b1,b2,b3,b4,wz1,wz0,igyro,ij,wp0,wt00,wt10,wp1,wt01,wt11)
     do m=1,mi
        e1=0.0
        e2=0.0
        e3=0.0
        
        b1=0.0
        b2=0.0
        b3=0.0
        b4=0.0
        
        wz1=wzion(m)
        wz0=1.0-wz1

        do igyro=1,ngyroi
           
           ij=jtion0(igyro,m)
           wp0=1.0-wpion(igyro,m)
           wt00=1.0-wtion0(igyro,m)
           e1=e1+wp0*wt00*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wp0*wt00*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wp0*wt00*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           
           b1=b1+wp0*wt00*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
           b2=b2+wp0*wt00*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
           b3=b3+wp0*wt00*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
           b4=b4+wp0*wt00*(wz0*gradind(3,0,ij)+wz1*gradind(3,1,ij))
           
           ij=ij+1
           wt10=1.0-wt00
           e1=e1+wp0*wt10*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wp0*wt10*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wp0*wt10*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           
           b1=b1+wp0*wt10*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
           b2=b2+wp0*wt10*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
           b3=b3+wp0*wt10*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
           b4=b4+wp0*wt10*(wz0*gradind(3,0,ij)+wz1*gradind(3,1,ij))

           ij=jtion1(igyro,m)
           wp1=1.0-wp0
           wt01=1.0-wtion1(igyro,m)
           e1=e1+wp1*wt01*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wp1*wt01*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wp1*wt01*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

           b1=b1+wp1*wt01*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
           b2=b2+wp1*wt01*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
           b3=b3+wp1*wt01*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
           b4=b4+wp1*wt01*(wz0*gradind(3,0,ij)+wz1*gradind(3,1,ij))
           
           ij=ij+1
           wt11=1.0-wt01
           e1=e1+wp1*wt11*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wp1*wt11*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wp1*wt11*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

           b1=b1+wp1*wt11*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
           b2=b2+wp1*wt11*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
           b3=b3+wp1*wt11*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
           b4=b4+wp1*wt11*(wz0*gradind(3,0,ij)+wz1*gradind(3,1,ij))
        enddo

        wpgc(1,m)=gyrodum*e1
        wpgc(2,m)=gyrodum*e2
        wpgc(3,m)=gyrodum*e3 
        
        wbgc(1,m)=gyrodum*b1
        wbgc(2,m)=gyrodum*b2
        wbgc(3,m)=gyrodum*b3     
        wbgc(4,m)=gyrodum*b4     
     enddo
  endif

! primary ion marker temperature and parallel flow velocity
  temp=1.0
  dtemp=0.0
  temp=1.0/(temp*rtemi*aion*vthi*vthi) !inverse local temperature
  ainv=1.0/aminor

! update GC position
!$omp parallel do private(m,r,rinv,ii,wp0,wp1,tem,q,qinv,cost,sint,cost0,&
!$omp& sint0,b,g,gp,ri,rip,dbdp,dbdt,dedb,deni,upara,energy,rfac,kappa,dptdp,&
!$omp& dptdt,dptdz,epara,dapdp,vap,vdr,wdrive,wpara,wdrift,wdot,pdot,tdot,zdot,rdot,&
!$omp& dapdp,dapdt,dapdz,vap)
  do m=1,mi
     r=sqrt(2.0_wp*zion(1,m))
     rinv=1.0_wp/r
     ii=max(0,min(mpsi-1,int((r-a0)*delr)))
     wp0=real(ii+1)-(r-a0)*delr
     wp1=1.0-wp0
     tem=wp0*temp(ii)+wp1*temp(ii+1)
     q=q0+q1*r*ainv+q2*r*r*ainv*ainv
     qinv=1.0/q
     cost=cos(zion(2,m))
     sint=sin(zion(2,m))
!     cost0=cos(zion(2,m)+r*sint)
!     sint0=sin(zion(2,m)+r*sint)
     b=1.0/(1.0+r*cost)
     g=1.0
     gp=0.0
!     ri=r*r*qinv
!     rip=(2.0*q0+q1*r*ainv)*qinv*qinv
     ri=0.0
     rip=0.0
     dbdp=-b*b*cost*rinv
     dbdt=b*b*r*sint
     dedb=cinv*(zion(4,m)*zion(4,m)*qion*b*cmratio+zion(6,m)*zion(6,m))
     deni=1.0/(g*q + ri + zion(4,m)*(g*rip-ri*gp))
     upara=zion(4,m)*b*cmratio
     energy=0.5*aion*upara*upara+zion(6,m)*zion(6,m)*b
     rfac=rw*(r-rc)
     rfac=rfac*rfac
     rfac=rfac*rfac*rfac
     rfac=exp(-rfac)
     kappa=rfac
!     kappa=((min(umax*umax,energy*tem)-1.5)*kappati+kappan)*kappa*rinv
     kappa=((energy*tem-1.5)*kappati+kappan)*kappa*rinv

! perturbed electric field
     dptdp=wpgc(1,m)
     dptdt=wpgc(2,m)
     dptdz=wpgc(3,m)-wpgc(2,m)*qinv

! perturbed magnetic field
     dapdp=wbgc(1,m)
     dapdt=wbgc(2,m)
     dapdz=wbgc(3,m)-wbgc(2,m)*qinv

! parallel electric field
     epara=-(wpgc(3,m)+wbgc(4,m))*b*q*deni

! ExB drift in radial direction for w-dot and flux diagnostics
     vdr=q*(ri*dptdz-g*dptdt)*deni
     vap=-q*(ri*dapdz-g*dapdt)*deni*upara

     wdrive=kappa*(vdr+vap)
     wpara=epara*(upara-dtemp(ii))*qion*tem
     wdrift=q*(g*dbdt*dptdp-g*dbdp*dptdt+ri*dbdp*dptdz)*deni*dedb*qion*tem
     wdot=(zion0(6,m)-paranl*zion(5,m))*(wdrive+wpara+wdrift)
     
! self-consistent and external electric field for marker orbits
     dptdp=dptdp*perturb
     dptdt=dptdt*perturb
     dptdz=dptdz*perturb

! particle velocity
     pdot = q*(-g*dedb*dbdt - g*dptdt + ri*dptdz)*deni-vdrtmp(ii)
     tdot = (upara*b*(1.0-q*gp*zion(4,m)) + q*g*(dedb*dbdp + dptdp))*deni
     zdot = (upara*b*q*(1.0+rip*zion(4,m)) - q*ri*(dedb*dbdp + dptdp))*deni
     rdot = ((gp*zion(4,m)-1.0)*(dedb*dbdt + paranl*dptdt)-&
          paranl*q*(1.0+rip*zion(4,m))*dptdz)*deni 
         
! update particle position
!     if(zion0(1,m) < paxis)then
! particles close to axis use (x,y) coordinates
!        pdum=sqrt(zion(1,m))
!        xdum   = pdum*cost  
!        ydum   = pdum*sint
!        pdum=1.0/zion(1,m)
!        xdot   = 0.5*pdot*xdum*pdum-ydum*tdot
!        ydot   = 0.5*pdot*ydum*pdum+xdum*tdot
!        pdum=sqrt(zion0(1,m))
!        xdum   = pdum*cos(zion0(2,m)) + dtime*xdot
!        ydum   = pdum*sin(zion0(2,m)) + dtime*ydot
!        zion(1,m) = max(1.0e-8_wp*psimax,xdum*xdum+ydum*ydum)
!        zion(2,m) = sign(1.0_wp,ydum)*acos(max(-1.0_wp,min(1.0_wp,xdum/sqrt(zion(1,m)))))
!     else
     zion(1,m) = max(1.0e-8_wp*psimax,zion0(1,m)+dtime*pdot)
     zion(2,m) = zion0(2,m)+dtime*tdot
!     endif

     zion(3,m) = zion0(3,m)+dtime*zdot
     zion(4,m) = zion0(4,m)+dtime*rdot
     zion(5,m) = zion0(5,m)+dtime*wdot
     
     zion(2,m)=modulo(zion(2,m),pi2)
     zion(3,m)=modulo(zion(3,m),pi2)

! store GC information for flux measurements
     wpgc(1,m)=vdr
     wpgc(2,m)=energy
     wpgc(3,m)=b
  enddo

  if(irk==2)then

! out of boundary particle
!$omp parallel do private(m)
     do m=1,mi
        if(zion(1,m) > psimax)then
           zion(1,m)=zion0(1,m)
           zion(2,m)=2.0*pi-zion0(2,m)
           zion(3,m)=zion0(3,m)
           zion(4,m)=zion0(4,m)
           zion(5,m)=zion0(5,m)
          
        elseif(zion(1,m) < psimin)then
           zion(1,m)=zion0(1,m)
           zion(2,m)=2.0*pi-zion0(2,m)
           zion(3,m)=zion0(3,m)
           zion(4,m)=zion0(4,m)
           zion(5,m)=zion0(5,m)
        endif
     enddo

! Restore temperature profile when running a nonlinear calculation
     if(nonlinear==1)then

!$omp parallel do private(m,r)
        do m=1,mi
           r=sqrt(2.0*zion(1,m))
           iir(m)=max(0,min(mpsi,int((r-a0)*delr)))
        enddo
       
        dtem=0.0
        dden=0.0
        dmark=0.0
        do m=1,mi
           ii=iir(m)
           dtem(ii)=dtem(ii)+wpgc(2,m)*zion(5,m)
           dmark(ii)=dmark(ii)+wpgc(1,m)
           dden(ii)=dden(ii)+1.0
        enddo
        
! S.Ethier 06/04/03  According to the MPI standard, the send and receive
! buffers cannot be the same for MPI_Reduce or MPI_Allreduce. 
        icount=mpsi+1
        ddum=0.0
        call MPI_ALLREDUCE(dtem,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
        dtem=ddum
        call MPI_ALLREDUCE(dmark,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
        dmark=ddum
        call MPI_ALLREDUCE(dden,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
        dden=ddum

        irsmooth=int(sqrt(real(mpsi)))
        dmark=dmark/max(1.0_wp,dden) !radial marker flux
        do i=1,irsmooth
           rdum=dmark(1)
           tdum=dmark(mpsi-1)
           dmark(1:mpsi-1)=0.5*dmark(1:mpsi-1)+0.25*(dmark(0:mpsi-2)+dmark(2:mpsi))
           dmark(0)=0.5*(dmark(0)+rdum)
           dmark(mpsi)=0.5*(dmark(mpsi)+tdum)
        enddo
        tdum=0.1
        pfluxi=(1.0-tdum)*pfluxi+tdum*dmark
        
! remove small scale temperature perturbation
        irsmooth=mpsi
        dtem=dtem*tem_inv/max(1.0_wp,dden)
        do i=1,irsmooth
           rdum=dtem(1)
           tdum=dtem(mpsi-1)
           dtem(1:mpsi-1)=0.5*dtem(1:mpsi-1)+0.25*(dtem(0:mpsi-2)+dtem(2:mpsi))
           dtem(0)=0.5*(dtem(0)+rdum)
           dtem(mpsi)=0.5*(dtem(mpsi)+tdum)
        enddo
        tdum=0.01
        rdtemi=(1.0-tdum)*rdtemi+tdum*dtem
       
!$omp parallel do private(m,ii)
        do m=1,mi
           ii=iir(m)
           zion(5,m)=zion(5,m)-(wpgc(2,m)*tem_inv-1.5)*rdtemi(ii)
        enddo
     endif
  endif
  
  if(idiag==0)then
! fluxes diagnose at irk=1
     diagion=0.0
     dden=0.0
     data1di=0.0
     do m=1,mi
        r=sqrt(2.0*zion0(1,m))
        rinv=1.0_wp/r
        ii=max(0,min(mpsi,int((r-a0)*delr+0.5)))
        dden(ii)=dden(ii)+1.0

        deltaf=zion0(5,m)
        upara=wpgc(3,m)*zion0(4,m)/rho0
        energy=wpgc(2,m)/(rho0*rho0)-1.5
        vdr=rinv*wpgc(1,m)

! radial profile of particle and energy flux
        data1di(ii,1)=data1di(ii,1)+vdr*deltaf
        data1di(ii,2)=data1di(ii,2)+vdr*deltaf*energy

!!! ion diagnosis: density,entropy,flow,energy,fluxes of particle,momentum,heat
        diagion(1)=diagion(1)+deltaf
        diagion(2)=diagion(2)+deltaf*deltaf
        diagion(3)=diagion(3)+upara
        diagion(4)=diagion(4)+upara*deltaf
        diagion(5)=diagion(5)+energy
        diagion(6)=diagion(6)+energy*deltaf
        diagion(7)=diagion(7)+vdr*deltaf
        diagion(8)=diagion(8)+vdr*upara*deltaf
        diagion(9)=diagion(9)+vdr*energy*deltaf
     enddo     
     diagion(10)=real(mi)

! sum over all MPI processes
     icount=mpdiag
     diagtmp=0.0
     call MPI_REDUCE(diagion,diagtmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     diagion=diagtmp
        
     icount=(mpsi+1)*mpdata1d
     data1dtmp=0.0
     call MPI_ALLREDUCE(data1di,data1dtmp,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)

     icount=mpsi+1
     ddum = 0. !---wj
     call MPI_ALLREDUCE(dden,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)

! radial profile data normailized by marker #
    if(mype==0) then !---wj
     do i=1,mpdata1d
        data1di(:,i)=data1dtmp(:,i)/ddum
     enddo
     endif !---wj
  endif
 
end subroutine pushi

subroutine pushe(icycle,irke)
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer i,im,ip,j,jt,m,ii,j11,j10,j01,j00,ierror,ij,icount,&
      ipjt,icycle,irke,irsmooth,jtgc0(memax),jtgc1(memax)
  real(wp) cost0,sint0,cost,sint,rdum,tdum,dtime,q,r,b,g,gp,ri,rip,dbdp,dbdt,&
      dedb,deni,deltaf,upara,energy,kappa,dptdp,dptdt,dptdz,vdr,pdot,tdot,zdot,rdot,&
      wdot,wdrift,vdrenergy,delr,delt(0:mpsi),delz,pi2_inv,pi2,xdum,ydum,xdot,ydot,&
      g_inv,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,temp(0:mpsi),dtemp(0:mpsi),tem,rfac,&
      qdum,wpara,psimax,psimin,paxis,pdum,wdrive,ainv,rinv,qinv,perturb,psitmp,&
      thetatmp,zetatmp,rhoi,e1,e2,e3,e4,cmratio,cinv,pptpt,dzonal,tem_inv,dtem(0:mpsi),&
      dden(0:mpsi),ddum(0:mpsi),vdrtmp(0:mpsi),dmark(0:mpsi),diagtmp(mpdiag),&
      wpgc(memax),wzgc(memax),wtgc0(memax),wtgc1(memax),data1dtmp(mpsi+1,mpdata1d),ev(4)
! need to use memax for array size declaration here since me changes in this subroutine
 
  delr=1.0/deltar
  delt=2.0*pi/deltat
  delz=1.0/deltaz
  pi2=2.0*pi
  pi2_inv=0.5/pi
  g_inv=1.0/rho0
  psimax=0.5*a1*a1
  psimin=0.5*a0*a0
  paxis=0.5*(8.0*rho0)**2
  cmratio=qelectron/aelectron
  cinv=1.0/qelectron
  tem_inv=1.0/(rho0*rho0)
  perturb=real(nonlinear)
  vdrtmp=0.0

  if(irke==1)then
! 1st step of Runge-Kutta method
     dtime=0.25*tstep/real(ncycle)

     if(icycle==1)then
! save electron info
        if(irk*ihybrid==1)then
           me1=me
!$omp parallel do private(m)
           do m=1,me
              zelectron1(1:nparam,m)=zelectron(1:nparam,m)
           enddo
        else

! electron start from initial position
           me=me1
!$omp parallel do private(m)
           do m=1,me
              zelectron(1:nparam,m)=zelectron1(1:nparam,m)
           enddo
        endif
     endif

!$omp parallel do private(m)
     do m=1,me
        zelectron0(1:5,m)=zelectron(1:5,m)
     enddo

     if(track_particles == 1)then
!$omp parallel do private(m)
        do m=1,me
           zelectron0(7:nparam,m)=zelectron(7:nparam,m)
        enddo
     endif

! 2nd step of Runge-Kutta method
  else
     dtime=0.5*tstep/real(ncycle)
     if(nonlinear==1 .and. ihybrid==nhybrid .and. irk==2)vdrtmp=pfluxe

  endif

!$omp parallel do private(m,psitmp,thetatmp,zetatmp,r,ip,jt,ipjt,wz1,&
!$omp& rdum,ii,wp1,im,tdum,j00,j01)
  do m=1,me
     psitmp=zelectron(1,m)
     thetatmp=zelectron(2,m)
     zetatmp=zelectron(3,m)
     
     r=sqrt(2.0*psitmp)
     ip=max(0,min(mpsi,int((r-a0)*delr+0.5)))
     jt=max(0,min(mtheta(ip),int(thetatmp*pi2_inv*delt(ip)+0.5)))
     ipjt=igrid(ip)+jt
     wz1=(zetatmp-zeta0)*delz
     wzgc(m)=wz1
     
     rdum=delr*max(0.0_wp,min(a1-a0,r-a0))
     ii=max(0,min(mpsi-1,int(rdum)))
     wp1=rdum-real(ii)
     wpgc(m)=wp1

! inner flux surface
     im=ii
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+10.0
     tdum=(tdum-aint(tdum))*delt(im)
     j00=max(0,min(mtheta(im)-1,int(tdum)))
     jtgc0(m)=igrid(im)+j00 !jtgc0 used here for poloidal grid label
     wtgc0(m)=tdum-real(j00)
     
! outer flux surface
     im=ii+1
     tdum=pi2_inv*(thetatmp-zetatmp*qtinv(im))+10.0
     tdum=(tdum-aint(tdum))*delt(im)
     j01=max(0,min(mtheta(im)-1,int(tdum)))
     jtgc1(m)=igrid(im)+j01
     wtgc1(m)=tdum-real(j01)
  enddo

! gather e_field for drift kinetic electron
!$omp parallel do private(m,e1,e2,e3,e4,ev,wz1,wz0,ij,wp0,wt00,wt10,wp1,wt01,wt11)
!dir$ nextscalar
  do m=1,me
!      e1=0.0
!      e2=0.0
!      e3=0.0
!      e4=0.0
     ev=0.0
     wz1=wzgc(m)
     wz0=1.0-wz1

     ij=jtgc0(m)
     wp0=1.0-wpgc(m)
     wt00=1.0-wtgc0(m)
!      e1=e1+wp0*wt00*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
!      e2=e2+wp0*wt00*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
!      e3=e3+wp0*wt00*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
!      e4=e4+wp0*wt00*(wz0*phit(0,ij)+wz1*phit(1,ij))
     ev=ev+wp0*wt00*(wz0*packphi(:,0,ij)+wz1*packphi(:,1,ij))
     
     ij=ij+1
     wt10=1.0-wt00
!      e1=e1+wp0*wt10*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
!      e2=e2+wp0*wt10*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
!      e3=e3+wp0*wt10*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
!      e4=e4+wp0*wt10*(wz0*phit(0,ij)+wz1*phit(1,ij))
     ev=ev+wp0*wt10*(wz0*packphi(:,0,ij)+wz1*packphi(:,1,ij))
        
     ij=jtgc1(m)
     wp1=1.0-wp0
     wt01=1.0-wtgc1(m)
!      e1=e1+wp1*wt01*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
!      e2=e2+wp1*wt01*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
!      e3=e3+wp1*wt01*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
!      e4=e4+wp1*wt01*(wz0*phit(0,ij)+wz1*phit(1,ij))
     ev=ev+wp1*wt01*(wz0*packphi(:,0,ij)+wz1*packphi(:,1,ij))
     
     ij=ij+1
     wt11=1.0-wt01
!      e1=e1+wp1*wt11*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
!      e2=e2+wp1*wt11*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
!      e3=e3+wp1*wt11*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
!      e4=e4+wp1*wt11*(wz0*phit(0,ij)+wz1*phit(1,ij))
     ev=ev+wp1*wt11*(wz0*packphi(:,0,ij)+wz1*packphi(:,1,ij))
     
     wpgc(m)=ev(1)
     wtgc0(m)=ev(2)
     wzgc(m)=ev(3)
     wtgc1(m)=ev(4)
  enddo

  temp=1.0
  dtemp=0.0
  temp=1.0/(temp*rteme)
  ainv=1.0/aminor
! update GC position
!$omp parallel do private(m,r,rinv,ii,wp0,wp1,tem,dzonal,q,qinv,&
!$omp& cost,sint,cost0,sint0,&
!$omp& b,g,gp,ri,rip,dbdp,dbdt,dedb,deni,upara,energy,rfac,kappa,dptdp,&
!$omp& dptdt,dptdz,pptpt,wpara,wdrift,vdr,wdrive,wdot,pdot,tdot,zdot,rdot)
  do m=1,me
     r=sqrt(2.0*zelectron(1,m))
     rinv=1.0/r
     ii=max(0,min(mpsi-1,int((r-a0)*delr)))
     wp0=real(ii+1)-(r-a0)*delr
     wp1=1.0-wp0
     tem=g_inv*g_inv*(wp0*temp(ii)+wp1*temp(ii+1))
     dzonal=wp0*phip00(ii)+wp1*phip00(ii+1)
     q=q0+q1*r*ainv+q2*r*r*ainv*ainv
     qinv=1.0/q
     cost=cos(zelectron(2,m))
     sint=sin(zelectron(2,m))
!     cost0=cos(zelectron(2,m)+r*sint)
!     sint0=sin(zelectron(2,m)+r*sint)
     b=1.0/(1.0+r*cost)
     g=1.0
     gp=0.0
!     ri=r*r*qinv
!     rip=(2.0*q0+q1*r*ainv)*qinv*qinv
     ri=0.0
     rip=0.0
     dbdp=-b*b*cost*rinv
     dbdt=b*b*r*sint
     dedb=cinv*(zelectron(4,m)*zelectron(4,m)*qelectron*b*cmratio+zelectron(6,m)*&
          zelectron(6,m))
     deni=1.0/(g*q + ri + zelectron(4,m)*(g*rip-ri*gp))
     upara=zelectron(4,m)*b*cmratio
     energy=0.5*aelectron*upara*upara+zelectron(6,m)*zelectron(6,m)*b
     rfac=rw*(r-rc)
     rfac=rfac*rfac
     rfac=rfac*rfac*rfac
     rfac=exp(-rfac)
     kappa=rfac
!     kappa=((min(umax*umax,energy*tem)-1.5)*kappate+kappan)*kappa*rinv
     kappa=((energy*tem-1.5)*kappate+kappan)*kappa*rinv

! perturbed quantities
     dptdp=wpgc(m)
     dptdt=wtgc0(m)
     dptdz=wzgc(m)-wtgc0(m)*qinv
     pptpt=wtgc1(m)
     
! ExB drift in radial direction for w-dot and flux diagnostics
     vdr=q*(ri*dptdz-g*dptdt)*deni
     wdrive=vdr*kappa
     wpara=pptpt*qelectron*tem
     wdrift=q*(g*dedb*dbdt + g*dptdt - ri*dptdz)*deni*dzonal*qelectron*tem*rinv  !!XY
     wdot=(zelectron0(6,m)-paranl*zelectron(5,m))*(wdrive+wpara+wdrift)

! self-consistent and external electric field for marker orbits
     dptdp=dptdp*perturb
     dptdt=dptdt*perturb
     dptdz=dptdz*perturb

! particle velocity
     pdot = q*(-g*dedb*dbdt - g*dptdt + ri*dptdz)*deni-vdrtmp(ii)
     tdot = (upara*b*(1.0-q*gp*zelectron(4,m)) + q*g*(dedb*dbdp + dptdp))*deni
     zdot = (upara*b*q*(1.0+rip*zelectron(4,m)) - q*ri*(dedb*dbdp + dptdp))*deni
     rdot = ((gp*zelectron(4,m)-1.0)*(dedb*dbdt + paranl*dptdt)-&
          paranl*q*(1.0+rip*zelectron(4,m))*dptdz)*deni 
         
! update particle position
!     if(zelectron0(1,m) < paxis)then
! particles close to axis use (x,y) coordinates
!        pdum=sqrt(zelectron(1,m))
!        xdum   = pdum*cost  
!        ydum   = pdum*sint
!        pdum=1.0/zelectron(1,m)
!        xdot   = 0.5*pdot*xdum*pdum-ydum*tdot
!        ydot   = 0.5*pdot*ydum*pdum+xdum*tdot
!        pdum=sqrt(zelectron0(1,m))
!        xdum   = pdum*cos(zelectron0(2,m)) + dtime*xdot
!        ydum   = pdum*sin(zelectron0(2,m)) + dtime*ydot
!        zelectron(1,m) = max(1.0e-8*psimax,xdum*xdum+ydum*ydum)
!        zelectron(2,m) = sign(1.0,ydum)*acos(max(-1.0,min(1.0,xdum/sqrt(zelectron(1,m)))))
!     else
     zelectron(1,m) = max(1.0e-8*psimax,zelectron0(1,m)+dtime*pdot)
     zelectron(2,m) = zelectron0(2,m)+dtime*tdot
!     endif

     zelectron(3,m) = zelectron0(3,m)+dtime*zdot
     zelectron(4,m) = zelectron0(4,m)+dtime*rdot
     zelectron(5,m) = zelectron0(5,m)+dtime*wdot
     
! S.Ethier 04/13/2004  The hand-coded procedure gives wrong answers in some
! cases so it is better to use the modulo() function.
!!!     zelectron(2,m)=zelectron(2,m)*pi2_inv+10.0 !period of 1
!!!     zelectron(2,m)=2.0*pi*(zelectron(2,m)-aint(zelectron(2,m))) ![0,2*pi)
!!!     zelectron(3,m)=zelectron(3,m)*pi2_inv+10.0
!!!     zelectron(3,m)=2.0*pi*(zelectron(3,m)-aint(zelectron(3,m)))

     zelectron(2,m)=modulo(zelectron(2,m),pi2)
     zelectron(3,m)=modulo(zelectron(3,m),pi2)

! store GC information for flux measurements
     wpgc(m)=vdr
     wtgc0(m)=energy
     wzgc(m)=b
  enddo
      
 if(irke==2)then
! out of boundary particle
!$omp parallel do private(m)
    do m=1,me
       if(zelectron(1,m) > psimax)then
          zelectron(1,m)=zelectron0(1,m)
          zelectron(2,m)=2.0*pi-zelectron0(2,m)
          zelectron(3,m)=zelectron0(3,m)
          zelectron(4,m)=zelectron0(4,m)
          zelectron(5,m)=zelectron0(5,m)
          
       elseif(zelectron(1,m) < psimin)then
          zelectron(1,m)=zelectron0(1,m)
          zelectron(2,m)=2.0*pi-zelectron0(2,m)
          zelectron(3,m)=zelectron0(3,m)
          zelectron(4,m)=zelectron0(4,m)
          zelectron(5,m)=zelectron0(5,m)
       endif
    enddo

! restore temperature profile
    if(nonlinear==1 .and. ihybrid==nhybrid .and. &
         icycle==2*ncycle .and. irk==2)then

!$omp parallel do private(m,r)
       do m=1,me
          r=sqrt(2.0*zelectron(1,m))
          jtgc0(m)=max(0,min(mpsi,int((r-a0)*delr))) !jtgc0 is now used for radial grid
       enddo

       dtem=0.0
       dden=0.0
       dmark=0.0
       do m=1,me
          ii=jtgc0(m)
          dtem(ii)=dtem(ii)+wtgc0(m)*zelectron(5,m)
          dmark(ii)=dmark(ii)+wpgc(m)
          dden(ii)=dden(ii)+1.0
       enddo

       icount=mpsi+1
       ddum=0.0
       call MPI_ALLREDUCE(dtem,ddum,icount,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
       dtem=ddum
       call MPI_ALLREDUCE(dmark,ddum,icount,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
       dmark=ddum
       call MPI_ALLREDUCE(dden,ddum,icount,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
       dden=ddum

       irsmooth=int(sqrt(real(mpsi)))
       dmark=dmark/max(1.0_wp,dden) !radial particle flow
       do i=1,irsmooth
          rdum=dmark(1)
          tdum=dmark(mpsi-1)
          dmark(1:mpsi-1)=0.5*dmark(1:mpsi-1)+0.25*(dmark(0:mpsi-2)+dmark(2:mpsi))
          dmark(0)=0.5*(dmark(0)+rdum)
          dmark(mpsi)=0.5*(dmark(mpsi)+tdum)
       enddo
       tdum=0.1
       pfluxe=(1.0-tdum)*pfluxe+tdum*dmark
       
       irsmooth=mpsi
       dtem=dtem*tem_inv/max(1.0_wp,dden) !perturbed thermal energy
       do i=1,irsmooth
          rdum=dtem(1)
          tdum=dtem(mpsi-1)
          dtem(1:mpsi-1)=0.5*dtem(1:mpsi-1)+0.25*(dtem(0:mpsi-2)+dtem(2:mpsi))
          dtem(0)=0.5*(dtem(0)+rdum)
          dtem(mpsi)=0.5*(dtem(mpsi)+tdum)
       enddo
       tdum=0.01
       rdteme=(1.0-tdum)*rdteme+tdum*dtem
       
!$omp parallel do private(m,ii)
       do m=1,me
          ii=jtgc0(m)
          zelectron(5,m)=zelectron(5,m)-(wtgc0(m)*tem_inv-1.5)*rdteme(ii)
       enddo
    endif
 endif

 if(idiag==0 .and. ihybrid*icycle*irke==1)then
! electron diagnosis
    dden=0.0
    data1de=0.0
    diagelectron=0.0
    do m=1,me
       r=sqrt(2.0*zelectron(1,m))
       rinv=1.0/r
       ii=max(0,min(mpsi,int((r-a0)*delr+0.5)))
       dden(ii)=dden(ii)+1.0
       
       deltaf=zelectron0(5,m)
       upara=wzgc(m)*zelectron0(4,m)/rho0
       energy=wtgc0(m)/(rho0*rho0)-1.5
       vdr=rinv*wpgc(m)

! radial profile of particle and energy fluxes       
       data1de(ii,1)=data1de(ii,1)+vdr*deltaf
       data1de(ii,2)=data1de(ii,2)+vdr*energy*deltaf
       
!!! volume averaged: density,entropy,flow,energy,fluxes of particle,momentum,heat
       diagelectron(1)=diagelectron(1)+deltaf
       diagelectron(2)=diagelectron(2)+deltaf*deltaf
       diagelectron(3)=diagelectron(3)+upara
       diagelectron(4)=diagelectron(4)+upara*deltaf
       diagelectron(5)=diagelectron(5)+energy
       diagelectron(6)=diagelectron(6)+energy*deltaf
       diagelectron(7)=diagelectron(7)+vdr*deltaf
       diagelectron(8)=diagelectron(8)+vdr*upara*deltaf
       diagelectron(9)=diagelectron(9)+vdr*energy*deltaf
    enddo
    diagelectron(10)=real(me)

! sum over all MPI processes
    icount=mpdiag
    diagtmp=0.0
    call MPI_REDUCE(diagelectron,diagtmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    diagelectron=diagtmp

    icount=(mpsi+1)*mpdata1d
    data1dtmp=0.0
    call MPI_ALLREDUCE(data1de,data1dtmp,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)

    icount=mpsi+1    
    ddum=0. !---wj
    call MPI_ALLREDUCE(dden,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror) !!!!XY

! particle data normalized by marker #
    if(mype==0) then !---wj
    do i=1,mpdata1d
       data1de(:,i)=data1dtmp(:,i)/ddum
    enddo
    endif !---wj
 endif
 
end subroutine pushe

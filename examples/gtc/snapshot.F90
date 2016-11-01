subroutine snapshot
  use global_parameters
  use particle_array
  use field_array
  implicit none
  
  integer,parameter :: nfield=3,nvgrid=65,iosnap=222
  integer mtgrid,nf,ns,m,ip,jenergy,kpitch,nsnap,i,ierror,j,jt,icount
  real(wp) upara,delf,rmu,fullf,ai_inv,vi_inv,ae_inv,ve_inv,emax_inv,delr,r,b,&
       energy,pitch,wt,tdum,pdum,profile(0:mpsi,6,nspecies),pdf(nvgrid,4,nspecies),&
       dfield(mgrid),proftmp(0:mpsi,6,nspecies),pdftmp(nvgrid,4,nspecies)
  real(wp),dimension(:),allocatable :: eachflux
  real(wp),dimension(:,:),allocatable :: allflux
  real(wp),dimension(:,:,:),allocatable :: poloidata,fluxdata
  character(len=64) cdum

! particle species: 1=ion, 2=electron, 3=EP
! radial profiles: density, flow, energy
  profile=0.0

! distribution function: energy, pitch angle
  pdf=0.0

! species ns=1: ion
  ns=1
  ai_inv=1.0/aion
  vi_inv=1.0/rho0
  delr=1.0/deltar
  emax_inv=0.2
  do m=1,mi
     r=sqrt(2.0*zion(1,m))
     b=1.0/(1.0+r*cos(zion(2,m)))
     delf=zion(5,m)
     rmu=zion(6,m)
     fullf=zion0(6,m)
     upara=zion(4,m)*b*qion*ai_inv*vi_inv 
     energy=0.5*upara*upara+rmu*rmu*b*ai_inv*vi_inv*vi_inv
     pitch=upara/sqrt(2.0*energy)

! particles sorted into bins in radius, energy, and pitch
     ip=max(0,min(mpsi,int((r-a0)*delr+0.5)))! radial bin
     jenergy=1+max(0,min(nvgrid-1,int(real(nvgrid-1)*energy*emax_inv)))! energy bin
     kpitch=1+max(0,min(nvgrid-1,int(real(nvgrid-1)*0.5*(pitch+1.0))))! pitch bin

     profile(ip,1,ns)=profile(ip,1,ns)+fullf
     profile(ip,2,ns)=profile(ip,2,ns)+delf
     profile(ip,3,ns)=profile(ip,3,ns)+fullf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+delf*upara
     profile(ip,5,ns)=profile(ip,5,ns)+fullf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+delf*energy
     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf
  enddo

! species ns=2: electron
  if(nhybrid>0)then
  ns=2
  ae_inv=1.0/aelectron
  ve_inv=sqrt(aelectron/aion)/rho0
  do m=1,me
     r=sqrt(2.0*zelectron(1,m))
     b=1.0/(1.0+r*cos(zelectron(2,m)))
     delf=zelectron(5,m)
     rmu=zelectron(6,m)
     fullf=zelectron0(6,m)
     upara=zelectron(4,m)*b*qelectron*ae_inv*ve_inv 
     energy=0.5*upara*upara+rmu*rmu*b*ae_inv*ve_inv*ve_inv
     pitch=upara/sqrt(2.0*energy)

     ip=max(0,min(mpsi,int((r-a0)*delr+0.5)))
     jenergy=1+max(0,min(nvgrid-1,int(real(nvgrid-1)*energy*emax_inv)))
     kpitch=1+max(0,min(nvgrid-1,int(real(nvgrid-1)*0.5*(pitch+1.0))))

     profile(ip,1,ns)=profile(ip,1,ns)+fullf
     profile(ip,2,ns)=profile(ip,2,ns)+delf
     profile(ip,3,ns)=profile(ip,3,ns)+fullf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+delf*upara
     profile(ip,5,ns)=profile(ip,5,ns)+fullf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+delf*energy
     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf
  enddo
  endif

! sum over MPI tasks
  icount=nspecies*6*(mpsi+1)
  call MPI_REDUCE(profile,proftmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  profile=proftmp

  icount=nspecies*4*(nvgrid)
  call MPI_REDUCE(pdf,pdftmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  pdf=pdftmp
  
! normalization by marker #
  if(mype==0)then !---wj
  do ns=1,nspecies
     do ip=0,mpsi
        profile(ip,2:6,ns)=profile(ip,2:6,ns)/profile(ip,1,ns)
     enddo
     do j=1,nvgrid
        pdf(j,2,ns)=pdf(j,2,ns)/max(1.0,pdf(j,1,ns))
        pdf(j,4,ns)=pdf(j,4,ns)/max(1.0,pdf(j,3,ns))
     enddo
  enddo
  endif !---wj

! poloidal resolution=poloidal grid on iflux surface
  mtgrid=mtheta(iflux)

  allocate(poloidata(0:mtgrid,0:mpsi,nfield+2),fluxdata(0:mtgrid,mtoroidal,nfield),&
       eachflux(mtgrid),allflux(mtgrid,mtoroidal),SOURCE=transfer(int(B'01111111100000100000000000000000'),real(1))) !---wj
  poloidata=0.0
  fluxdata=0.0

! field quantities: phi, a_para, fluidne. Last two coloumn of poloidal for coordinates
  do nf=1,nfield
     if(nf==1)dfield=phi(0,:)
     if(nf==2)dfield=apara(0,:)
     if(nf==3)dfield=fluidne(0,:)

! gather potential on flux surface
     allflux=0.0
     do j=1,mtgrid
        eachflux(j)=dfield(igrid(iflux)+j)
     enddo

     icount=mtgrid
     call MPI_GATHER(eachflux,icount,mpi_Rsize,allflux,icount,mpi_Rsize,0,toroidal_comm,ierror) 
     fluxdata(1:mtgrid,:,nf)=allflux

! poloidal BC
     fluxdata(0,:,nf)=fluxdata(mtgrid,:,nf)

! potential data on poloidal plain uses polar coordinates
     do j=0,mtgrid
        tdum=2.0*pi*real(j)/real(mtgrid)
        do i=0,mpsi
           jt=max(0,min(mtheta(i),1+int(tdum/deltat(i))))
           wt=tdum/deltat(i)-real(jt-1)
           poloidata(j,i,nf)=(wt*dfield(igrid(i)+jt)+(1.0-wt)*dfield(igrid(i)+jt-1))/rho0**2
        enddo
     enddo
  enddo

! poloidal grid position in polar coordinates
  do j=0,mtgrid
     tdum=2.0*pi*real(j)/real(mtgrid)
     do i=0,mpsi
        r=a0+deltar*real(i)
        poloidata(j,i,nfield+1)=1.0+r*cos(tdum)
        poloidata(j,i,nfield+2)=r*sin(tdum)
     enddo
  enddo

! open snapshot output file
  if(mype==0)then
     nsnap=mstepall+istep
     write(cdum,'(i5.5,".out")')nsnap
     cdum='snap'//trim(cdum)
     open(iosnap,file=cdum,status='replace')

! parameters: # of species, fields, and grids in velocity, radius, poloidal, toroidal; T_up
     write(iosnap,101)nspecies,nfield,nvgrid,mpsi+1,mtgrid+1,mtoroidal
     write(iosnap,102)1.0/emax_inv

! write out particle radial profile and pdf, and 2D field
     write(iosnap,102)profile,pdf,poloidata,fluxdata

! close snapshot file
     close(iosnap)
  endif

101 format(i6)
102 format(e10.4)
  
end subroutine snapshot

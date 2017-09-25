#key parameters:
OPENMP=y
DOUBLE_PRECISION=n
ADIOS=n
RESTART_TIME=y
PETSc=n
64BITS=y
XT3=n
DEBUG=n

#machine name:
dragon=n

############################################################################
#             Makefile to build the GTC code
#           ==================================
#
# You only need to type "gmake" to build the code on the platforms
# defined below. The makefile runs the "uname -s" command to detect
# the operating system automatically. Other options are:
#
#  % gmake OPENMP=y       Builds the code with OpenMP support
#  % gmake OPENMP=n       Builds the code WITHOUT OpenMP support
#
#  % gmake DOUBLE_PRECISION=y  Builds with 8-byte floating point precision
#
#  % gmake DEBUG=y        Compiles the files with debug option (-g)
#                         The default is no debug option
#
#  % gmake 64BITS=y       To compile a 64-bit version on AIX
#  % gmake 64BITS=n       To compile a 32-bit version on AIX
#                         The default is 32 bits on Linux clusters
#
#  % gmake PGI=y          Use the PGI compiler (pgf90) on Linux. The default
#                         is to use the Lahey-Fujitsu compiler lf95.
#
#  % gmake ALTIX=y        Compiles with Intel compilers on the Altix
#                         using ifort ... -lmpi
#
# Special targets:
#
#  % gmake clean      Removes the executable and all object files (*.o)
#
#  % gmake cleanomp   Removes the executable and the object files
#                     containing OpenMP directives
#
#  % gmake doc        Rebuilds the documentation.
#
#############################################################################
# Default executable name
CMD:=gtc
ALL:${CMD}
LIB:=

# Common file for fft routine using the Glassman FFT source code.
FFT:=fft_gl.o

ifeq (${PETSc},y)
   PETSC_solver = petsc.o
   PETSC_OPT= -D__PETSc
   include ${PETSC_DIR}/bmake/common/base
endif

# In the next declaration, "uname" is a standard Unix command which prints the name of the current system.
os:= $(shell uname -s)

# inquire the hostname to judge the system name
#jaguar	-> Jaguar
#nid -> Franklin 
#b -> Bassi
SYSTEMS = $(shell hostname | sed -e 's/[[:digit:]]\+//')


# We initialize to an empty string the variable that contains the compiler option for the OpenMP directives.
OMPOPT:=

ifeq ($(os),IRIX64)
    # Flags for the SGI MIPSpro compilers (hecate)
    # If we use "gmake OPENMP=y", the -mp flag is added to the compilation
    # options to take into account the OpenMP directives in the code.
    # If we specify "gmake SPEEDSHOP=y ...", the code is linked with the
    # SPEEDSHOP library libss.so to recognize the "calipers" inserted in the
    # code to do precise profiling of certain loops. Those calipers are
    # used by the SPEEDSHOPpro profiling tool on IRIX (see: man ssrun).
    CMP:=f90
    F90C:=f90
    OPT:=-64 -mips4 -Ofast -freeform -I/usr/pppl/include
    LIB:=-lmpi
    ifeq ($(OPENMP),y)
      OMPOPT:=-mp
      CMD:=gtc
    endif
    ifeq ($(DOUBLE_PRECISION),y)
        OPT:=$(OPT) -DDOUBLE_PRECISION
    endif

    ifeq ($(SPEEDSHOP),y)
      OPT:=$(OPT) -D__SPEEDSHOP
      LIB:=$(LIB) -lss
    endif
    ifeq ($(DEBUG),y)
      OPT:=-g $(OPT)
    endif
endif

ifeq ($(os),AIX)
    # Flags for the IBM AIX compiler (Machine IBM SP seaborg)
    # stacksize:=0.5GB, datasize:=1.75GB
    # We also define the "__NERSC" symbol to be passed to the preprocessor
    # through the compiler option "-WF,-D__NERSC"
    # -g -pg
    CMP:=mpxlf90_r
    OPT:= -qsuffix=cpp=F90 -WF,-D__NERSC \
        -qsuffix=f=f90 -qfree=f90 -qinitauto \
        -qarch=auto -qtune=auto -qcache=auto -qstrict -O3 -Q -u
    # -O3 is recommended by IBM;  without -qstrict the compiler 
    # might bend some IEEE rules and one has to check carefully the result;
    # -qhot  improves nested loops 
    ifdef 64BITS
#       OPT:=-qdpc -qautodbl=dbl4 -q64 -WF,-DMPI_REAL=MPI_DOUBLE $(OPT)
      OPT:=-q64 $(OPT)
#       HOME_ETHIER:=/usr/common/homes/e/ethier
      LIB:=~/timers/system_stats.o $(NAG) $(MASS) $(NETCDF)
    else
      OPT:=-bmaxstack:0x20000000 -bmaxdata:0x7000000 $(OPT)
      LIB:=/usr/common/usg/timers/system_stats.o $(NAG) $(MASS) $(NETCDF)
    endif
    # If we use "gmake OPENMP=y" then we compile with mpxlf90_r and add
    # the -qsmp=omp option to take into account the OpenMP directives.
    ifdef OPENMP
      CMP:=mpxlf90_r
      OMPOPT:=-qsmp=omp -qnosave
    endif
    ifdef DEBUG
      OPT:=-g $(OPT)
    endif
    ifeq ($(DOUBLE_PRECISION),y)
        OPT:=$(OPT) -WF,-DDOUBLE_PRECISION
    endif
endif

# Settings for Linux platform. The default is to use the LF95 compiler
ifeq ($(os),Linux)
  # Flags for the Linux system
  # Default compiler and options: Lahey-Fujitsu compiler
    CMP:=mpifort
    F90C:=pgf95

    ##OPT:=-O --ap --tpp --ntrace --staticlink -I/usr/local/lff95/include
    OPT:=-O `pkg-config netcdf-fortran --cflags` 
	LIB:=`pkg-config netcdf-fortran --libs`
  ifeq ($(PGI),y)
    MPIMODULE:=/usr/pppl/pgi/5.2-1/mpich-1.2.6/include/f90base
    F90C:=pgf90
    OPT:=-O -D__PGF90 -Mfree -Kieee
    LIB:=
  endif
  ifeq ($(XT3),y)
    PORTALS=n
    CMP:=mpif90
    F90C:=pgf90
    OPT:=-fastsse 
    ##-I$(NETCDFDIR)/include ##-DADIOS=1
    #OPT:=-fastsse ${HDF5_FLIB}
    #LIB:= -L$(NETCDFDIR)/lib -lnetcdf
  endif

  ifeq ($(SYSTEMS),jaguar)
    CMP:=ftn
    F90C:=ftn
  endif
  
  ifeq ($(SYSTEMS),jaguarpf-login)
    CMP:=ftn 
    F90C:=ftn
#    OPT:=-Mbounds
  endif

  ifeq ($(SYSTEMS),nid)
    CMP:=ftn
    F90C:=ftn
  endif

  ifeq ($(SYSTEMS),ewok)
    CMP:=mpif90
    F90C:=mpif90
  endif

  ifeq ($(SYSTEMS),pike)
    OPT:=-ffree-line-length-none
  endif

  ifeq ($(INTEL),y)
    F90C:=ifort
    OPT:=-O
  endif
  ifeq ($(PATHSCALE),y)
    F90C:=pathf90
    #OPT:=-O3 -static
    OPT:=-O3
  endif
  ifeq ($(ALTIX),y)
    CMP:=ifort
    ###OPT:=-O3 -ipo -ftz -stack_temps
    OPT:=-O3 -ftz-
    #OPT:=-O3 -ipo -ftz -stack_temps
    ###OPT:=-O -g -ftz
    LIB:=-lmpi
  endif

  ifeq ($(DOUBLE_PRECISION),y)
      OPT:=-DDOUBLE_PRECISION $(OPT)
  endif
  ifeq ($(ADIOS),y)
    OPT:=-DADIOS=1 $(OPT)
  endif
  ifeq ($(RESTART_TIME),y)
    OPT:=-DRESTART_TIME=1 $(OPT)
  endif


  ifeq ($(OPENMP),y)
      OMPOPT:=
  endif

  ifeq ($(DEBUG),y)
    OPT:=-g -C $(OPT)
  endif
endif

# Settings for the NEC SX-6 machine and Earth Simulator. For these machines,
# a cross-compiler on a different platform is usually used. In that case,
# the name returned by the $(os) command will not be SUPER-UX. The cross-
# compiler for the SX machine is usually called "sxmpif90", and "esmpif90"
# for the Earth Simulator.
ifeq ($(os),SUPER-UX)
    CMP:=sxmpif90
    OPT:=-f4 -C vopt -Wf'-pvctl loopcnt=8000000 vwork=stack -L fmtlist mrgmsg transform source' -ftrace -R2 -D_SX
    ###OPT:=-f4 -C vsafe -Wf,-pvctl loopcnt=10000000 -R2 -ftrace -D_SX
  # The next line is for David Skinner's MPI profiling library
    ##LIB:=-L/S/n003/home003/m0052/ipm/newer/ipm -lipm -mpiprof
    ###LIB:=-Wl,-Z8G,-m
    ifeq ($(DOUBLE_PRECISION),y)
        OPT:=-DDOUBLE_PRECISION $(OPT)
    endif
    ifeq ($(FFTSX6),y)
      FFT:=fft_sx6.o
      LIB:=$(LIB) -lfft
    endif
    ifeq ($(DEBUG),y)
      OPT:=-Cdebug $(OPT)
    endif
    SETUP:=setup_vec.o
    CHARGEI:=chargei_vec.o
    PUSHI:=pushi_vec.o
    POISSON:=poisson_vec.o
    SHIFTI:=shifti_vec.o
endif

# Settings for the CRAY-X1
ifeq ($(os),UNICOS/mp)
    CMP:=ftn
    OPT:= -D_CRAYX1 -Ostream2 -Otask0 -rm
    ##OPT:= -Ostream2 -Otask0 -rm -rd
    ifeq ($(DEBUG),y)
      OPT:=-g $(OPT)
    endif
    ifeq ($(DOUBLE_PRECISION),y)
        OPT:=-DDOUBLE_PRECISION $(OPT)
    endif
    # Set options for FFT
    ifeq ($(FFTCRAY),y)
      FFT:=fft_cray.o
    endif
    SETUP:=setup_vec.o
    CHARGEI:=chargei_vec.o
    PUSHI:=pushi_vec.o
    POISSON:=poisson_vec.o
    SHIFTI:=shifti_vec.o
endif

# to use ADIOS, we need this defined
#    LIBD=/autofs/na1_home/ge1/work/jaguar/lib
#    INCD=/autofs/na1_home/ge1/work/jaguar/include
     #LIBD=-L/ccs/proj/e2e/adios/jaguar/lib 
     #INCD=-I/ccs/proj/e2e/adios/jaguar/include
#LIBD=
#INCD=

#ifeq ($(SYSTEMS),ewok)
    LIBD=/ccs/home/hardes/AdioDist/ewok/lib
    INCD=/ccs/home/hardes/AdioDist/ewok/include
    LIBD=${ADIOS_LIB}
    INCD=${ADIOS_INC}
#endif

ifeq ($(SYSTEMS),jaguar)
    LIBD=${ADIOS_LIB} 
    INCD=${ADIOS_INC}	
#    LIBD=/ccs/home/sklasky/e2e/adios-svn/jaguar/lib
#    INCD=/ccs/home/sklasky/e2e/adios-svn/jaguar/include
endif
ifeq ($(SYSTEMS),jaguarpf-login)
    LIBD=${ADIOS_LIB} 
    INCD=${ADIOS_INC}	
     #LIBD=-L/ccs/proj/e2e/adios/jaguar/lib 
     #INCD=-I/ccs/proj/e2e/adios/jaguar/include
    #LIBD=${ADIOS_LIB} 
    #INCD=${ADIOS_INC}	
    #LIBD=/ccs/home/sklasky/adios/trunk/jaguarpf/lib
    #INCD=/ccs/home/sklasky/adios/trunk/jaguarpf/include
endif

ifeq ($(SYSTEMS),nid)
    LIBD=/u0/x/xiao/AdioFranklin/lib
    INCD=/u0/x/xiao/AdioFranklin/include
endif

ifeq ($(ADIOS),y)
    #LIBADIOS=$(LIBD) -ladios -ladios_phdf5 -lmxml ${HDF5_CLIB} 
    LIBD=-L/ccs/proj/e2e/sklasky/adios-svn/jaguar.pgi/lib
    INCD=-I/ccs/proj/e2e/sklasky/adios-svn/jaguar.pgi/include
    LIBADIOS=$(LIBD) -ladios -lmxml ${HDF5_CLIB} 
    INCADIOS=$(INCD)
    NETCDF_LIB=
else
    LIBADIOS=
    INCADIOS=
    ifeq ($(dragon),y)
      NETCDF_LIB=-I/usr/pgi/linux86-64/6.1/include  -I/usr/lib/netcdf-3.6.1/include \
      -L/usr/lib/netcdf-3.6.1/lib -lnetcdf
    endif
    ##ifeq ($(SYSTEMS),ewok)
    ##  NETCDF_LIB=-I/apps/netcdf/3.6.1/ewok_x86_64_pgi625/lib \
    ##  -I/apps/netcdf/3.6.1/ewok_x86_64_pgi625/include -L/apps/netcdf/3.6.1/ewok_x86_64_pgi625/lib -lnetcdf
    ##endif
     NETCDF_LIB=${NETCDF_FLIB}#
    ifeq ($(SYSTEMS),jaguar)
      NETCDF_LIB=${NETCDF_CLIB}
    endif
    ifeq ($(SYSTEMS),jaguarpf-login)
      NETCDF_LIB=${NETCDF_CLIB}
    endif
    ifeq ($(SYSTEMS),nid)
      NETCDF_LIB=${NETCDF}
    endif
    ifeq ($(SYSTEMS),pike)
      NETCDF_LIB=-I/usr/include -L/usr/lib/ -lnetcdff
    endif
endif

ifeq ($(dragon),y)
    NETCDF_LIB=-I/usr/pgi/linux86-64/6.1/include  -I/usr/lib/netcdf-3.6.1/include \
    -L/usr/lib/netcdf-3.6.1/lib -lnetcdf
endif

ifeq ($(PORTALS),y)
    PORTALS_LIBS=-lportals
    PORTALS_OBJ=restart_ptl.o
    USE_PORTALS=1
else
    PORTALS_LIBS=
    PORTALS_OBJ=
    USE_PORTALS=0
endif

OPT:=$(INCADIOS) $(PETSC_OPT) $(OPT)
LIB:=$(LIB) $(LIBADIOS) $(NETCDF_LIB) $(PORTALS_LIBS)

##################################################################
# We add ".F90" to the list of suffixes to allow source files on which the
# co-processor will be run automatically.
.SUFFIXES: .o .F90 .F

# List of all the object files needed to build the code
OBJ:=module.o main.o function.o setup.o load.o restart.o diagnosis.o snapshot.o chargei.o\
     poisson.o smooth.o field.o pushi.o shifti.o chargee.o pushe.o shifte.o tracking.o\
     dataout3d.o pushfield.o eqdata.o eqplot.o packup_phi.o sort.o $(PETSC_solver) $(FFT)

$(CMD): $(OBJ)
	$(CMP) $(OMPOPT) $(OPT) -o $(CMD) $(OBJ) $(LIB) ${PETSC_FORTRAN_LIB} ${PETSC_KSP_LIB}

module.o : module.F90
	$(CMP) $(OMPOPT) $(OPT) -c  module.F90

$(OBJ): module.o

.F90.o : module.o
	$(CMP) $(OMPOPT) $(OPT) $(LIB) -c  $<
.c.o : 
	mpicc -O2 $(OMPOPT) $(LIB) -c  $<


# The following tag is meant to "clean" the directory by removing the
# executable along with all the object files created by the compilation 
# One only has to run:  gmake clean

clean:
	rm -f $(CMD) $(OBJ) *.mod


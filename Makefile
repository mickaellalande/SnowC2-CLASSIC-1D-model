# Makefile for CLASSIC.

# Usage:

# make mode=<mode> debug=<debug> r4=<r4>

#   or, for the default mode, just:

# make

# where <mode> is:
#
# serial           - compiles code using the serial GNU compiler on the local platform (default)
# parallel         - compiles code using the parallel GNU compiler on the local platform
#
# serial_intel     - compiles code using the serial Intel compiler on ECCC systems
# parallel_intel   - compiles code using the parallel Intel compiler on ECCC systems
#
# serial_gnu       - compiles code using the serial GNU compiler on ECCC systems 

# where <debug> is:
# on               - enables compiler switches for debugging, otherwise "off" if not defined (default)
#                    *** NB: Intel compiler may catch more uninitialized variables than gfortran since it 
#                            supports the initialization of allocatable arrays to NaNs.

# where <r4> is:
# on               - compiles the code with default real type set to 4-bytes long, 
#                    otherwise "off" if not defined (i.e. 8-bytes is the default)

# For ECCC systems, please use the "make_classic.sh" script instead of using "make" directly with this Makefile.
# The script sets the specific environment for each of the supported modes.

# ==================================================================================================================

# Modifications: 
#   Ed Chan, Jun 2022.
#     - Added support for serial GNU compiler on ECCC U2 platforms.
#   Ed Chan, Apr 2022.
#     - Modified for ECCC U2 systems.
#   Ed Chan, Feb 2021.
#     - Added option to compile reals in 32-bit precision instead of 64-bit.
#     - Stopped initializing variables to zero by the compiler.
#   Ed Chan, Nov 2020.
#     - Added debugging options.
#     - Added serial mode for ppp Intel compiler. 
#   Ed Chan, Mar 2020.
#     - Added separate serial modes for ppp and singularity container.
#     - Updated supercomputer mode to use skylake specific compiler option.
#   Ed Chan, Feb 2020.
#     - Updated serial mode to use netCDF4/HDF5 libraries as for the other modes.
#   Ed Chan, Dec 2019.
#     - Revised supercomputer/ppp compiler/linker options, which are now the same on both platforms.
#     - On ppp, support for the parallel netCDF library must now be enabled PRIOR to running make.
#   Ed Chan, Sep 2018.
#     - Updated various compiler options.
#     - Added "cray" and "ppp" modes.
#     - All objects go into directories specific to the mode used.
#     - All executable names are labelled with the mode appended.

# Object files
OBJ = classicParams.o ctemStateVars.o classStateVars.o generalUtils.o peatlandsMod.o \
	fileIOModule.o ctemUtilities.o calcLandSurfParams.o groundAlbedo.o mvidx.o snowProcesses.o \
	snowAlbedoTransmissMod.o soilHeatFluxPrep.o waterInfiltrateUnsat.o \
	canopyInterception.o waterFlowNonInfiltrate.o energyBudgetPrep.o waterInfiltrateSat.o canopyAlbedoTransmiss.o \
	classGatherScatter.o waterFlowInfiltrate.o snowSublimation.o energBalVegSolve.o waterCalcPrep.o methaneProcesses.o \
	ctemGatherScatter.o applyAllometry.o photosynCanopyConduct.o canopyWaterUpdate.o \
	atmosphericVarsCalc.o canopyPhaseChange.o iceSheetBalance.o waterUnderCanopy.o energBalNoVegSolve.o errorHandler.o classGrowthIndex.o \
	DIASURFZ2.o screenRelativeHumidity.o pondedWaterFreeze.o snowTempUpdate.o checkWaterBudget.o energyBudgetDriver.o DRCOEF2.o SLDIAG2.o waterUpdates.o snowHeatCond.o \
	radiationDriver.o waterBudgetDriver.o FLXSURFZ2.o snowInfiltrateRipen.o snowMelt.o soilWaterPhaseChg.o soilProperties.o energyWaterBalanceCheck.o snowAddNew.o \
	soilHeatFluxCleanup.o waterBaseflow.o  balanceCarbon.o autotrophicRespiration.o phenology.o  \
	turnoverMod.o mortality.o disturbance.o competitionMod.o landuseChangeMod.o \
	allocateCarbon.o heterotrophicRespirationMod.o soilCProcesses.o \
	tracer.o readFromJobOptions.o outputManager.o modelStateDrivers.o dynamicTiling.o tiledDisturbance.o ctemDriver.o \
	prepareOutputs.o metModule.o mainCore.o main.o xmlParser.o \
	xmlManager.o n_processes.o bnf.o nvolatil.o nitrification.o denitrification.o nleach.o \
	nuptake.o nallocate.o nlitter.o nhumific.o nMineralImm.o balnitro.o CLASSIC.o

# Object files (.o and .mod) are placed into directories named according to the mode. In the serial
# compiler case, these are placed into the "objectFiles" directory.
ODIR = objectFiles_$(mode)

# Compile with 32 bit precision, if switch is enabled.
ifeq ($(r4), on)
	label=_32bit
	ODIR=objectFiles_$(mode)_32bit
else
	ifeq ($(mode), parallel)
                R8 = -fdefault-real-8
        else ifeq ($(mode), serial)
                R8 = -fdefault-real-8
        else ifeq ($(mode), serial_gnu)
		R8 = -fdefault-real-8
	else
		R8 = -r8
	endif
endif

# Set variables based on the "mode" specified.

# *** NB: The first 3 options are specifically for ECCC platforms. 
#         Other users will need to find their own recipes for the last 2 options (i.e."parallel" and "serial").

ifeq ($(mode), parallel_intel)
	# Parallel Intel compiler for ECCC systems.
	# NB: Support for the Intel compiler and parallel netCDF4/HDF5 libraries must first be enabled prior to running make via:
	#     . ssmuse-sh -x /fs/ssm/main/opt/intelcomp/inteloneapi-2022.1.2/intelcomp+mpi+mkl
	#     . ssmuse-sh -x /fs/ssm/main/opt/hdf5-netcdf4/parallel/intelmpi-2022.1.2/static/inteloneapi-2022.1.2/01
	# MPI wrapper to the default parallel Fortran compiler, which will be Intel after the SSM packages are loaded.
	COMPILER = mpif90
        # Fortran Flags.
	FFLAGS = -DPARALLEL -Dwithout_agcm_ $(R8) -O2 -mp1 -xCORE-AVX512 -align array64byte -static-intel -static-libgcc -traceback -module $(ODIR)
	ifeq ($(debug), on)
		FFLAGS = -DPARALLEL -Dwithout_agcm_ $(R8) -g -O0 -init=snan,arrays -static-intel -static-libgcc -traceback -fpe0 -module $(ODIR)
	endif
        # Library Flags.
	LFLAGS = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl
else ifeq ($(mode), serial_intel)
	# Serial Intel compiler for ECCC systems.
	# NB: Support for the Intel compiler and serial netCDF4/HDF5 libraries must first be enabled prior to running make via:
	#     . ssmuse-sh -x /fs/ssm/main/opt/intelcomp/inteloneapi-2022.1.2/intelcomp+mpi+mkl
	#     . ssmuse-sh -x /fs/ssm/main/opt/hdf5-netcdf4/serial/static/inteloneapi-2022.1.2/01
	COMPILER = ifort
	# Fortran Flags.
	FFLAGS = -Dwithout_agcm_ $(R8) -g -O2 -mp1 -xCORE-AVX512 -align array64byte -static-intel -static-libgcc -traceback -module $(ODIR)
	ifeq ($(debug), on)
		FFLAGS = -Dwithout_agcm_ $(R8) -g -O0 -init=snan,arrays -static-intel -static-libgcc -traceback -fpe0 -module $(ODIR)
	endif
	# Library Flags.
	LFLAGS = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl
else ifeq ($(mode), serial_gnu)
	# Serial gfortran compiler for ECCC systems.
	COMPILER = gfortran
	ODIR = objectFiles$(label)
	# Fortran Flags.
	FFLAGS = -Dwithout_agcm_ -O3 -g $(R8) -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR)
	ifeq ($(debug), on)
		FFLAGS = -Dwithout_agcm_ -O0 -g $(R8) -ffree-line-length-none -finit-real=snan -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR) #-Wall -Wextra
	endif
	# Include Flags.
	IFLAGS = -I/usr/lib64/gfortran/modules
	# Library Flags.
	LFLAGS = -lnetcdff
else ifeq ($(mode), parallel)
	# MPI wrapper to the default parallel Fortran compiler.
	COMPILER = mpif90
        # Fortran Flags. The following is specific to the gfortran compiler.
	FFLAGS = -DPARALLEL -Dwithout_agcm_ -O3 -g $(R8) -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR)
	ifeq ($(debug), on)
		FFLAGS = -DPARALLEL -Dwithout_agcm_ -O0 -g $(R8) -ffree-line-length-none -finit-real=snan -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR)
	endif
        # Include Flags.
	IFLAGS = -I$(shell nf-config --includedir)
        # Library Flags.
	LFLAGS = -L$(shell nf-config --flibs)
else
	# Serial gfortran compiler.
	COMPILER = gfortran
	mode = serial
	ODIR = objectFiles$(label)
	# Fortran Flags.
	FFLAGS = -Dwithout_agcm_ -O3 -g $(R8) -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR)
	ifeq ($(debug), on)
		FFLAGS = -Dwithout_agcm_ -O0 -g $(R8) -ffree-line-length-none -finit-real=snan -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR) #-Wall -Wextra
	endif
	# Include Flags.
	IFLAGS = -I$(shell nf-config --includedir)
	# Library Flags.
	LFLAGS = -L$(shell nf-config --flibs)
endif

# Create required directory/.gitignore file, if missing.
VOID := $(shell mkdir -p $(ODIR))

# RECIPES
# Compile object files from .F90 sources
$(ODIR)/%.o: src/%.F90
	$(COMPILER) $(FFLAGS) $(IFLAGS) -c $< -o $@

# Compile object files from .f90 sources
$(ODIR)/%.o: src/%.f90
	$(COMPILER) $(FFLAGS) $(IFLAGS) -c $< -o $@

# Compile object files from .f (Fortran 77) sources
$(ODIR)/%.o: src/%.f
	$(COMPILER) $(FFLAGS) $(IFLAGS) -c $< -o $@

# Properly reference the ODIR for the linking
OBJD = $(patsubst %,$(ODIR)/%,$(OBJ))

# Link objects together and put executable in the bin/ directory
CLASSIC: $(OBJD)
	$(COMPILER) $(FFLAGS) $(IFLAGS) -o bin/CLASSIC_$(mode)$(label) $(OBJD) $(LFLAGS)

# "make clean mode=supercomputer" removes all object files in objectFiles_supercomputer
clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.mod bin/CLASSIC_$(mode)

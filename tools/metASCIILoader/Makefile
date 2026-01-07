.PHONY: all clean
# Serial compiler
#COMPILER = gfortran
# Serial Include Flags
#IFLAGS = -O3 -g -I/usr/include
# Serial Library Flags
#LFLAGS = -lnetcdff
# Additional flags
#LFLAGS += -ldl -lz -lm -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow  -fbounds-check -Wall -Wextra

#IFLAGS = -O3 -g -I/fs/ssm/hpco/exp/hdf5-netcdf4/serial/static/gnu-4.8.5/01/ubuntu-18.04-amd64-64/include
#LFLAGS = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -ldl -lz -lcurl
#LFLAGS += -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow  -fbounds-check -Wall -Wextra

COMPILER = ifort

# Fortran Flags.
#FFLAGS = -r8 -g -O2 -mp1 -xCORE-AVX512 -align array64byte -static-intel -static-libgcc -traceback -module $(ODIR)
#ifeq ($(debug), on)
FFLAGS = -r8 -g -O0 -init=snan,arrays -static-intel -static-libgcc -traceback -fpe0 #-module $(ODIR)
#endif
# Library Flags.
LFLAGS = -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl

# RECIPES
# Compile object files from .F90 sources
%.o: src/%.F90
	$(COMPILER) $(FFLAGS) -c $< -o $@ $(LFLAGS)

# Compile object files from .f90 sources
%.o: src/%.f90
	$(COMPILER) $(FFLAGS) -c $< -o $@ $(LFLAGS)

all: fileIOModule.o metASCIILoader.o
	$(COMPILER) $(FFLAGS) -o bin/metASCIILoader fileIOModule.o metASCIILoader.o $(LFLAGS)

# "make clean" removes all object files and executables
clean:
	rm -f *.o *.mod *~ core bin/metASCIILoader
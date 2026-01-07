#!/bin/bash
  
# Script to compile CLASSIC on ECCC Ubuntu 22.04 systems.
# NB: tested only on gpscc3 (the default mode on collab).

# Usage:

#   make_classic_collab.sh <mode> [debug=on] [r4=on]
#   make_classic_collab.sh mode=<mode> [debug=on] [r4=on]

# If <mode> is not specified, the default is to compile for the local host.

# Otherwise, if <mode> is:
#   parallel_intel    - uses parallel netCDF libraries with the Intel compiler
#   serial_intel      - uses serial netCDF libraries with the Intel compiler
#   serial            - uses serial netCDF libraries with the GNU compiler

# The debug/r4 switches are optional and are off by default.

# Modifications:

# Initial Implementation: Ed Chan, Feb. 2024 based on script for ECCC U2 systems.

# --------------------------------------------------------------------------------------------

# Set the default mode to the local host and the debug/r4 switches to off.

debug=off
r4=off

case "$ORDENV_TRUEHOST" in
  gpscc3)           localmode=parallel_intel ;;
  *)                localmode=serial_gnu ;;
esac

mode=$localmode

# Override the defaults if set on the command line.

for arg in $@
do
  case $arg in
    *=*) eval $arg ;;
      *) mode=$arg ;;
  esac
done

# Run required commands prior to make (if any), depending on the platform/compiler.

case $mode in
  parallel_intel)
    # Load Intel compiler and parallel netCDF4/HDF5 libraries.
    . ssmuse-sh -x main/opt/intelcomp/inteloneapi-2023.2.0/intelcomp+mpi+mkl >/dev/null 2>&1
    . ssmuse-sh -x main/opt/hdf5-netcdf4/parallel/alllib/inteloneapi-2023.2.0/01 2>&1
  ;;
  serial_intel)
    # Load Intel compiler and serial netCDF4/HDF5 libraries.
    . ssmuse-sh -x main/opt/intelcomp/inteloneapi-2023.2.0/intelcomp+mpi+mkl >/dev/null 2>&1
    . ssmuse-sh -x main/opt/hdf5-netcdf4/serial/alllib/inteloneapi-2023.2.0/01 >/dev/null 2>&1 
  ;;
  serial)
    # Load GNU compiler and serial netCDF4/HDF5 libraries (assumes these are already installed on the system).
    mode=serial
  ;;
  *)
    echo "Option $mode is not supported. Please specify: parallel_intel or serial_intel."
    echo "If no parameter is supplied, the local default mode is $localmode."
    exit
  ;;
esac

# Remove unsupported compiler flag from Makefile, adjust libraries, and save result to temporary file.
sed 's/-xCORE-AVX512//; s/-lhdf5_hl *-lhdf5 *-lz *-lcurl/-lm/' Makefile > makefile

/usr/bin/make -f makefile mode=$mode debug=$debug r4=$r4
rm makefile

exit

# Use your own paths here for testing this script within a built container.
#cd /home/jskye99/ECCC2022ContractWork/container_ops/CLASSIC_singularity_container/para_netcdf_hdf-4.6.3/setup
#cd ..
rm -rf tarballs
mkdir tarballs

### Intel Fortran Compiler stuff:
# rm -rf compiler
# mkdir compiler
###

curl -L --output tarballs/hdf5-1.10.9.tar.bz2 https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.9/src/hdf5-1.10.9.tar.bz2
curl -L --output tarballs/netcdf-c-4.9.0.tar.gz https://github.com/Unidata/netcdf-c/archive/v4.9.0.tar.gz
curl -L --output tarballs/netcdf-fortran-4.5.4.tar.gz https://downloads.unidata.ucar.edu/netcdf-fortran/4.5.4/netcdf-fortran-4.5.4.tar.gz
curl -L --output tarballs/zlib-1.2.12.tar.gz https://zlib.net/fossils/zlib-1.2.12.tar.gz

### Intel Fortran Compiler stuff:
#echo "setup_netcdf STEP 2.4"
#curl -L --output compiler/l_fortran-compiler_p_2022.1.0.134.sh https://registrationcenter-download.intel.com/akdlm/irc_nas/18703/l_fortran-compiler_p_2022.1.0.134.sh
###

rm -rf src
mkdir src
cd src

tar -xvjf ../tarballs/hdf5-1.10.9.tar.bz2
tar -xvzf ../tarballs/netcdf-c-4.9.0.tar.gz
tar -xvzf ../tarballs/netcdf-fortran-4.5.4.tar.gz
tar -xvzf ../tarballs/zlib-1.2.12.tar.gz

pwd

cp --verbose ../hdf5/dobuildsrc hdf5-1.10.9
cp --verbose ../netc/dobuildsrc netcdf-c-4.9.0
cp --verbose ../netf/dobuildsrc netcdf-fortran-4.5.4
cp --verbose ../zlib/dobuildsrc zlib-1.2.12

cd ..
cp do* ..

### Intel Fortran Compiler stuff:
# echo "setup_netcdf STEP 7"
# cd compiler
# chmod +x l_fortran-compiler_p_2022.1.0.134.sh
# mkdir ifort_package
###

cd ..
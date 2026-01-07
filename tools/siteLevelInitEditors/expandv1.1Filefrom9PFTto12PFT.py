import xarray as xr
import os, fnmatch
import shutil
import numpy as np
import sys


# M. Brady. Jan 2021

# This script assumes you are trying to convert CLASSIC v 1.01.x initialization
# files to be able to run with two shrub and one sedge PFTs.

#Coordinate variables that are to be changed:
#Increase icp1 coordinate by 1
# - It is the number of CLASS (physics) PFTs plus bareground so in CLASS 1.0
#   this is 5. CLASS 1.1 has 5 PFTs plus bare so becomes 6.
# - insert the new position at the second-last index
#Increase ic coordinate by 1
# - Same as above but it doesn't include the bareground.
# - insert the new position at the last index
#Increase icc coordinate by 3
# - This is the number of CTEM (biogeochem) PFTs . In CLASSIC 1.0
#   this is 9 but we add two shrubs and sedges in CLASSIC 1.1 so this becomes 12.
#Increase iccp2 coordinate by 3
# - Same as above but this one includes 2 extra spots. One for bareground
#   and the other for land use change C.

# NOTE: it expects the init files to have the form ??-???_init.nc, e.g. AU-Tum_init.nc
# All the init files should be using the same naming convention.

# Set the upper directory for the original CLASSIC v1.1 (9 PFT) init files.
dir_path = str(sys.argv[1])

def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


for filename in find_files(dir_path, '??-???_init.nc'):

    print (filename)

    ds = xr.open_dataset(filename)

    # maintain same variable order as source file
    variable_order = list(ds.data_vars)
    new_das = {k:None for k in variable_order}

    icp1_vars = [var for var in ds.data_vars if 'icp1' in ds[var].coords]
    #print(icp1_vars)
    # create the new icp1 array
    new_icp1 = xr.DataArray(np.arange(1, 7, dtype='int32'), dims='icp1', attrs=ds['icp1'].attrs)
    for var in icp1_vars:
        data = ds[var]
        # preallocate the correct-size array for variable data
        new_data = np.zeros((1,6,1,1))
        # copy data from original array for the first 4 PFTs
        new_data[:,:4,:,:] = data.data[:,:4,:,:]
        # copy the data for bareground into the last index
        new_data[:,5,:,:] = data.data[:,4,:,:]
        # instantiate a new DataArray using as much from the old one as possible
        new_da = xr.DataArray(
            data=new_data,
            dims=['tile', 'icp1','lat','lon'],
            coords={
                'tile': data['tile'],
                'icp1': new_icp1,
                'lat': data['lat'],
                'lon': data['lon']
            },
            attrs=data.attrs
        )
        # copy encoding from source
        new_da.encoding = data.encoding
        # append to our dict of "new" dataarrays
        new_das[var] = new_da.to_dataset(name=var)

    # same as icp1 but goes from 4 to 5
    ic_vars = [var for var in ds.data_vars if 'ic' in ds[var].coords]
    #print(ic_vars)
    # create the new ic array
    new_ic = xr.DataArray(np.arange(1, 6, dtype='int32'), dims='ic', attrs=ds['ic'].attrs)
    for var in ic_vars:
        data = ds[var]
        # preallocate the correct-size array for variable data
        new_data = np.zeros((1,5,1,1))
        # copy data from original array, leaving a zero in the last index
        new_data[:,:4,:,:] = data.data
        # instantiate a new DataArray using as much from the old one as possible
        new_da = xr.DataArray(
            data=new_data,
            dims=['tile', 'ic','lat','lon'],
            coords={
                'tile': data['tile'],
                'ic': new_ic,
                'lat': data['lat'],
                'lon': data['lon']
            },
            attrs=data.attrs
        )
        # copy encoding from source
        new_da.encoding = data.encoding
        # append to our dict of "new" dataarrays
        new_das[var] = new_da.to_dataset(name=var)

    icc_vars = [var for var in ds.data_vars if 'icc' in ds[var].coords]
    #print(icc_vars)
    # create the new icc array
    new_icc = xr.DataArray(np.arange(1, 13, dtype='int32'), dims='icc', attrs=ds['icc'].attrs)
    for var in icc_vars:
        data = ds[var]
        # preallocate the correct-size array for variable data
        new_data = np.zeros((1,12,1,1))
        # copy data from original array, leaving a zero in the last 3 indexes
        new_data[:,:9,:,:] = data.data
        # instantiate a new DataArray using as much from the old one as possible
        new_da = xr.DataArray(
            data=new_data,
            dims=['tile', 'icc','lat','lon'],
            coords={
                'tile': data['tile'],
                'icc': new_icc,
                'lat': data['lat'],
                'lon': data['lon']
            },
            attrs=data.attrs
        )
        # copy encoding from source
        new_da.encoding = data.encoding
        # append to our dict of "new" dataarrays
        new_das[var] = new_da.to_dataset(name=var)
        
    iccp2_vars = [var for var in ds.data_vars if 'iccp2' in ds[var].coords]
    #print(iccp2_vars)
    # create the new iccp2 array
    new_iccp2 = xr.DataArray(np.arange(1,15,dtype='int32'), dims='iccp2', attrs=ds['iccp2'].attrs)
    for var in iccp2_vars:
        data = ds[var]
        # preallocate the correct-size array for variable data
        new_data = np.zeros((1,20,14,1,1))
        # copy data from original array, first 9 indexes
        new_data[:,:,:9,:,:] = data.data[:,:,:9,:,:]
        # copy data, 10th and 11th
        new_data[:,:,12:14:,:,:] = data.data[:,:,9:,:,:]
        # instantiate a new DataArray using as much from the old one as possible
        new_da = xr.DataArray(
            data=new_data,
            dims=['tile', 'layer', 'iccp2','lat','lon'],
            coords={
                'tile': data['tile'],
                'layer': data['layer'],
                'iccp2': new_iccp2,
                'lat': data['lat'],
                'lon': data['lon']
            },
            attrs=data.attrs,
        )
        # copy encoding from source
        new_da.encoding = data.encoding
        # append to our dict of "new" dataarrays
        new_das[var] = new_da.to_dataset(name=var)

    # now just copy the variables that don't need adjustment
    other_vars = [var for var in ds.data_vars if var not in icp1_vars + icc_vars + ic_vars + iccp2_vars]
    #print(other_vars)
    for var in other_vars:
        new_das[var] = ds[var]

    # create the new Dataset with variables ordered the same as source Dataset
    new_ds = xr.merge(new_das.values())
    new_ds.attrs = ds.attrs
    #new_ds

    # Rename the old file  to show that it is the 9 PFT version.
    os.rename(filename, filename.replace('.nc', '_orig_nine_pfts.nc'))

    # Write this to netcdf. We keep the same name so it is not necessary to update all the
    # joboptions files.
    new_ds.to_netcdf(filename)

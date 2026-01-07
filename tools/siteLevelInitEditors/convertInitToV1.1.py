##!/space/hall3/sitestore/eccc/crd/ccrp/mib001/usr_conda/envs/py3-base/bin/python
import xarray as xr
import os, fnmatch
import shutil
import numpy as np
import sys

# J. Melton. Jan 2021

# This script assumes you are trying to convert CLASSIC v 1.0.x initialization
# files to be able to run with CLASSIC v.1.1.x. The major changes then are 
# to add the non-structural and structural C pools, the N-related variables, make the litter and soil C
# per layer instead of bulk, and to add peatsoilC (which is ignored in non-peatland runs
# anyway). NOTE: these files will need to be respun up by the model since the pools
# are set to somewhat arbitrary new values. The model switch spinfast should then
# be set to >1 for that spinup to avoid triggering the C balance check.

# Also NOTE: it expects the init files to have the form ??-???_init.nc, e.g. AU-Tum_init.nc
# All the init files should be using the same naming convention.

# Set the upper directory for the original CLASSIC v1.0 init files.
dir_path = str(sys.argv[1])

def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


for filename in find_files(dir_path, '??-???_init.nc'):
    
    print (filename)
    
    #Open the original init file
    orig = xr.open_dataset(filename)
    
    # Rename it to show that it is the v1.0 version. 
    os.rename(filename, filename.replace('.nc', '_v1.0.x.nc'))

    # rename coordinate
    orig = orig.rename({'icctem' : 'icc'})

    # read in some coordinates:
    lat = orig['lat']
    lon = orig['lon']
    tile = orig['tile']
    layer = orig['layer']
    icc = orig['icc']
    iccp2 = orig['iccp2']

    # read in the gleaf, stem, and root mass to find an estimate of the _ns that we can start with.
    gleafmas = orig['gleafmas']
    gleafmas_ns = 0.1 * gleafmas
    gleafmas_s = 0.9* gleafmas

    stemmass = orig['stemmass']
    stemmass_ns = 0.1 * stemmass
    stemmass_s = 0.9* stemmass

    rootmass = orig['rootmass']
    rootmass_ns = 0.1 * rootmass
    rootmass_s = 0.9* rootmass

    # now add the new variables to the data structure:
    arrgleafmas_ns = xr.DataArray(gleafmas_ns, coords=[tile, icc, lat, lon ], dims=('tile','icc','lat', 'lon'))
    orig['gleafmas_ns'] = arrgleafmas_ns
    orig['gleafmas_ns'].attrs['long_name']='Non-structural green leaf mass'
    orig['gleafmas_ns'].attrs['_FillValue']=-999
    orig['gleafmas_ns'].attrs['units']='kg C/m^2'

    arrgleafmas_s = xr.DataArray(gleafmas_s, coords=[tile, icc, lat, lon ], dims=('tile','icc','lat', 'lon'))
    orig['gleafmas_s'] = arrgleafmas_s
    orig['gleafmas_s'].attrs['long_name']='Structural green leaf mass'
    orig['gleafmas_s'].attrs['_FillValue']=-999
    orig['gleafmas_s'].attrs['units']='kg C/m^2'

    arrstemmass_ns = xr.DataArray(stemmass_ns, coords=[tile, icc, lat, lon ], dims=('tile','icc','lat', 'lon'))
    orig['stemmass_ns'] = arrstemmass_ns
    orig['stemmass_ns'].attrs['long_name']='Non-structural stem mass'
    orig['stemmass_ns'].attrs['_FillValue']=-999
    orig['stemmass_ns'].attrs['units']='kg C/m^2'

    arrstemmass_s = xr.DataArray(stemmass_s, coords=[tile, icc, lat, lon ], dims=('tile','icc','lat', 'lon'))
    orig['stemmass_s'] = arrstemmass_s
    orig['stemmass_s'].attrs['long_name']='Structural stem mass'
    orig['stemmass_s'].attrs['_FillValue']=-999
    orig['stemmass_s'].attrs['units']='kg C/m^2'

    arrrootmass_ns = xr.DataArray(rootmass_ns, coords=[tile, icc, lat, lon ], dims=('tile','icc','lat', 'lon'))
    orig['rootmass_ns'] = arrrootmass_ns
    orig['rootmass_ns'].attrs['long_name']='Non-structural root mass'
    orig['rootmass_ns'].attrs['_FillValue']=-999
    orig['rootmass_ns'].attrs['units']='kg C/m^2'

    arrrootmass_s = xr.DataArray(rootmass_s, coords=[tile, icc, lat, lon ], dims=('tile','icc','lat', 'lon'))
    orig['rootmass_s'] = arrrootmass_s
    orig['rootmass_s'].attrs['long_name']='Structural root mass'
    orig['rootmass_s'].attrs['_FillValue']=-999
    orig['rootmass_s'].attrs['units']='kg C/m^2'

    # Add peat soil C
    # Choose a variable with the right dimensions, it will be overwritten with 0 anyway since this is a non-peatland
    peatC = orig['SDEP'] #np.zeros((len(tile),len(lat),len(lon)))
    peatC = 0 * peatC
    peatSoilC = xr.DataArray(peatC, coords=[tile, lat, lon ], dims=('tile','lat', 'lon'))
    orig['peatSoilC'] = peatC #peatSoilC * 100.
    orig['peatSoilC'].attrs['long_name']='Peat soil C for peatland tiles'
    orig['peatSoilC'].attrs['_FillValue']=-999
    orig['peatSoilC'].attrs['units']='kg C/m^2'

    # Now make the litter and soil C be per layer:

    ignd = orig.layer
    iccp2num = int(len(iccp2)) 
    litat = orig['litrmass'].attrs
    scat = orig['soilcmas'].attrs
    litzero = np.zeros((len(tile),len(ignd),iccp2num,len(lat),len(lon)))

    lit = xr.DataArray(litzero, coords=[tile, layer, iccp2, lat, lon ], dims=('tile','layer','iccp2', 'lat', 'lon'))
    sc = xr.DataArray(litzero, coords=[tile, layer, iccp2, lat, lon ], dims=('tile','layer','iccp2', 'lat', 'lon'))
    lit.attrs = litat
    sc.attrs = scat
    orig['litrmass'] = lit
    orig['soilcmas'] = sc
    
    # Add in the N variables:       
              # ngleafmasrow   = ncGet3DVar(initid, 'ngleafmas', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # ngleafmas_NSrow = ncGet3DVar(initid, 'ngleafmas_ns', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # ngleafmassrow  = ncGet3DVar(initid, 'ngleafmas_s', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # nbleafmasrow   = ncGet3DVar(initid, 'nbleafmas', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # nstemmassrow   = ncGet3DVar(initid, 'nstemmass', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # nstemmass_NSrow = ncGet3DVar(initid, 'nstemmass_ns', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # nstemmasssrow  = ncGet3DVar(initid, 'nstemmass_s', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # nrootmassrow   = ncGet3DVar(initid, 'nrootmass', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # nrootmass_NSrow = ncGet3DVar(initid, 'nrootmass_ns', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # nrootmasssrow  = ncGet3DVar(initid, 'nrootmass_s', start = [lonIndex, latIndex, 1, 1], count = [1, 1, icc, nmos], format = [nlat, nmos,icc])
              # nh4_massrow    = ncGet3DVar(initid, 'nh4_mass', start = [lonIndex, latIndex, 1, 1,1], count = [1, 1, iccp1, nmos], format = [nlat, nmos,iccp1])
              # no3_massrow    = ncGet3DVar(initid, 'no3_mass', start = [lonIndex, latIndex, 1, 1,1], count = [1, 1, iccp1, nmos], format = [nlat, nmos,iccp1])
              # nlitrmassrow   = ncGet3DVar(initid, 'nlitrmass', start = [lonIndex, latIndex, 1, 1, 1], count = [1, 1, iccp1, nmos], format = [nlat, nmos, iccp1])
              # soilnmasrow    = ncGet3DVar(initid, 'soilnmas', start = [lonIndex, latIndex, 1, 1,1], count = [1, 1, iccp1, nmos], format = [nlat, nmos,iccp1])
              # soilPHrow      = ncGet2DVar(initid, 'soilPH', start = [lonIndex, latIndex, 1], count = [1, 1, nmos], format = [nlat, nmos])


    # Write this to netcdf. We keep the same name so it is not necessary to update all the 
    # joboptions files.
    orig.to_netcdf(filename)
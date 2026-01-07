#!/space/hall5/sitestore/eccc/crd/ccrp/users/mib001/conda/envs/py3-cforge
import rioxarray
import xarray
import matplotlib.pyplot as plt
import numpy as np

#S.R. Curasi June 2023
#This script takes classic output files and applies circular shifts to 
#the timesteps to account for the earth's rotation
#it was originally used to make an inversion prior using CLASSIC but may have other applications

#user input
base_path='/path/path/'
input_file='in.nc'
output_file='out.nc'
var_name='nepCMIP'
timestep=3 #timestep of the run files in hours or fractions of hours
run_typ_rot=False #True or False if the run is on a rotated grid

#read and check the data
data = xarray.open_dataset(base_path+input_file,engine='netcdf4')
#print(data)

#replicate the calculations done in modeModule.f90
#to calculate the time zone offsets this was written
if(run_typ_rot==True):
    lon=data['longitude'][:,:]+(180*2)
else:
    lon=data['longitude'][:]

    #these would be used for the encis inversion grid
    #lon=data['longitude'][:]+180
    #lon=np.roll(lon, 180, axis=0)

timeZone = (np.round(lon/15*60))/(1800/60)
metInputTimeStep=21600 #time step of the met input in seconds (0.25 of a day)
numberMetInputinDay=int(86400/metInputTimeStep)
numberPhysInMet=round(metInputTimeStep/1800)
shortSteps=numberPhysInMet*numberMetInputinDay
timeZone = np.where((lon >= 180) | (lon < 0),shortSteps - timeZone,-timeZone)
timeZone = timeZone/((timestep*60)/30) #adjust the offset for non-half-hourly files
timeZone = timeZone.astype(int)

#used for testing
# print(timeZone*0.5)
# pshif=data['gpp'][0,:,:]
# for i in range(len(lon)):
#     pshif[:,i] = pshif[:,i]*0+timeZone[i]
# (pshif*0.5).plot(levels=24,cmap="husl")

#plt.imshow( timeZone*0.5 , cmap = 'magma' )
#print((timeZone*0.5).mean())
#plt.colorbar()

#apply circular shifts as done in modeModule.f90
if(run_typ_rot==True):
    for i in range(len(data['lon'])):
        for j in range(len(data['lat'])):
            print("index = "+str(i)+" "+str(j)+" lon = "+str(data['longitude'][j,i].values)+" shift(steps) = " +str(timeZone[j,i])+" shift(hrs) = "+str(timeZone[j,i]*timestep))
            data[var_name][:,j,i] = np.roll(data[var_name][:,j,i], timeZone[j,i], axis=0)  
else:
    for i in range(len(lon)):
        print("index = "+str(i)+" lon band = "+str(data['longitude'][i].values)+" shift(steps) = " +str(timeZone[i])+" shift(hrs) = "+str(timeZone[i]*timestep))
        data[var_name][:,:,i] = np.roll(data[var_name][:,:,i], timeZone[i], axis=0)

data.to_netcdf(path=base_path+output_file, mode='w')


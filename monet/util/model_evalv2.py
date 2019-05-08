#model_evalv2.py
#Using MONET tools to remap hysplit and volcat data to same grid
#Will ideally execute statistical routines to compare the datasets

import monet
import xarray as xr
import numpy as np
import numpy.ma as ma
from numpy import sort
from glob import glob
import pandas as pd
from monet.util import tools
from monet.models import hysplit
import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from netCDF4 import Dataset
import volcat

hr=['19']

#Reading in HYSPLIT data
direct='/pub/ECMWF/JPSS/reventador/'
fname1 = 'cdump.004b'
flist = [direct+fname1]
d1=datetime.datetime(2019, 2, 25, 19, 0)
d2=datetime.datetime(2019, 2, 25, 20, 0)
#All HYSPLIT data
for fname in flist:
    print(fname)
    hxr = hysplit.open_dataset(fname,drange=[d1,d2])
    alts=hxr.attrs['Level top heights (m)']
    hlat=hxr.coords['latitude']
    hlon=hxr.coords['longitude']
    time=hxr.coords['time']
#Adding 4 variable types together
par1=hxr.p006
par2=hxr.par2
par3=hxr.par3
par4=hxr.par4
total_par=par1+par2+par3+par4
#altitude array (z coordinate) 
alt_str=alts.astype(str)

######Creating top height array
#Creating array of 0 and 1 (0 where nan, 1 where value)
heights=(total_par > 0).astype(int)
#Multiplying each level of total particulate by the
#altitude (z) to yield cloud top heights (km)
heights[:,0,:,:]=heights[:,0,:,:] * (alts[0]/1000)
heights[:,1,:,:]=heights[:,1,:,:] * (alts[1]/1000)
heights[:,2,:,:]=heights[:,2,:,:] * (alts[2]/1000)
heights[:,3,:,:]=heights[:,3,:,:] * (alts[3]/1000)
heights[:,4,:,:]=heights[:,4,:,:] * (alts[4]/1000)
#Creating top height array by selecting the max value along the z axis
top_height=np.amax(heights,axis=1)
#Assigning nan to 0 values in top_height
top_height=top_height.where(top_height != 0)
#Creating mask of top_height array
#mask_hysp=np.isnan(top_height)
#Creating 2-D HYSPLIT array
HYSP_height=top_height[0,:,:]

#Creating netcdf file of HYSPLIT output
#hxr.to_netcdf('file.nc')

#Reading in VOLCAT files
directory="/pub/ECMWF/JPSS/volcat/reventador/"
i=0
while i < len(hr):
    print(i)
    fname2='geocatL2.GOES-16.Full_Disk.2019056.'+hr[i]+'*.hdf'
    files=sort(glob(directory+fname2))
    das=[]
#Reading in multiple files
    for j in files:
        das.append(volcat.open_dataset(j))
        dset=xr.concat(das,dim='time')
#Making mask array for plotting
    missing=dset.vcat_ashprop_13_15_16_ash_top_height._FillValue
#Extracting ash top height
    mask_volc=dset.vcat_ashprop_13_15_16_ash_top_height.where(dset.vcat_ashprop_13_15_16_ash_top_height != missing)
#Creating 1hr average over time dimension
    VOLC_height=mask_volc.mean(axis=0)
    vlat=dset.latitude
    vlon=dset.longitude
    i += 1

#Regridding  HYSPLIT array to VOLCAT latlon grid
#remapped = VOLC_height.monet.remap_xesmf(HYSP_height,method='nearest_s2d', reuse_weights=True)
#remapped = VOLC_height.monet.remap_data(HYSP_height,radius_of_influence=1e3)
#q = VOLC_height.rename({"lines":'y',"elements":'x'})
remapped = VOLC_height.monet.remap_nearest(HYSP_height, radius_of_influence=1e3)
mask_HYSP = remapped.where(remapped != 0.)

#Plotting data at various HYSPLIT levels
fig=plt.figure('Ash_Top_Height')
img_proj=ccrs.PlateCarree()
m=plt.axes(projection=img_proj)
m.add_feature(cfeat.LAND)
m.add_feature(cfeat.COASTLINE)
m.add_feature(cfeat.BORDERS)
#Plots volcat data
plt.pcolormesh(vlon, vlat, VOLC_height, cmap='autumn',transform=ccrs.PlateCarree())
plt.colorbar()
#Plots hysplit data
plt.pcolormesh(vlon, vlat, mask_HYSP,cmap='winter')
#plt.pcolormesh(hlon, hlat, HYSP_height,cmap='winter')
plt.colorbar()
#Plots Reventador Volcano symbol
plt.scatter(-77.6558,-0.0775,c='black',marker='^')
#Zoomed in region
plt.xlim(-78, -76.8)
plt.ylim(-0.9,1.1)
plt.title('Volocat & HYSPLIT at hr '+hr[0]+'\n Ash Top Height (km)')
print('saving','Images/VOLCATHYSPLIT_AshTopHeight_hr'+hr[0]+'avg_'+fname1+'_remapped.png')
#monet.plots.savefig('Images/VOLCATHYSPLIT_AshTopHeight_hr'+hr[0]+'avg_'+fname1+'.png')
monet.plots.savefig('Images/VOLCATHYSPLIT_AshTopHeight_hr'+hr[0]+'avg_'+fname1+'_remapped.png')
plt.show()

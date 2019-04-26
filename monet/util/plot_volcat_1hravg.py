#plot_volcat_1hravg.py
#Hard coded plotting of volcat data from Reventador eruption
#Creating 1hr averages from 15min output
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
from glob import glob
from numpy import sort
import pandas as pd
import volcat

hr='21'
directory="/pub/ECMWF/JPSS/volcat/reventador/"
fname="geocatL2.GOES-16.Full_Disk.2019056."+hr+"*.hdf"
files=sort(glob(directory+fname))
das=[]
#Reading in multiple files
for j in files:
    das.append(volcat.open_dataset(j))
dset=xr.concat(das,dim='time')

#Making mask array for plotting
missing=dset.vcat_ashprop_13_15_16_ash_top_height._FillValue

#Extracting average ash top height, ash effective radius, ash mass loading
masked_height=ma.masked_equal(dset.vcat_ashprop_13_15_16_ash_top_height, missing)
masked_radius=ma.masked_equal(dset.vcat_ashprop_13_15_16_ash_effective_radius, missing)
masked_mass=ma.masked_equal(dset.vcat_ashprop_13_15_16_ash_mass_loading, missing)

#Creating 1hr average over time dimension
height_avg=masked_height.mean(axis=0)
radius_avg=masked_radius.mean(axis=0)
mass_avg=masked_mass.mean(axis=0)
lattit=dset.pixel_latitude[:,:] * dset.pixel_latitude.scale_factor
longit=dset.pixel_longitude[:,:] * dset.pixel_longitude.scale_factor
lat=lattit.mean(axis=0)
lon=longit.mean(axis=0)

#Plotting data
fig=plt.figure('Ash_Top_Height')
img_proj=ccrs.PlateCarree()
m=plt.axes(projection=img_proj)
m.add_feature(cfeat.LAND)
m.add_feature(cfeat.COASTLINE)
m.add_feature(cfeat.BORDERS)
plt.pcolormesh(lon, lat, height_avg,cmap='autumn')
plt.colorbar()
plt.scatter(-77.6558,-0.0775,c='black',marker='^')
plt.xlim(-78.5, -77)
plt.ylim(-1,0.5)
plt.title('1-hr Avg Ash Top Height (km)')
plt.savefig('Images/Ash_top_height_hr'+hr+'avg.png')
#plt.show()

fig=plt.figure('Ash_Effective_radius')
img_proj=ccrs.PlateCarree()
m=plt.axes(projection=img_proj)
m.add_feature(cfeat.LAND)
m.add_feature(cfeat.COASTLINE)
m.add_feature(cfeat.BORDERS)
plt.pcolormesh(lon, lat, radius_avg, cmap='autumn')
plt.colorbar()
plt.scatter(-77.6558,-0.0775,c='black',marker='^')
plt.xlim(-78.5, -77)
plt.ylim(-1,0.5)
plt.title('1-hr Avg Ash Effective Radius (um)')
plt.savefig('Images/Ash_effective_radius_hr'+hr+'avg.png')
#plt.show()

fig=plt.figure('Ash_Mass_Loading')
img_proj=ccrs.PlateCarree()
m=plt.axes(projection=img_proj)
m.add_feature(cfeat.LAND)
m.add_feature(cfeat.COASTLINE)
m.add_feature(cfeat.BORDERS)
plt.pcolormesh(lon, lat, mass_avg, cmap='autumn')
plt.colorbar()
plt.scatter(-77.6558,-0.0775,c='black',marker='^')
plt.xlim(-78.5, -77)
plt.ylim(-1,0.5)
plt.title('1-hr Avg Ash Mass Loading (ug/m^2)')
plt.savefig('Images/Ash_mass_loading_hr'+hr+'avg.png')
#plt.show()

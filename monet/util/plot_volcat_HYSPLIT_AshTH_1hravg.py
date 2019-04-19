#plot_volcat_HYSPLIT_AshTH_1hravg.py
#Hard coded plotting of volcat data and HYSPLIT data from Reventador eruption
#Creating 1hr averages from 15min output
import xarray as xr
import matplotlib.pyplot as plt
import textwrap
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
from glob import glob
from numpy import sort
import pandas as pd
import volcat
import tools
import hysplit
import datetime
import seaborn as sns

hr=['17','18','19','20','21']
directory="/pub/ECMWF/JPSS/volcat/reventador/"
#Reading in HYSPLIT data
direct='/pub/ECMWF/JPSS/reventador/'
fname1 = 'cdump.emt2'
flist = [direct+fname1]

#All data
for fname in flist:
    print(fname)
    hxr = hysplit.open_dataset(fname)
    alts=hxr.attrs['Level top heights (m)']
    hlat=hxr.coords['latitude']
    hlon=hxr.coords['longitude']
    time=hxr.coords['time']

#Adding 4 variable types together
par1=hxr.p006
#par2=hxr.par2
#par3=hxr.par3
#par4=hxr.par4
total_par=par1#+par2+par3+par4

#altitude array (z coordinate) 
alt_str=alts.astype(str)
######Creating top height array
#Creating array of 0 and 1 (0 where nan, 1 where value)
heights=(total_par > 0).astype(int)
#Multiplying each level of heights by the altitude (z) value (km)
#heights[:,0,:,:]=heights[:,0,:,:] * (alts[0]/1000)
#heights[:,1,:,:]=heights[:,1,:,:] * (alts[1]/1000)
#heights[:,2,:,:]=heights[:,2,:,:] * (alts[2]/1000)
#heights[:,3,:,:]=heights[:,3,:,:] * (alts[3]/1000)
#heights[:,4,:,:]=heights[:,4,:,:] * (alts[4]/1000)
heights[:,0,:,:]=heights[:,0,:,:] * (alts[2]/1000)
heights[:,1,:,:]=heights[:,1,:,:] * (alts[3]/1000)
heights[:,2,:,:]=heights[:,2,:,:] * (alts[4]/1000)
#Creating top height array by selecting the max value along the z axis
top_height=np.amax(heights,axis=1)
#Assigning nan to 0 values in top_height
top_height=top_height.where(top_height != 0)

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

#Extracting ash mass loading
    masked_height=ma.masked_equal(dset.vcat_ashprop_13_15_16_ash_top_height, missing)

#Creating 1hr average over time dimension
    height_avg=masked_height.mean(axis=0)
    lattit=dset.pixel_latitude[:,:] * dset.pixel_latitude.scale_factor
    longit=dset.pixel_longitude[:,:] * dset.pixel_longitude.scale_factor
    lat=lattit.mean(axis=0)
    lon=longit.mean(axis=0)

#Plotting data at various HYSPLIT levels
    fig=plt.figure('Ash_Top_Height')
    img_proj=ccrs.PlateCarree()
    m=plt.axes(projection=img_proj)
    m.add_feature(cfeat.LAND)
    m.add_feature(cfeat.COASTLINE)
    m.add_feature(cfeat.BORDERS)
#Plots volcat data
    plt.pcolormesh(lon, lat, height_avg, cmap='autumn',transform=ccrs.PlateCarree())
    plt.colorbar()
#Plots hysplit data
    plt.pcolormesh(hlon, hlat, top_height[i,:,:],cmap='winter')
    plt.colorbar()
#Plots Reventador Volcano symbol
    plt.scatter(-77.6558,-0.0775,c='black',marker='^')
#Zoomed in region
    plt.xlim(-78, -76.8)
    plt.ylim(-0.9,01.1)
    plt.title('Volocat & HYSPLIT at hr '+hr[i]+'\n Ash Top Height (km)')
    plt.savefig('Images/VOLCATHYSPLIT_AshTopHeight_hr'+hr[i]+'avg_'+fname1+'.png')
    plt.show()
    plt.close()
    i += 1

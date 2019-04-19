#plot_volcat_HYSPLIT_1hravg.py
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

hr=['17','18','19','20']
directory="/pub/ECMWF/JPSS/volcat/reventador/"
#Reading in HYSPLIT data
direct='/pub/ECMWF/JPSS/reventador/'
fname1 = 'cdump.003'
flist = [direct+fname1]

d1 = datetime.datetime(2019,2,25,17)
d2 = datetime.datetime(2019,2,26,3)

#All data
for fname in flist:
    print(fname)
    hxr = hysplit.open_dataset(fname, drange=[d1,d2])
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

#Hours: 18:22 (5 hours)
#Z: 2000, 4000, 6000, 8000, 10000 (5 levels)
alt_str=alts.astype(str)
#Creating arrays based on levels and multiplying by the depth
tot_lvl0=total_par[:,0,:,:] * alts[0]
tot_lvl1=total_par[:,1,:,:] * alts[1]
tot_lvl2=total_par[:,2,:,:] * alts[2]
tot_lvl3=total_par[:,3,:,:] * alts[3]
tot_lvl4=total_par[:,4,:,:] * alts[4]

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

#Making mask array for plotting - fill value is the same across variables
    missing=dset.vcat_ashprop_13_15_16_ash_top_height._FillValue

#Extracting ash mass loading
    masked_mass=ma.masked_equal(dset.vcat_ashprop_13_15_16_ash_mass_loading, missing)

#Creating 1hr average over time dimension
    mass_avg=masked_mass.mean(axis=0)
    lattit=dset.pixel_latitude[:,:] * dset.pixel_latitude.scale_factor
    longit=dset.pixel_longitude[:,:] * dset.pixel_longitude.scale_factor
    lat=lattit.mean(axis=0)
    lon=longit.mean(axis=0)

#Plotting data at various HYSPLIT levels
    fig=plt.figure('Ash_Mass_Loading')
    img_proj=ccrs.PlateCarree()
    m=plt.axes(projection=img_proj)
    m.add_feature(cfeat.LAND)
    m.add_feature(cfeat.COASTLINE)
    m.add_feature(cfeat.BORDERS)
#Plots Reventador Volcano symbol
    plt.scatter(-77.6558,-0.0775,c='black',marker='^')
#Plots volcat data
    plt.pcolormesh(lon, lat, mass_avg, transform=ccrs.PlateCarree())
    plt.colorbar()
#Plots hysplit data
    plt.pcolormesh(hlon, hlat, tot_lvl0[i,:,:],cmap='seismic')
    plt.colorbar()
#Zoomed in region
    plt.xlim(-78, -76.8)
    plt.ylim(-0.9,0.9)
    plt.title('Volocat & HYSPLIT at hr '+hr[i]+'\n Ash Mass Loading (g/m^2) \n at '+alt_str[0]+'m')
    plt.savefig('Images/VOLCATHYSPLIT_Ash_mass_loading_hr'+hr[i]+'avg_'+alt_str[0]+'m.png')
    plt.show()
    #plt.close()

    fig=plt.figure('Ash_Mass_Loading')
    img_proj=ccrs.PlateCarree()
    m=plt.axes(projection=img_proj)
    m.add_feature(cfeat.LAND)
    m.add_feature(cfeat.COASTLINE)
    m.add_feature(cfeat.BORDERS)
#Plots Reventador Volcano symbol
    plt.scatter(-77.6558,-0.0775,c='black',marker='^')
#Plots volcat data
    plt.pcolormesh(lon, lat, mass_avg, transform=ccrs.PlateCarree())
    plt.colorbar()
#Plots hysplit data
    plt.pcolormesh(hlon, hlat, tot_lvl1[i,:,:],cmap='seismic')
    plt.colorbar()
#Zoomed in region
    #plt.xlim(-78.5, -77)
    #plt.ylim(-1,0.5)
    plt.xlim(-78, -76.8)
    plt.ylim(-0.9,0.9)
    plt.title('Volocat & HYSPLIT at hr '+hr[i]+'\n Ash Mass Loading (ug/m^2) \n at '+alt_str[1]+'m')
    plt.savefig('Images/VOLCATHYSPLIT_Ash_mass_loading_hr'+hr[i]+'avg_'+alt_str[1]+'m.png')
    plt.show()
    #plt.close()

    fig=plt.figure('Ash_Mass_Loading')
    img_proj=ccrs.PlateCarree()
    m=plt.axes(projection=img_proj)
    m.add_feature(cfeat.LAND)
    m.add_feature(cfeat.COASTLINE)
    m.add_feature(cfeat.BORDERS)
#Plots Reventador Volcano symbol
    plt.scatter(-77.6558,-0.0775,c='black',marker='^')
#Plots volcat data
    plt.pcolormesh(lon, lat, mass_avg, transform=ccrs.PlateCarree())
    plt.colorbar()
#Plots hysplit data
    plt.pcolormesh(hlon, hlat, tot_lvl2[i,:,:],cmap='seismic')
    plt.colorbar()
#Zoomed in region
    plt.xlim(-78, -76.8)
    plt.ylim(-0.9,0.9)
    plt.title('Volocat & HYSPLIT at hr '+hr[i]+'\n Ash Mass Loading (ug/m^2) \n at '+alt_str[2]+'m')
    plt.savefig('Images/VOLCATHYSPLIT_Ash_mass_loading_hr'+hr[i]+'avg_'+alt_str[2]+'m.png')
    plt.show()
    #plt.close()

    fig=plt.figure('Ash_Mass_Loading')
    img_proj=ccrs.PlateCarree()
    m=plt.axes(projection=img_proj)
    m.add_feature(cfeat.LAND)
    m.add_feature(cfeat.COASTLINE)
    m.add_feature(cfeat.BORDERS)
#Plots Reventador Volcano symbol
    plt.scatter(-77.6558,-0.0775,c='black',marker='^')
#Plots volcat data
    plt.pcolormesh(lon, lat, mass_avg, transform=ccrs.PlateCarree())
    plt.colorbar()
#Plots hysplit data
    plt.pcolormesh(hlon, hlat, tot_lvl3[i,:,:],cmap='seismic')
    plt.colorbar()
#Zoomed in region
    plt.xlim(-78, -76.8)
    plt.ylim(-0.9,0.9)
    plt.title('Volocat & HYSPLIT at hr '+hr[i]+'\n Ash Mass Loading (ug/m^2) \n at '+alt_str[3]+'m')
    plt.savefig('Images/VOLCATHYSPLIT_Ash_mass_loading_hr'+hr[i]+'avg_'+alt_str[3]+'m.png')
    plt.show()
    #plt.close()

    fig=plt.figure('Ash_Mass_Loading')
    img_proj=ccrs.PlateCarree()
    m=plt.axes(projection=img_proj)
    m.add_feature(cfeat.LAND)
    m.add_feature(cfeat.COASTLINE)
    m.add_feature(cfeat.BORDERS)
#Plots Reventador Volcano symbol
    plt.scatter(-77.6558,-0.0775,c='black',marker='^')
#Plots volcat data
    plt.pcolormesh(lon, lat, mass_avg, transform=ccrs.PlateCarree())
    plt.colorbar()
#Plots hysplit data
    plt.pcolormesh(hlon, hlat, tot_lvl4[i,:,:],cmap='seismic')
    plt.colorbar()
#Zoomed in region
    plt.xlim(-78, -76.8)
    plt.ylim(-0.9,0.9)
    plt.title('Volocat & HYSPLIT at hr '+hr[i]+'\n Ash Mass Loading (ug/m^2) \n at '+alt_str[4]+'m')
    plt.savefig('Images/VOLCATHYSPLIT_Ash_mass_loading_hr'+hr[i]+'avg_'+alt_str[4]+'m.png')
    plt.show()
    #plt.close()
    i+= 1

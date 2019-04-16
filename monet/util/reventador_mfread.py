#Reading in multiple volcat files
import volcat
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat

flist=["geocatL2.GOES-16.Full_Disk.2019056.170030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.171530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.173030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.174530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.180030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.181530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.183030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.184530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.190030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.191530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.193030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.194530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.200030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.201530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.203030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.204530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.210030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.211530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.213030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.214530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.220030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.221530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.223030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.224530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.230030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.231530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.233030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.234530.hdf"]

flist2=["geocatL2.GOES-16.Full_Disk.2019056.190030.hdf",
        "geocatL2.GOES-16.Full_Disk.2019056.191530.hdf",
        "geocatL2.GOES-16.Full_Disk.2019056.193030.hdf",
        "geocatL2.GOES-16.Full_Disk.2019056.194530.hdf"]

fname="geocatL2.GOES-16.Full_Disk.2019056.190030.hdf"

directory="/pub/ECMWF/JPSS/volcat/reventador/"

#vdata=volcat.open_mfdataset(directory+'geocatL2.GOES-16.Full_Disk.2019056.19*.hdf')
vdata=volcat.open_mfdataset(directory+'geocatL2.GOES-16.Full_Disk.2019056.*.hdf')
#vdata=volcat.open_dataset(directory+fname)
#volcat.plot_height(vdata)
#volcat.plot_radius(vdata)
#volcat.plot_mass(vdata)

#Extract desired variables
height=volcat.get_height(vdata)
height_test=height.resample(time='1H').mean()
#radius=volcat.get_radius(vdata)
#mass=volcat.get_mass_loading(vdata)

#Create means of masked vdata
#mean_height=np.mean(height,axis=0)
#mean_radius=np.mean(radius,axis=0)
#mean_mass=np.mean(mass,axis=0)
#mean_mass=mass.resample(time='1H').mean()

#Create lat lon array for plotting means
#lat = vdata.pixel_latitude[:,:] * vdata.pixel_latitude.scale_factor
#lon = vdata.pixel_longitude[:,:] * vdata.pixel_longitude.scale_factor
#mean_lat=np.mean(lat,axis=0)
#mean_lon=np.mean(lon,axis=0)

#Create Plots of mean data
#fig = plt.figure('Ash_Top_Height_mean')
#m = plt.axes(projection=ccrs.PlateCarree())
#m.add_feature(cfeat.LAND)
#m.add_feature(cfeat.COASTLINE)
#m.add_feature(cfeat.BORDERS)
#plt.pcolormesh(mean_lon, mean_lat, mean_height, transform=ccrs.PlateCarree())
#plt.colorbar()
#plt.title('Mean Ash Top Height (km)')
#plt.show()

#for i in flist2:
#    volcat.plot_height(vdata,height[i,:,:])

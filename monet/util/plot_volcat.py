#plot_volcat.py
#Hard coded plotting of volcat data from Reventador eruption
#15min output files
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import volcat

hr=['17','18','19','20','21','22']
mm=['00','15','30','45']
directory="/pub/ECMWF/JPSS/volcat/reventador/"

i=0
while i < len(hr):
    j=0
    while j < len(mm):
        fname='geocatL2.GOES-16.Full_Disk.2019056.'+hr[i]+mm[j]+'30.hdf'
        dset=volcat.open_dataset(directory+fname)
        print(fname)
        #Identifying fill value for array
        missing=dset.vcat_ashprop_13_15_16_ash_top_height._FillValue

        #Extracting ash mass loading
        masked_height=ma.masked_equal(dset.vcat_ashprop_13_15_16_ash_top_height, missing)

        #Creating lat and lon arrays 
        lat=dset.pixel_latitude[:,:] * dset.pixel_latitude.scale_factor
        lon=dset.pixel_longitude[:,:] * dset.pixel_longitude.scale_factor
    
        #Plotting data
        fig=plt.figure('Ash_Top_Height')
        img_proj=ccrs.PlateCarree()
        m=plt.axes(projection=img_proj)
        m.add_feature(cfeat.LAND)
        m.add_feature(cfeat.COASTLINE)
        m.add_feature(cfeat.BORDERS)
        #Plots volcat data
        plt.pcolormesh(lon, lat, masked_height, vmin=3.5, vmax=10, cmap='autumn', transform=ccrs.PlateCarree())
        plt.colorbar()
        #Plots Reventador Volcano symbol
        plt.scatter(-77.6558,-0.0775,c='black',marker='^')
        #Zoomed in region
        plt.xlim(-78.2, -76.8)
        plt.ylim(-1.,0.5)
        plt.title('Ash Top Height (km) \n from VOLCAT at '+hr[i]+':'+mm[j])
        plt.savefig('Images/VOLCAT_AshTopHeight_'+hr[i]+mm[j]+'30.png')
       #plt.show()
        plt.close()
        j += 1
    i += 1

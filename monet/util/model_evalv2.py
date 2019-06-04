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
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from netCDF4 import Dataset
import volcat
import csi
import csi_updated
from itertools import permutations

#Designate hour for VOLCAT files
#hr = ['17', '18', '19', '20', '21']
hr = ['17']
i = 0
#Statistics to be calculated
CSI = []
FSS = []
FSS2 = []
directory="/pub/ECMWF/JPSS/volcat/reventador/"
direct = '/pub/ECMWF/JPSS/reventador/'
fname1 = 'cdump.004b'
#Looping over all hours
while i < len(hr):
    print(i)
    hr2 = int(hr[i])
    hr3 = hr2 + 1
    d1 = datetime.datetime(2019, 2, 25, hr2, 0)
    d2 = datetime.datetime(2019, 2, 25, hr3, 0)
    #HYSPLIT    
    flist = [direct+fname1]
    for fname in flist:
        print(fname)
        hxr = hysplit.open_dataset(fname,drange=[d1,d2])
        alts = hxr.attrs['Level top heights (m)']
        hlat = hxr.coords['latitude']
        hlon = hxr.coords['longitude']
        time = hxr.coords['time']
    #Adding 4 variable types together
    par1 = hxr.p006
    par2 = hxr.par2
    par3 = hxr.par3
    par4 = hxr.par4
    total_par = par1+par2+par3+par4
    #altitude array (z coordinate) 
    alt_str = alts.astype(str)

    ##### Creating top height array #####
    #Creating array of 0 and 1 (0 where nan, 1 where value)
    heights = (total_par > 0).astype(int)
    #Multiplying each level of total particulate by the
    #altitude (z) to yield cloud top heights (km)
    heights[:,0,:,:] = heights[:,0,:,:] * (alts[0]/1000)
    heights[:,1,:,:] = heights[:,1,:,:] * (alts[1]/1000)
    heights[:,2,:,:] = heights[:,2,:,:] * (alts[2]/1000)
    heights[:,3,:,:] = heights[:,3,:,:] * (alts[3]/1000)
    heights[:,4,:,:] = heights[:,4,:,:] * (alts[4]/1000)
    #Creating top height array by selecting the max value along the z axis
    top_height = np.amax(heights,axis=1)
    #Assigning nan to 0 values in top_height
    top_height = top_height.where(top_height != 0)
    #Creating 2-D HYSPLIT array
    HYSP_height = top_height[0,:,:]

    #Reading in VOLCAT files
    fname2 = 'geocatL2.GOES-16.Full_Disk.2019056.'+hr[i]+'*.hdf'
    files = sort(glob(directory+fname2))
    das = []
    #Reading in multiple files
    for j in files:
        das.append(volcat.open_dataset(j))
        dset = xr.concat( das, dim = 'time')
    #Making mask array for plotting
    missing = dset.vcat_ashprop_13_15_16_ash_top_height._FillValue
    #Extracting ash top height
    mask_volc = dset.vcat_ashprop_13_15_16_ash_top_height.where( dset.vcat_ashprop_13_15_16_ash_top_height != missing)
    #Creating 1hr average over time dimension
    VOLC_height = mask_volc.mean(axis=0)
    vlat = dset.latitude
    vlon = dset.longitude
    #Regridding  HYSPLIT array to VOLCAT latlon grid
    remapped = VOLC_height.monet.remap_nearest( HYSP_height, radius_of_influence = 1e3)
    mask_HYSP = remapped.where( remapped != 0.)

    ### STATISTICS CALCULATIONS FOLLOW ###
    #Calculating CSI for all hours 
    #CSI.append(csi.calc_csi(VOLC_height, mask_HYSP, threshold = 0, verbose = 0))

    #Calculating fss
    #FSS.append(csi.calc_fss(VOLC_height, mask_HYSP, threshold = 0, verbose = 0, sn = [0,3]))
    FSS2.append(csi_updated.calc_fss(VOLC_height, mask_HYSP, threshold = 0, verbose = 0, sn = [0,3]))
    #Creating binary arrays
    #Vma = csi.mask_threshold(VOLC_height)
    #Hma = csi.mask_threshold(mask_HYSP)

    #Vtrim, Htrim, Vtrimboth, Htrimboth = csi.trim_arrays(Vma, Hma, sn=[0,3])
    

    i += 1  #End of loop for each hour

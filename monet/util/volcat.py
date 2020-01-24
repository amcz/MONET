#volcat.py
#A reader for VOLCAT data using xarray
#For use with MONET
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import numpy as np
import numpy.ma as ma
import pandas as pd

"""
This script contains routines that open/read VOLCAT data in xarray format, 
manipulate arrays as necessary, and plots desirable variables.
-------------
Functions:
-------------
open_dataset: opens single VOLCAT file
open_mfdataset: opens multiple VOLCAT files
_get_time: set time dimension for VOLCAT data
_get_latlon: rename lat/lon, set coordinates of VOLCAT data
get_height: returns array of ash top height from VOLCAT
get_radius: returns array of ash effective radius from VOLCAT
get_mass:  returns array of ash mass loading from VOLCAT
get_atherr: returns array of ash top height error from VOLCAT
plot_height: plots ash top height from VOLCAT
plot_radius: plots ash effective radius from VOLCAT
plot_mass: plots ash mass loading from VOLCAT
------------
"""

def open_dataset(fname):
      """Opens single VOLCAT file"""
      print(fname)
      dset = xr.open_dataset(fname,mask_and_scale=False,decode_times=False)
      dset = dset.rename({"Dim1":'y',"Dim0":'x'})
      dset = _get_latlon(dset)
      dset = _get_time(dset)
      return dset

def open_mfdataset(fname):
      """Opens multiple VOLCAT files"""
      print(fname)
      dset = xr.open_mfdataset(fname,concat_dim='time',decode_times=False,mask_and_scale=False)
      from glob import glob
      from numpy import sort
      files = sort(glob(fname))
      das = []
      for i in files:
            das.append(open_dataset(i))
      dset = xr.concat(das,dim='time') 
      dset = _get_latlon(dset)
      dset = dset.rename({"lines":'y',"elements":'x'})
      return dset

def _get_latlon(dset):
      dset = dset.set_coords(['latitude', 'longitude'])
      return dset

def _get_time(dset):
      import pandas as pd
      temp = dset.attrs['time_coverage_start']
      time = pd.to_datetime(temp)
      dset['time'] = time
      dset = dset.expand_dims(dim='time')
      dset = dset.set_coords(['time'])
      return dset

#Extracting variables  
def get_height(dset):      
      """Returns array with retrieved height of the highest layer of ash."""
      """Default units are km above sea-level"""
      height = dset.ash_cth
      h_missing = height._FillValue
      height = height.where(height != h_missing, drop=True)
      return height

def get_radius(dset):
      """Returns 2d array of ash effective radius"""
      """Default units are micrometer"""
      radius = dset.ash_r_eff
      radius = radius.where(radius != radius._FillValue, drop=True)
      return radius

def get_mass(dset):
      """Returns 2d array of ash mass loading"""
      """Default units are grams / meter^2"""
      mass = dset.ash_mass
      mass = mass.where(mass != mass._FillValue, drop=True)
      return mass

def mass_sum(dset):
      mass = get_mass(dset)
      mass2 = mass.where(mass > 0., 0.0).values
      mass_sum = np.sum(mass2)
      return mass_sum

def get_time(dset):
      time = dset.time_coverage_start
      return time

def get_atherr(dset):      
      """Returns array with uncertainty in ash top height from VOLCAT."""
      """Default units are km above sea-level"""
      height_err = dset.ash_cth_uncertainty
      height_err = height_err.where(height_err != height_err._FillValue, drop=True)
      return height_err

#Trim VOLCAT array
def trim_arrray(dset):
      """Trim the VOLCAT array around data
      Make smaller for comparison to HYSPLIT """

#Plotting variables
def plot_height(dset):
      """Plots ash top height from VOLCAT
      Does not save figure - quick image creation"""
      fig = plt.figure('Ash_Top_Height')
      lat=dset.latitude
      lon=dset.longitude
      masked_height=get_height(dset)
      m = plt.axes(projection=ccrs.PlateCarree())
      m.add_feature(cfeat.LAND)
      m.add_feature(cfeat.COASTLINE)
      m.add_feature(cfeat.BORDERS)
      plt.pcolormesh( lon, lat, masked_height, transform=ccrs.PlateCarree())
      plt.colorbar()
      plt.title('Ash Top Height (km)')
      plt.show()

def plot_radius(dset):
      """Plots ash effective radius from VOLCAT
      Does not save figure - quick image creation"""
      fig = plt.figure('Ash_Effective_Radius')
      lat=dset.latitude
      lon=dset.longitude
      masked_radius=get_radius(dset)
      m=plt.axes(projection=ccrs.PlateCarree())
      m.add_feature(cfeat.LAND)
      m.add_feature(cfeat.COASTLINE)
      m.add_feature(cfeat.BORDERS)
      plt.pcolormesh( lon, lat, masked_radius, transform=ccrs.PlateCarree())
      plt.colorbar()
      plt.title('Ash effective radius (um)')
      plt.show()

def plot_mass(dset):
      """Plot ash mass loading from VOLCAT
      Does not save figure - quick image creation"""
      fig = plt.figure('Ash_Mass_Loading')
      lat=dset.latitude
      lon=dset.longitude
      masked_mass=get_mass(dset)
      m=plt.axes(projection=ccrs.PlateCarree())
      m.add_feature(cfeat.LAND)
      m.add_feature(cfeat.COASTLINE)
      m.add_feature(cfeat.BORDERS)
      plt.pcolormesh( lon, lat, masked_mass, transform=ccrs.PlateCarree())
      plt.colorbar()
      plt.title('Ash mass loading (g/m^2)')
      pltshow()
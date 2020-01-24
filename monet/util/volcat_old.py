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
get_mass_loading:  returns array of ash mass loading from VOLCAT
plot_height: plots ash top height from VOLCAT
plot_radius: plots ash effective radius from VOLCAT
plot_mass: plots ash mass loading from VOLCAT
------------
"""

def open_dataset(fname):
      """Opens single VOLCAT file"""
      print(fname)
      dset = xr.open_dataset(fname,mask_and_scale=False,decode_times=False)
      dset = _get_latlon(dset)
      #dset = _get_time(dset)
      #print(dset.latitude.scale_factor)
      dset = dset.rename({"lines":'y',"elements":'x'})
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

def _get_time(dset):
      import pandas as pd
      test = str(dset.Image_Date)[1:] + str(dset.Image_Time)
      time = pd.to_datetime(test,format='%y%j%H%M%S')
      dset = dset.set_coords(['time'])
      dset['time'] = time
      dset = dset.expand_dims(dim='time')
      dset.dims['time'] = time
      return dset

def _get_latlon(dset):
      dset = dset.rename({'pixel_latitude':'latitude'})
      dset = dset.rename({'pixel_longitude':'longitude'})
      dset = dset.set_coords(['latitude', 'longitude'])
      lat = dset.latitude[:,:] * dset.latitude.scale_factor
      lon = dset.longitude[:,:] * dset.longitude.scale_factor
      dset['latitude'] = lat
      dset['longitude'] = lon
      return dset
                  
def get_height(dset):      
      """Returns array with retrieved height of the highest layer of ash."""
      """Default units are km above sea-level"""
      var_header = dset.attrs['Default_Name_ash_ret']
      h_missing = dset[var_header + '_ash_top_height']._FillValue
      masked_height = ma.masked_equal(dset[var_header + '_ash_top_height'],h_missing)
      return masked_height

def get_radius(dset):
      """Returns 2d array of ash effective radius"""
      """Default units are micrometer"""
      var_header = dset.attrs['Default_Name_ash_ret']
      r_missing = dset[var_header + '_ash_effective_radius']._FillValue
      masked_radius = ma.masked_equal(dset[var_header + '_ash_effective_radius'],r_missing)
      return masked_radius

def get_mass_loading(dset):
      """Returns 2d array of ash effective radius"""
      """Default units are micrometer"""
      var_header = dset.attrs['Default_Name_ash_ret']
      m_missing = dset[var_header + '_ash_mass_loading']._FillValue
      masked_mass = ma.masked_equal(dset[var_header + '_ash_mass_loading'],m_missing)
      return masked_mass

def plot_height(dset):
      """Plots ash top height from VOLCAT
      Does not save figure - quick image creation"""
      fig = plt.figure('Ash_Top_Height')
      lat=dset.latitude
      lon=dset.longitude
      #lat = dset.pixel_latitude[:,:] * dset.pixel_latitude.scale_factor
      #lon = dset.pixel_longitude[:,:]*dset.pixel_longitude.scale_factor
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

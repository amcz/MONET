#hysp_func.py
#Functions for manipulating HYSPLIT data
#For use with MONET
"""Functions for manipulating HYSPLIT data.
-------------
Functions:
-------------
hysp_heights: determines ash top height from HYSPLIT
hysp_massload: determines total mass loading from HYSPLIT
"""
from monet.models import hysplit
import xarray as xr
import numpy as np

def hysp_heights(dset):
    """ Calculate ash top-height from HYSPLIT xarray """

    total_par = _add_species(dset) 
    #Create array of 0 and 1 (1 where data exists)
    #NEED TO APPLY ASH THRESHOLD HERE
    heights = (total_par > 0.).astype(int)    
    #Multiply each level by the altitude
    #Yields top heights in (km)
    #alts = total_par.attrs['Level top heights (m)']
    alts = total_par.coords['z']
    x = 0
    while x < (len(alts)):
        heights[:,x,:,:] = heights[:,x,:,:] * (alts[x]/1000.)
        x += 1                 #End of loop calculating heights
    #Determine top height: take max of heights array along z axis
    top_height = heights.max(dim = 'z')
    return top_height

def hysp_massload(dset):
    """ Calculate mass loading from HYSPLIT xarray """
    total_par = _add_species(dset)
    #Replace 0. with nan for plotting purposes
    aml = total_par.where(total_par != 0.)
    return aml
    
def _add_species(dset):
    #Calculate sum of particles
    species = dset.attrs["Species ID"]
    s = 0
    tmp = []
    #Looping through all species in dataset
    while s < len(species):
        tmp.append(dset[species[s]].fillna(0))
        s += 1            #End of loop through species
    total_par = tmp[0]
    p = 1
    #Adding all species together
    while p < len(tmp):
        total_par = total_par + tmp[p]
        p += 1             #End of loop adding all species    
    return total_par

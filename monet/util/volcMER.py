# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
#Calculates mass eruption rate in a variety of ways
#import matplotlib.pyplot  as plt
#from math import *


def mastin(H , DRE=2500):
    """Calculates MER using Mastin equation.
    Input is height (km) . Output is MER (kg/s). Optional input is DRE"""
    p=1/0.241
    dV=(H/2.00)**p
    dM=dV*DRE
    return(dM)

def sparks(H , DRE=2500):
    """Calculates MER using Sparks equation.
    Input is height. Output is MER. Optional input is DRE"""
    p=1/0.259
    dV=(H/1.67)**p
    dM=dV*DRE
    return(dM)

def MER2unit(MER, M63=0.1):
    """MER in kg/s. Output in grams. Assume model output is one unit mass per hour."""
    unit_mass=MER*M63*3600*1000
    #print 'MER %0.3e kg/s , M63 %0.2f , unit mass=%0.3e g.' % (MER, M63 , unit_mass) 
    return unit_mass


def HT2unit(ht, M63=0.1, verbose=True):
    """
    ht  : float  height in kilometers
    M63 : float  mass fraction of fine ash

    RETURNS:
  
    unit_mass : float  
    Assume model output is one unit mass per hour.
    this factor will convert the unit mass to grams. 
    """
    MER = mastin(ht)
    unit_mass = MER*M63*3600*1000
    if verbose: print('HEIGHT %0.1f km,  MER %0.3e kg/s , M63 %0.2f , unit mass=%0.3e g.' %
                                                       (ht, MER, M63 ,
                                                        unit_mass) )
    return unit_mass 



#!/n-home/alicec/anaconda/bin/python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
#from math import *

# import sys
# from scipy.io import netcdf
# from pylab import *
#import numpy as np

# import numpy.ma as ma
#import matplotlib.dates as dates
#import matplotlib.pyplot as plt
#import string
#import datetime

# from matplotlib.path import Path
# import shapely.geometry as sgeo
#from scipy import interpolate

# from read_arason import *
#from plume_models import *
from optparse import OptionParser
from usgstable import USGStable


##Parse command line arguments and set defaults#####
parser = OptionParser()
parser.add_option(
    "-v",
    type="string",
    dest="volcano_name",
    default="kasatochi",
    help="name of volcano. Must match name in USGS table",
)
parser.add_option(
    "--d1",
    type="string",
    dest="d1str",
    default="",
    help="Start date in format YYYYMMDDHH",
)
parser.add_option(
    "--d2",
    type="string",
    dest="d2str",
    default="",
    help="Start date in format YYYYMMDDHH",
)
parser.add_option(
    "-t", type="float", dest="dt", default=3, help="Time increment in hours"
)

(options, args) = parser.parse_args()

volcano_name = options.volcano_name
dir = "../data/"
usgsname = "usgs_table.csv"


# usgs = usgstable(dir=dir, fname=usgsname)
# volc = usgs.volcano(volcano_name, accents=1)
# if volc == -1:
#   sys.exit(-1)
#
# print volc


##Input volcano name
## get volcano information from USGS table.
usgs = USGStable(dir=dir, fname=usgsname)
if len(volcano_name) == 1:
    volc = usgs.volclist(volcano_name)
volc = usgs.volcano(volcano_name, accents=1)
iloop = 0
while volc == -1 and iloop < 10:
    volcano_name = input(
        "That volcano was not found in the USGS database. Please give a different name. Enter 2 to exit."
    )
    print(volcano_name)
    if volcano_name == "2":
        sys.exit(-1)
    volc = usgs.volcano(volcano_name, accents=1)
    iloop += 1
print(volc)


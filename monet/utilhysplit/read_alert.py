#read_alert.py
#reads volcat alert xml file, creates hysplit traj. setup and control files
from monet.util import volcalert
from monet.utilhysplit import hcontrol
from monet.utilhysplit import forecast_data as fd
from datetime import timedelta as td

#Path and name of VOLCAT alert xml file
fpath = '/hysplit-users/allisonr/VOLCANO/'
fname = 'geocatASHOBJ.GOES-16.Full.xml'
#Met model for trajectory
met_model = 'GDAS0p5'
#Number of days for met fields needed
ndays = 2
#Project ID and working directory
pid = 'VALRT1'
wdir = '/hysplit-users/allisonr/VOLCANO/Trajectory/'

#Opens the xml file
alerts = volcalert.open_xml(fpath+fname)
#Returns the alert_header, alert_type, status, confidence, method, lat of radiative center
#lon of radiative center, and datetime object of eruption start
# An 8 dimensional array
data = volcalert.get_alert_values(alerts)
#Returns the nearest volcanoes to the alert as separate  arrays of:
#Volcano name, latitude, longitude, distance from alert, thermal signature, and volc. ID
#6 arrays, listed 
start_lat = data[5]
nearby = volcalert.get_nearby(alerts)
#Returns the name, lat, lon, therm signature, and volc. ID of closest volcano 
#to the detected radiative signature returned from get_alert_values
#Calculated based on minimum distance - requires all nearby values as input
closest = volcalert.get_closest(nearby[0],nearby[1],nearby[2],nearby[3],nearby[4],nearby[5])

#Find met fields necessary for the hysplit run
#findcycles_archive(datetime object start, datetime object end, met model,'Forward' or 'Back')
#findcycles_forecast(datetime object start, met model)
sdate = data[7]
edate = data[7] + td(days = ndays)
met = fd.findcycles_archive(sdate,edate,met_model,'Forward')

#Write SETUP and CONTROL files
#control = 



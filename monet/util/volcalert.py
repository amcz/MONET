#read_volcat_alert.py
#Reads XML alert files
#Extracts: date, time, lat, lon, volcano name
#max height (km), max height (ft), total mass (Tg), median effective radius (um)
import lxml.etree as let
import lxml.objectify as obj
import numpy as np
from datetime import datetime

#Opening xml file
def open_xml(xmlFile):
    print(xmlFile)
    with open(xmlFile) as f:
        xml = f.read()
        root = obj.fromstring(xml)    
    return root

#Extracting necessary values from alert
def get_alert_values(root):
    alert_header = root.alert.alert_header.attrib.get('value').strip()
    alert_type = root.alert.alert_type.attrib.get('value')
    status = root.alert.status.attrib.get('value')
    confidence = root.alert.confidence.attrib.get('value')
    method = root.alert.method.attrib.get('value')
    lat_rc = root.alert.lat_rc.attrib.get('value')
    lon_rc = root.alert.lon_rc.attrib.get('value')
    date_time = root.alert.object_date_time.attrib.get('value') #Put into datetime object
    date_time = datetime.strptime(date_time, '%Y-%m-%d %H:%M:%S').replace(second = 0)  #Strips seconds, leaving year, month, day, hour, minute as datetime object
    return alert_header, alert_type, status, confidence, method, lat_rc, lon_rc, date_time

#Extracting array of nearby possible volcanoes
def get_nearby(root):
    values=[]
    for data in root.alert.volcanoes.getchildren():
        for volc in data.getchildren():
            values.append( volc.attrib.get('value'))
    i=0
    name=[]
    lat=[]
    lon=[]
    dist=[]
    therm=[]
    vid=[]
    #.strip() removes the whitespace in the string
    while i < len(values):
        name.append(values[i].strip())
        lat.append(values[i+1].strip())
        lon.append(values[i+2].strip())
        dist.append(values[i+3].strip())
        therm.append(values[i+4].strip())
        vid.append(values[i+5].strip())
        i+=6
    return name, lat, lon, dist, therm, vid

#Finding closest volcano (minimum distance)
def get_closest(name, lat, lon, dist, therm, vid):
    closest = dist.index(min(dist))
    closest_name = name[closest]
    closest_lat = lat[closest]
    closest_lon = lon[closest]
    closest_therm = therm[closest]
    closest_vid = vid[closest] 
    return closest_name, closest_lat, closest_lon, closest_therm, closest_vid






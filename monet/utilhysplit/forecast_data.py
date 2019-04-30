import os

def findcycles_forecast(dstart, metdata):
    """dstart : datetime object. start date
       metdata : str 
       GFS, GFS0p5, GFS0p25, NAM 12 km, NAMAK, NAMHI 
    """
    FCTDIR='/pub/forecast/'
    cycles=["t00z","t06z","t12z","t18z"]
    ctimes=[0,6,12,18]
    metdirlist = []
    metfilelist=[]
    days = [""]
    if(metdata=="GFS"):
       meta = "gfsa"
       met = "gfsf"
    elif(metdata=="GFS0p5"):
       met = "gfs0p5f"
    elif(metdata=="GFS0p25"):
       met = "gfs0p25"
    elif(metdata=="NAM 12 km"):
       met = "namf"
    elif(metdata=="NAMAK"):
       met = "namsf.AK"
       meta = "namsa.AK"
    elif(metdata=="NAMHI"):
       met = "namsf.HI"
       meta = "namsa.HI"
    metdir1=FCTDIR + '/' + dstart.strftime("%Y%m%d") + '/'
    #print('<br>' + metdir1)
    cyclename = 'None'
    metnamefinal='No data found'
    for cyc, tms in zip(cycles, ctimes):
        metfilename = 'hysplit.' + cyc + '.' + met
        metname = metdir1 + metfilename
        if os.path.isfile(metname):     
            if tms < hour:
               metnamefinal = metfilename
               cyclename = str(tms) + ' UTC'
    for fs in days:
        metfilelist.append(metnamefinal + fs)         
        metdirlist.append(metdir1)         
    return metdirlist, metfilelist

def findcycles_archive(dstart, dend, metdata, direction):
     """dstart : datetime object start date
         dend : datetime object end date
         metdata : string : GDAS0p5, GDAS1
         direction: string : Forward, Back
    """
     from datetime import timedelta as td
     DIR='/pub/archives/'
     metdirlist = []
     metfilelist=[]
     datelist = []
     if (direction == 'Forward'):
         ndays = dend.day - dstart.day
         for x in range (0, ndays + 1):
             datelist.append(dstart + td(days = x))
     elif (direction == 'Back'):
         ndays = dstart.day - dend.day
         for  x in range (0, ndays + 1):
             datelist.append(dstart - td(days = x))
     if (metdata=="GDAS0p5"):
         met = "gdas0p5"
     elif (metdata=="GDAS1"):
         met = "gdas1"
     metdir1=DIR + met + '/'
     y = 0
     while y < len(datelist):
         metfilelist.append(datelist[y].strftime("%Y%m%d") + '_' + met)         
         metdirlist.append(metdir1)
         y += 1
     return metdirlist, metfilelist

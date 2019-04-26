
def findcycles(dstart, dend, metdata):
    """dstart : datetime object. start date
       dend : datetime, end date (not used)
       metdata : str 
       GFS, GFS0p5, GFS0p25, NAM 12 km, NAMAK, NAMHI 
    """
    FCTDIR='/pub/forecast/'
    cycles=["t00z","t06z","t12z","t18z"]
    ctimes=[0,6,12,18]
    metdirlist = []
    metfilelist=[]
    if(metdata=="GFS"):
       meta = "gfsa"
       met = "gfsf"
       days=[""]
    elif(metdata=="GFS0p5"):
       met = "gfs0p5f"
       days=[""]
    elif(metdata=="GFS0p25"):
       met = "gfs0p25"
       days=["000","024","048","072","096","120","144","168"]
    elif(metdata=="NAM 12 km"):
       met = "namf"
       days=[""]
    elif(metdata=="NAMAK"):
       met = "namsf.AK"
       meta = "namsa.AK"
       days=[""]
    elif(metdata=="NAMHI"):
       met = "namsf.HI"
       meta = "namsa.HI"
       days=[""]
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


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import numpy as np
import datetime
import time
import os
from os import path, chdir
from subprocess import call
import pandas as pd

"""
NAME: svhy.py
PRGMMR: Alice Crawford  ORG: ARL  
This code written at the NOAA  Air Resources Laboratory
ABSTRACT: choosing met files for HYSPLIT control file


"""



def getmetfiles(strfmt, sdate, runtime):
    mf = MetFiles(strfmt)
    return mf.get_files(sdate, runtime)



class MetFiles:

    def __init__(self, strfmt, hours=None, verbose=False):
        self.verbose = verbose
        self.strfmt = strfmt
        if not hours: 
            self.mdt = self.find_mdt()
        else:
            self.mdt = datetime.timedelta(hours=hours)

    def get_files(self, sdate, runtime):
        """
        sdate : datetime object
        runtime : integer. hours of runtime.
        """
        nlist = self.make_file_list(sdate, runtime)
        return self.process(nlist)

    def find_mdt(self):
        # finds time spacing between met files by
        # seeing which spacing produces a new file name.
        testdate = datetime.datetime(2010,10,12)
        mdtlist = [1,24,24*7,24*7*4, 24*356]
        
        file1 = testdate.strftime(self.strfmt)
        done = False
        iii=0
        while not done:
            dt = datetime.timedelta(hours=mdtlist[iii])
            d2 = testdate   + dt
            file2 = d2.strftime(self.strfmt)
            if d2 != testdate: done=True
            iii += 1
            if iii >= len(mdtlist): done=True
        return  dt

    def parse_week(self, edate):
        # used if week is in the strfmt (mostly for gdas1)
        temp = edate.strftime(self.strfmt)
        day = int(edate.strftime('%d'))
        if day < 7: 
           temp=temp.replace('week','w1') 
        elif day < 14: 
           temp=temp.replace('week','w2') 
        elif day < 21: 
           temp=temp.replace('week','w3') 
        elif day < 28: 
           temp=temp.replace('week','w4') 
        else:
           temp=temp.replace('week','w5') 
        return temp

    def make_file_list(self, sdate, runtime):
        nlist = []
        sdate = sdate.replace(tzinfo=None)
        if runtime < 0:
            runtime = abs(runtime)
            end_date = sdate
            sdate = end_date - datetime.timedelta(hours=runtime)
        else:
            end_date = sdate + datetime.timedelta(hours=runtime)
        edate = sdate
        done = False
        if self.verbose:
            print("GETMET", sdate, edate, end_date, runtime)
        zzz = 0
        while not done:
            if 'week' in self.strfmt:
               temp = self.parse_week(edate)
            else:
               temp = edate.strftime(self.strfmt)
            edate = edate + self.mdt
            if not path.isfile(temp):
              temp = temp.lower()
            if not path.isfile(temp):
                  print("WARNING", temp, " meteorological file does not exist")
            else:
              if temp not in nlist:
                 nlist.append(temp)
            if edate > end_date:
                done = True
            if zzz > 50: done=True 
            zzz+=1
        return nlist

    def process(self, nlist):
            # convert list of filenames with full path to 
            # list of directories and list of filenames
            # and then zips the lists to return list of tuples.
            mfiles = []
            mdirlist = []
            for temp in nlist:
              si = [x for x, char in enumerate(temp) if char=='/']
              si = si[-1]
              fname = temp[si+1:]
              mdir = temp[0:si+1] 
              mfiles.append(fname)
              mdirlist.append(mdir)
            return list(zip(mdirlist, mfiles))


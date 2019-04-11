#import volcat
import volcat_monet

flist=["geocatL2.GOES-16.Full_Disk.2019056.164530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.170030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.171530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.173030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.174530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.180030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.181530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.183030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.184530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.190030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.191530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.193030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.194530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.200030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.201530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.203030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.204530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.210030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.211530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.213030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.214530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.220030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.221530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.223030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.224530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.230030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.231530.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.233030.hdf",
       "geocatL2.GOES-16.Full_Disk.2019056.234530.hdf"]

flist=["geocatL2.GOES-16.Full_Disk.2019056.190030.hdf"]

direct="/pub/ECMWF/JPSS/volcat/reventador/"

for fname in flist:
    print(fname)
    #vdata=volcat.VolcatHDF(direct+fname, verbose=0)
    vdata=volcat_monet.VolcatHDF(direct+fname, verbose=0)
    #vdata.quickplot()
    vdata.geomap_plot()

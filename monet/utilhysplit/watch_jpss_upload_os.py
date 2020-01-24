#!/hysplit_users/allisonr/anaconda3/bin/python
#watch_jpss_upload_os.py
#Monitors ftp folder (/pub/jpsss_upload/) for new xml file
#If new file, then run file (run_hysp_alert.py)
#New file is found in directory: /hysplit-users/allisonr/Alice/MONET/monet/utilhysplit/
#Should create a hysplit control file, and run hysp trajectory whenever a new xml file 
#is detected in the ftp folder

import time
import os
import re
import json
import numpy as np
from monet.utilhysplit import run_hysptraj_alert as rha

#Path of ftp folder
path_to_watch = '/pub/jpsss_upload/'
#Creating list of current files in ftp folder
current = list(fi for fi in os.listdir(path_to_watch) if fi.endswith('.xml'))
#Reading the list of files from the last iteration
direct = '/hysplit-users/allisonr/VOLCANO/'
with open(direct+'original_list.txt', 'r') as fs:
    original = json.loads(fs.read())
#Determining what was added and what was removed
added = [fs for fs in current if not fs in original]
removed = [fs for fs in original if not fs in current]

original1 = original
#If files were added, loop for all files added 
if added:
    print('Added:', ',\n '.join(added),'\n')
    print('Adding '+str(len(added))+' files')
    f = 0
    while f < len(added):
        add = ''.join(added[f])   #Convert to string from list
        print(add)
        out = rha.make_files(path_to_watch,add)
        rha.run_hysp(out[0], out[1])
        rha.make_traj(out[0], out[1])
        rha.compress_kml(out[0], out[1])
        rha.make_gif(out[0], out[1])
        
        #Check if trajectory figure is generated
        trajfile = rha.check_for_figs(out[0],out[1])
        if trajfile == True:
            print('Trajectory file does exist!\n')
            original1.append(add)
        f += 1

if removed:
    print('Files removed:',', \n'.join(removed),'\n')
    g = 0
    while g < len(removed):
        remov = ''.join(removed[g])
        original1.remove(remov)
        g += 1

if added or removed:
    #Write current list of files that successfully created trajectory figs to original_list2.txt
    print('Updates recorded to file!\n')
    with open(direct+'/original_list2.txt', 'w') as fis:
                fis.write(json.dumps(original1))
    os.system('mv '+direct+'original_list2.txt '+direct+'original_list.txt')
else:
    print('No updates to '+path_to_watch+' folder!\n')







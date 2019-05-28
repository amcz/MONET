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
from monet.utilhysplit import run_hysp_alert as rha

#Path of ftp folder
path_to_watch = '/pub/jpsss_upload/'
#Creating list of current files in ftp folder
current = os.listdir(path_to_watch)
#Reading the list of files from the last iteration
with open('/hysplit-users/allisonr/VOLCANO/original_list.txt', 'r') as fs:
    original = json.loads(fs.read())
#Determining what was added and what was removed
added = [fs for fs in current if not fs in original]
removed = [fs for fs in original if not fs in current]

#If files were added, loop for all files added 
if added:
    print('Added:', ', '.join(added))
    f = 0
    while f < len(added):
        add = ''.join(added[f])   #Convert to string from list
        #See file run_hysp_alert.py for information about functions
        out = rha.make_files(path_to_watch,add)
        rha.run_hysp(out[0],out[1])
        rha.make_traj(out[0],out[1])
        f += 1

if removed:
    print('Files removed:',', '.join(removed))

if current == original:
    print('No updates to '+path_to_watch+' folder!')

#Write current list of files to original_list.txt - update for next iteration
with open('/hysplit-users/allisonr/VOLCANO/original_list.txt', 'w') as fis:
    fis.write(json.dumps(current))

#os.system('mv original_list2.txt original_list.txt')

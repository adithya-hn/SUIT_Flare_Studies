
import os
import timeit
import numpy as np
import shutil
import matplotlib.pyplot as plt
import drms
import sunpy.map
from sunpy.net import attrs as a  
from datetime import datetime, timedelta
import requests
import re
from sunpy.net import Fido
from astropy import units as u
import sunpy.net.jsoc as jsoc
from sunpy.net import attrs as a
from itertools import repeat

import requests


client = drms.Client()

save_dir = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June02_Flare/SHARP_data_2/'
email_id = 'adithya1@atulbhats.com'


result = Fido.search(a.Time("2024-06-02 02:00:00", "2024-06-02 06:00:00"),
                     a.jsoc.Series("hmi.sharp_cea_720s"),
                     a.jsoc.PrimeKey("HARPNUM", 11297),
                     a.jsoc.Notify(email_id))
print(result)
file = Fido.fetch(result,path=save_dir)


print(f"Task Completed")

#Remove From here if running over SSH or if Beep annoys you. ;)  

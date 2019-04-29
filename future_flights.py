from astropy.io import ascii
import numpy as np
import warnings
import pyb_traj
import get_gfs
import pyb_io
import time
import sys
import os

warnings.filterwarnings("ignore")

def future_flights(weather_file=None, datestr=None, utc_hour=None, lat0=None, lon0=None, alt0=None, descent_only=True, latest_data=False, next_point='0'):

	time0 = time.time()

	in_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

	lat0 = float(lat0)
	lon0 = float(lon0)
	alt0 = float(alt0)

	if descent_only == 'True':
		descent_only = True
	else:
		descent_only = False

	if next_point == '2':
		ext = 'next2_points/'
	elif next_point == '1':
		ext = 'next_points/'		
	else:
		ext = 'all_points/'

	if weather_file is None:
		hhhh, hhh = get_gfs.get_closest_gfs_file(datestr=datestr, utc_hour=utc_hour, verbose=True)
		file = 'gfs_4_' + datestr +'_' + str(hhhh) + '_' +str(hhh)
	else:
		hhhh, hhh = get_gfs.get_gfs_file(weather_file=weather_file, verbose=True)
		file, datestr = weather_file[:-5], weather_file[6:14]

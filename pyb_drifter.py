""" method to determine trajectories for one starting position with different times """

from astropy.io import ascii
import datetime as dt
import numpy as np
import time
import sys
import os

import info_searcher as inf
import param_file as p
import pyb_runner
import pyb_io

def get_ini(datestr=None, descent_only=False, next_point='1'):

	dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

	if descent_only:
		fname = 'descent_only_start' + next_point + '.txt'
	else:
		print('No file with specified flight data!')
		sys.exit()

	lines = [line.rstrip('\n').split(' ') for line in open(dir + fname)]
	for i in range(len(lines)):
		if lines[i][0] == datestr:
			break

	return lines[i][1], lines[i][2], lines[i][3], lines[i][4]

def drifter(run=None, datestr=None, utc_hour=None, lat0=None, lon0=None, alt0=None, params=None, balloon=None, verbose=False):

	time0 = time.time()

	out_dir = '/home/ellen/Desktop/SuperBIT/Output'

	if run == None:

		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

		files = [filename for filename in os.listdir(out_dir) if now_str in filename]
		run = now_str + '_' + str(len(files))

	dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

	if params == None:

		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		else:
			next_point = '0'
		interpolate = p.interpolate

	else:

		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		else:
			next_point = '0'
		interpolate = bool(params[2])

	if balloon == None:
		balloon = p.balloon

	if utc_hour == None or lat0 == None or lon0 == None or alt0 == None:
		if descent_only:
			utc_hour, lat0, lon0, alt0 = get_ini(datestr=datestr, descent_only=descent_only, next_point=next_point)
		else:
			utc_hour, lat0, lon0, alt0 = get_ini(datestr=datestr, descent_only=descent_only)

	if verbose:

		print('\nGeneral Parameters')
		print('----------')
		print('descent_only: ' + str(descent_only))
		if descent_only:
			print('starting point: ' + next_point)
		print('interpolate: ' + str(interpolate))
		print('----------')
		print('\nBalloon/Parachute Parameters')
		print('----------')
		print('altitude step: ' + str(balloon['altitude_step']) + ' m')
		print('equipment mass: ' + str(balloon['equip_mass']) + ' kg')
		print('parachute Cd: ' + str(balloon['Cd_parachute']))
		print('parachute radius: ' + str(round(np.sqrt(balloon['parachute_areas'][0]/np.pi), 2)) + ' m')
		print('parachute area: ' + str(round(balloon['parachute_areas'][0], 2)) + ' m^2')
		print('----------')
		print('\nRunning date: ' + datestr)
		print('Starting point: ' + str(lat0) + ' lat., ' + str(lon0) + ' lon., ' + str(alt0) + ' m\n')
		print('----------\n')

	drift_times = np.arange(0., 60., 30.)

	for drift_time in drift_times:

		print('Running drift time: ' + str(drift_time) + ' minutes')

		pyb_runner.runner(datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, params=[descent_only, next_point, interpolate, drift_time], balloon=balloon)

		print('\n----------\n')

	folders = inf.seach_info_drift_times(drift_times=drift_times)
	pyb_io.merge_kml(datestr=datestr, folders=folders, params=params, balloon=balloon, drift_times=drift_times)

	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	drifter(datestr='20181215', verbose=True)
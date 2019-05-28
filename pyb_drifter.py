""" method to determine trajectories for one starting position with different times """

from astropy.io import ascii
import numpy as np
import time
import sys

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

def drifter(datestr=None, utc_hour=None, lat0=None, lon0=None, alt0=None, params=None, balloon=None, verbose=False):

	time0 = time.time()

	dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

	if params == None:

		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		interpolate = p.interpolate

	else:

		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		interpolate = bool(params[-2])

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
		print('parachute area: ' + str(round(balloon['parachute_areas'][0], 2)) + ' m^2')
		print('----------')
		print('\nRunning date: ' + datestr)
		print('Starting point: ' + str(lat0) + ' lat., ' + str(lon0) + ' lon., ' + str(alt0) + ' m\n')
		print('----------\n')

	drift_times = np.arange(0., 120., 30.)

	for drift_time in drift_times:

		print('Running drift time: ' + str(drift_time) + ' minutes')

		if descent_only:
			pyb_runner.run(datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, params=[descent_only, next_point, interpolate, drift_time], balloon=balloon)
		else:
			pyb_runner.run(datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, params=[descent_only, interpolate, drift_time], balloon=balloon)

		print('\n----------\n')

	pyb_io.merge_kml(datestr=datestr, params=params, drift_times=drift_times)

	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	drifter(datestr='20181215', verbose=True)
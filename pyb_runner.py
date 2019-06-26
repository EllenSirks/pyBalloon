""" 
Script for running a single trajectory.
Certain parameters are read from the param_file, or can be manually fed in.
"""

from astropy.io import ascii
import datetime as dt
import numpy as np
import warnings
import time
import sys
import os

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl

import pyb_traj
import get_gfs
import pyb_io

import test

import param_file as p

warnings.filterwarnings("ignore")

# method to run entire simulation of flight
# requires the starting location & starting point (0 for highest), whether or not its descent only, and
# either the date & time of the starting point OR the relevent weather file

def runner(datestr=None, utc_hour=None, lat0=None, lon0=None, alt0=None, balloon=None, params=None, run=None, make_folder=True, verbose=False):

	print('')
	time0 = time.time()

	# starting location
	lat0, lon0, alt0 = float(lat0), float(lon0), float(alt0)
	loc0 = (lat0, lon0, alt0)

	# set parameters
	if params == None:
		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		else:
			next_point = '0'
		interpolate = p.interpolate
		drift_time = float(p.drift_time)
		resolution = float(p.resolution)
	else:
		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		else:
			next_point = '0'
		interpolate = bool(params[2])
		drift_time = float(params[3])
		resolution = float(params[4])

	if balloon == None:
		balloon = p.balloon

	if verbose:
		print('General Parameters')
		print('----------')
		print('descent_only: ' + str(descent_only))
		if descent_only:
			print('starting point: ' + next_point)
		print('interpolate: ' + str(interpolate))
		print('drift time: ' + str(drift_time) + ' minutes')
		print('resolution of forecasts: ' + str(resolution) + ' degrees')
		print('----------')
		print('\nBalloon/Parachute Parameters')
		print('----------')
		print('altitude step: ' + str(balloon['altitude_step']) + ' m')
		print('equipment mass: ' + str(balloon['equip_mass']) + ' kg')
		print('parachute Cd: ' + str(balloon['Cd_parachute']))
		print('parachute radius: ' + str(round(np.sqrt(balloon['parachute_areas'][0]/np.pi), 2)) + ' m')
		print('parachute area: ' + str(round(balloon['parachute_areas'][0], 2)) + ' m^2')
		print('----------\n')

		print('Running date/time: ' + datestr + ', ' + utc_hour + ' hr')
		print('Starting point: ' + str(lat0) + ' lat., ' + str(lon0) + ' lon., ' + str(alt0) + ' m\n')

	# create relevant paths/folders
	dir_base = '/home/ellen/Desktop/SuperBIT/'

	in_dir = dir_base + '/Weather_data/GFS/'
	traj_dir = dir_base + 'Output/'

	if run != None:
		folder = run + '/'
	else:
		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)
		files = [filename for filename in os.listdir(traj_dir) if now_str in filename]
		run = now_str + '_' + str(len(files))
		folder = run  + '/'		

	kml_dir = traj_dir + folder + 'kml_files/'
	params_dir = traj_dir + folder
	traj_dir += folder + 'trajectories/'

	# check if the paths exist/make them
	if not os.path.exists(kml_dir):
		os.makedirs(kml_dir)
	if not os.path.exists(traj_dir):
		os.makedirs(traj_dir)

	# check if want interpolation
	if interpolate:
		files = test.get_interpolation_gfs_files(datestr=datestr, utc_hour=utc_hour, resolution=resolution)
	else:
		files = test.get_closest_gfs_file(datestr=datestr, utc_hour=utc_hour, resolution=resolution)

	# calculate the trajectory of the balloon
	trajectories = pyb_traj.run_traj(weather_files=files, loc0=loc0, datestr=datestr, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only, interpolate=interpolate, drift_time=drift_time, resolution=resolution)

	# write out trajectory to file
	traj_file = traj_dir + 'trajectory_' + files[0][6:14] + '_' + str(utc_hour) + '_' + str(loc0) + '.dat'
	ascii.write([trajectories['lats'], trajectories['lons'], trajectories['alts'], trajectories['dists'], trajectories['times'], trajectories['speeds'], trajectories['z_speeds'], trajectories['omegas'], trajectories['temperatures']], \
	 traj_file, names=['lats', 'lons', 'alts', 'dists', 'times', 'speeds', 'z_speeds', 'omegas', 'temps'], overwrite=True)

	# write parameters of this run to file
	f = open(params_dir + 'params.txt', 'w+')
	f.write('General parameters\n')
	f.write('----------------------\n')
	f.write('descent_only: ' + str(descent_only) + '\n')
	if descent_only:
		f.write('starting point: ' + str(next_point) + '\n')
	f.write('interpolate: ' + str(interpolate) + '\n')
	f.write('drift time: ' + str(drift_time) + ' min\n')
	f.write('resolution of forecasts: ' + str(resolution) + ' degrees\n')
	f.write('----------------------\n')
	f.write('\n')
	f.write('Balloon/Parachute parameters\n')
	f.write('----------------------\n')
	f.write('altitude step: ' + str(balloon['altitude_step']) + ' m\n')
	f.write('equipment mass: ' + str(balloon['equip_mass']) + ' kg\n')
	f.write('parachute Cd: ' + str(balloon['Cd_parachute']) + '\n')
	f.write('parachute radius: ' + str(round(np.sqrt(balloon['parachute_areas'][0]/np.pi), 2)) + ' m\n')
	f.write('parachute area: ' + str(round(balloon['parachute_areas'][0], 2)) + ' m^2\n')
	f.close()

	# add info to run information file
	lines = [line.rstrip('\n').split(' ') for line in open(dir_base + 'Output/runs_info.txt')]
	runs = [lines[i][0] for i in range(len(lines)) if i != 0]

	f = open(dir_base + 'Output/runs_info.txt', 'a+')
	if run not in runs:
		f.write('\n' + str(folder[:-1]) + ' ' + str(descent_only) + ' ' + str(next_point) + ' ' +  str(interpolate) + ' ' + str(drift_time) + ' ' + str(resolution) + ' ' + str(balloon['Cd_parachute']) + ' ' + str(np.sqrt(balloon['parachute_areas'][0]/np.pi)) + ' ' + \
		 str(balloon['altitude_step'])  + ' ' + str(balloon['equip_mass']) + ' ' + str(balloon['balloon_mass']) + ' ' +  str(balloon['fill_radius']) + ' ' + str(balloon['radius_empty']) + ' ' + \
		 str(balloon['burst_radius']) + ' ' + str(balloon['thickness_empty']) + ' ' + str(balloon['Cd_balloon']) + ' ' + str(balloon['simple_ascent_rate']) + ' ' + str(balloon['parachute_change_altitude']))

	f.close()

	# highest point in main-run trajectory (bit obsolete for descent only)
	idx, = np.where(trajectories['alts'] == np.max(trajectories['alts']))
	latx = trajectories['lats'][idx][0]
	lonx = trajectories['lons'][idx][0]
	altx = trajectories['alts'][idx][0]
	timex = trajectories['times'][idx][0]

	print('Starting location: (' + str(float(latx)) + ', ' + str(float(lonx)) + '), altitude ' +  str(float(altx)) + ', at %.0f minutes' % (timex))

	other_info = [(latx, lonx, altx, 'Burst point', '%.0f minutes, %.0f meters' % (timex, altx))]

	# write out file for google-earth
	kml_fname = kml_dir + files[0][6:14] + '_' + str(utc_hour) + '_' + str(loc0) + '.kml'
	pyb_io.save_kml(kml_fname, trajectories, other_info=other_info, params=params)

	print('\nProgram finished in %.1f s' % (time.time() - time0))

if __name__ == "__main__":

	runner(datestr=sys.argv[1], utc_hour=sys.argv[2], lat0=sys.argv[3], lon0=sys.argv[4], alt0=sys.argv[5], verbose=True)

""" 
Script for running a single trajectory.
Certain parameters are read from the param_file, or can be manually fed in.
"""

from astropy.io import ascii
import numpy as np
import warnings
import time
import sys
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import pyb_traj
import get_gfs
import pyb_io

import param_file as p

warnings.filterwarnings("ignore")

# method to run entire simulation of flight
# requires the starting location & starting point (0 for highest), whether or not its descent only, and
# either the date & time of the starting point OR the relevent weather file

def run(datestr=None, utc_hour=None, lat0=None, lon0=None, alt0=None, balloon=None, params=None, verbose=False):

	print('')
	time0 = time.time()

	# starting location
	lat0, lon0, alt0 = float(lat0), float(lon0), float(alt0)
	loc0 = (lat0, lon0, alt0)

	# set parameters
	if params == None:
		descent_only = p.descent_only
		interpolate = p.interpolate
		if descent_only:
			next_point = p.next_point
		drift_time = float(p.drift_time)
	else:
		descent_only = bool(params[0])
		interpolate = bool(params[-2])
		if descent_only:
			next_point = str(params[1])
		drift_time = float(params[-1])

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
		print('----------')
		print('\nBalloon/Parachute Parameters')
		print('----------')
		print('altitude step: ' + str(balloon['altitude_step']) + ' m')
		print('equipment mass: ' + str(balloon['equip_mass']) + ' kg')
		print('parachute Cd: ' + str(balloon['Cd_parachute']))
		print('parachute area: ' + str(round(balloon['parachute_areas'][0], 2)) + ' m^2')
		print('----------\n')

		print('Running date/time: ' + datestr + ', ' + utc_hour + ' hr')
		print('Starting point: ' + str(lat0) + ' lat., ' + str(lon0) + ' lon., ' + str(alt0) + ' m\n')

	# create relevant paths/folders
	dir_base = '/home/ellen/Desktop/SuperBIT/Weather_data/'

	if descent_only:
		ext = 'descent_only/start_point' + next_point + '/'
	else:
		ext = 'ascent+descent/'

	in_dir = dir_base + 'GFS/'
	kml_dir = dir_base + 'kml_files/' + ext
	traj_dir = dir_base + 'Trajectories/' + ext
	end_dir = dir_base + 'Endpoints/' + ext

	# check if the paths exist/make them
	if not os.path.exists(kml_dir):
		os.makedirs(kml_dir)
	if not os.path.exists(traj_dir):
		os.makedirs(traj_dir)
	if not os.path.exists(end_dir):
		os.makedirs(end_dir)

	# check if want interpolation
	if interpolate:
		interp = '_interpolated_'
		# get all the files needed for the interpolation
		files = get_gfs.get_interpolation_gfs_files(datestr=datestr, utc_hour=utc_hour)
		file = files[0]

	else:
		interp = '_'
		# retrieve relevant (initial) weather_file
		files = get_gfs.get_closest_gfs_file(datestr=datestr, utc_hour=utc_hour)
		file = files[0]

	# calculate the trajectory of the balloon
	trajectories = pyb_traj.run_traj(weather_files=files, loc0=loc0, datestr=datestr, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only, interpolate=interpolate, drift_time=drift_time)

	# write out trajectory to file
	traj_file = traj_dir + 'trajectory_' + file[6:14] + '_' + str(utc_hour) + '_' + str(loc0) + interp + '+' + str(int(drift_time)).zfill(4) + 'min.dat'
	ascii.write([trajectories['lats'], trajectories['lons'], trajectories['alts'], trajectories['dists'], trajectories['times'], trajectories['speeds'], trajectories['temperatures']], traj_file, names=['lats', 'lons', 'alts', 'dists', 'times', 'speeds', 'temps'], overwrite=True)

	# highest point in main-run trajectory (bit obsolete for descent only)
	idx, = np.where(trajectories['alts'] == np.max(trajectories['alts']))
	latx = trajectories['lats'][idx][0]
	lonx = trajectories['lons'][idx][0]
	altx = trajectories['alts'][idx][0]
	timex = trajectories['times'][idx][0]

	print('Starting location: (' + str(float(latx)) + ', ' + str(float(lonx)) + '), altitude ' +  str(float(altx)) + ', at %.0f minutes' % (timex))

	other_info = [(latx, lonx, altx, 'Burst point', '%.0f minutes, %.0f meters' % (timex, altx))]

	# write out endpoint to file
	f = open(end_dir + 'endpoint_' + file[6:14] + '_' + str(utc_hour) + '_' + str(loc0) + interp + '+' + str(int(drift_time)).zfill(4) + 'min.dat','w+')
	f.write('lat lon alt\n')
	f.write(str(lat0) + ' ' + str(lon0) + ' ' + str(alt0) + '\n')
	f.write(str(trajectories['lats'][-1]) + ' ' + str(trajectories['lons'][-1]) + ' ' + str(trajectories['alts'][-1]))
	f.close()

	# write out file for google-earth
	# radius = np.sqrt((np.sqrt(np.sum(trajectories['sigmas_u']**2)))**2 + (np.sqrt(np.sum(trajectories['sigmas_v']**2)))**2)*1000.
	# no_steps = len(trajectories['times'])
	# radius = radius*200/no_steps
	radius=5000

	kml_fname = kml_dir + file[6:14] + '_' + str(utc_hour) + '_' + str(loc0) + interp + '+' + str(int(drift_time)).zfill(4) + 'min.kml'
	pyb_io.save_kml(kml_fname, trajectories, other_info=other_info, params=params, radius=radius)

	print('\nProgram finished in %.1f s' % (time.time() - time0))

	# return trajectories['distance']

if __name__ == "__main__":

	run(datestr=sys.argv[1], utc_hour=sys.argv[2], lat0=sys.argv[3], lon0=sys.argv[4], alt0=sys.argv[5], verbose=True)

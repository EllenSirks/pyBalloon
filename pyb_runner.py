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

import pyb_traj
import get_gfs
import pyb_io

import param_file as p

warnings.filterwarnings("ignore")

# method to run entire simulation of flight
# requires the starting location & starting point (0 for highest), whether or not its descent only, and
# either the date & time of the starting point OR the relevent weather file

def run(datestr=None, utc_hour=None, lat0=None, lon0=None, alt0=None, params=None):

	print('')
	time0 = time.time()

	if params == None:
		descent_only = p.descent_only
		interpolate = p.interpolate
		if descent_only:
			next_point = p.next_point
		drift_time = p.drift_time
	else:
		descent_only = bool(params[0])
		interpolate = bool(params[-2])
		if descent_only:
			next_point = str(params[1])
		drift_time = float(params[-1])

	# create relevant paths/folders
	dir_base = '/home/ellen/Desktop/SuperBIT/Weather_data/'

	if descent_only:
		ext = 'descent_only/start_point' + next_point + '/'
	else:
		ext = 'ascent+descent/'

	in_dir = dir_base + 'grb_files/'
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

	# check if want to interpolate
	if interpolate:
		interp = '_interpolated_'
	else:
		interp = '_'

	# starting location
	lat0, lon0, alt0 = float(lat0), float(lon0), float(alt0)
	loc0 = (lat0, lon0, alt0)

	# drift time
	drift_time = float(drift_time)

	# set balloon parameters
	balloon = {}
	balloon['altitude_step'] = 100.0 # meters
	balloon['equip_mass'] = 1.608 # kg
	balloon['balloon_mass'] = 1.50 # kg
	balloon['fill_radius'] = 2.122/2 # meters
	balloon['radius_empty'] = 2.41/2 # meters (flaccid body length - neck length)
	balloon['burst_radius'] = 9.44/2 # meters
	balloon['thickness_empty'] = 0.2 * 10**-3 # mm -> meters
	balloon['Cd_balloon'] = 0.5
	balloon['simple_ascent_rate'] = 5.0 # m/s

	# parachute parameters
	balloon['Cd_parachute'] = 0.5
	radius = 1.0
	balloon['parachute_areas'] = np.pi * np.array([radius])**2 # m^2
	balloon['parachute_change_altitude'] = None # meters

	# check if want interpolation
	if interpolate:
		# get all the files needed for the interpolation
		files = get_gfs.get_interpolation_gfs_files(datestr=datestr, utc_hour=utc_hour)
		file = files[0]

	else:
		# retrieve relevant weather_file
		files = get_gfs.get_closest_gfs_file(datestr=datestr, utc_hour=utc_hour, verbose=True)
		file = files[0]

	# calculate the trajectory of the balloon
	trajectories = pyb_traj.run_traj(weather_files=files, loc0=loc0, datestr=datestr, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only, interpolate=interpolate, drift_time=drift_time)

	traj_file = traj_dir + 'trajectory_' + file[6:14] + '_' + str(utc_hour) + '_' + str(loc0) + interp + '+' + str(int(drift_time)).zfill(4) + 'min.dat'
	ascii.write([trajectories['lats'], trajectories['lons'], trajectories['alts'], trajectories['dists'], trajectories['times'], trajectories['speeds']], traj_file, names=['lats', 'lons', 'alts', 'dists', 'times', 'speeds'], overwrite=True)

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
	kml_fname = kml_dir + file[6:14] + '_' + str(utc_hour) + '_' + str(loc0) + interp + '+' + str(int(drift_time)).zfill(4) + 'min.kml'
	pyb_io.save_kml(kml_fname, trajectories, other_info=other_info, params=params)

	print('\nProgram finished in %.1f s' % (time.time() - time0))

if __name__ == "__main__":

	run(datestr=sys.argv[1], utc_hour=sys.argv[2], lat0=sys.argv[3], lon0=sys.argv[4], alt0=sys.argv[5])

"""Script for running a single trajectory. Certain can be manually fed in or read from the param_file."""

from astropy.io import ascii
import datetime as dt
import numpy as np
import warnings
import sys, os
import time

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl

import pyb_traj
import get_gfs
import pyb_io

import param_file as p

warnings.filterwarnings("ignore")

# method to run entire simulation of flight
# requires the starting location & starting point (0 for highest), whether or not its descent only, and the date & time of the starting point
# if parameters are not supplied, it takes it from the param file

def runner(datestr=None, utc_hour=None, lat0=None, lon0=None, alt0=None, balloon=None, params=None, run=None, print_verbose=False, write_verbose=False, add_run_info=True):

	print('')
	time0 = time.time()
	now = dt.datetime.now()

	############################################################################################################ <---- set trajectory parameters

	# starting location
	lat0, lon0, alt0 = float(lat0), float(lon0), float(alt0)
	loc0 = (lat0, lon0, alt0)

	# general parameters
	next_point = '0'
	if params == None:
		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		interpolate = p.interpolate
		drift_time = float(p.drift_time)
		resolution = float(p.resolution)
		vz_correct = bool(p.vz_correct)
		hr_diff = int(p.hr_diff)
		
	else:
		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		interpolate = bool(params[2])
		drift_time = float(params[3])
		resolution = float(params[4])
		vz_correct = bool(params[5])
		hr_diff = int(params[6])

	# balloon parameters
	if balloon == None:
		balloon = p.balloon

	# print out parameters to terminal
	if print_verbose:
		pyb_io.print_verbose(datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, descent_only=descent_only, next_point=next_point, interpolate=interpolate,\
		 drift_time=drift_time, resolution=resolution, vz_correct=vz_correct, hr_diff=hr_diff, balloon=balloon)

	############################################################################################################ <---- set/create paths

	# initialise paths
	dir_base = '/home/ellen/Desktop/SuperBIT/'
	in_dir = dir_base + '/Weather_data/GFS/'
	traj_dir = dir_base + 'Output/'
	fig_dir = dir_base + 'Output/'

	# determine run number (datestr + no.)
	if run == None:
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)
		files = [filename for filename in os.listdir(traj_dir) if now_str in filename]
		run = now_str + '_' + str(len(files))

	# create path names with correct run
	kml_dir = traj_dir + run + '/Kml_files/'
	params_dir = traj_dir + run + '/'
	traj_dir += run + '/Trajectories/'
	fig_dir += run + '/Figs/Checks/'

	# check if the paths exist/make them
	if not os.path.exists(kml_dir):
		os.makedirs(kml_dir)
	if not os.path.exists(traj_dir):
		os.makedirs(traj_dir)
	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	############################################################################################################ <---- calculation trajectories

	# get weather files depending on if we want interpolation
	if interpolate:
		files = get_gfs.get_interpolation_gfs_files(datestr=datestr, utc_hour=utc_hour, resolution=resolution)
	else:
		files = get_gfs.get_closest_gfs_file(datestr=datestr, utc_hour=utc_hour, resolution=resolution, hr_diff=hr_diff)

	# calculate the trajectory of the balloon
	trajectories, fig_dicts = pyb_traj.run_traj(weather_files=files, loc0=loc0, datestr=datestr, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only, interpolate=interpolate, \
		drift_time=drift_time, resolution=resolution, vz_correct=vz_correct, hr_diff=hr_diff)

	############################################################################################################ <---- create/write output 

	# initialise trajectory and kml file names
	traj_file = traj_dir + 'trajectory_' + datestr + '_' + str(utc_hour) + '_' + str(loc0)
	kml_fname = kml_dir + datestr + '_' + str(utc_hour) + '_' + str(loc0)

	# check for double file names and add count if there are (usefull for testing drifting/forecasts)
	count = 2
	while True:
		if not os.path.isfile(traj_file + '.dat'):
			traj_file += '.dat'
			kml_fname += '.kml'
			break
		else:
			if count == 2:
				traj_file += '_' + str(count)
				kml_fname += '_' + str(count)
			else:
				ind1 = [pos for pos, char in enumerate(traj_file) if char == '_'][-1]
				ind2 = [pos for pos, char in enumerate(kml_fname) if char == '_'][-1]
				traj_file = traj_file[:ind1+1] + str(count)
				kml_fname = kml_fname[:ind2+1] + str(count)
			count += 1

	# write out trajectory file
	ascii.write([trajectories['lats'], trajectories['lons'], trajectories['alts'], trajectories['dists'], trajectories['times'], trajectories['speeds'], trajectories['z_speeds'], \
		trajectories['omegas'], trajectories['temperatures']], traj_file, names=['lats', 'lons', 'alts', 'dists', 'times', 'speeds', 'z_speeds', 'omegas', 'temps'], overwrite=True)

	# write parameters of this run to file
	if write_verbose:
		pyb_io.write_verbose(params_dir=params_dir, datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, descent_only=descent_only, next_point=next_point, \
			interpolate=interpolate, drift_time=drift_time, resolution=resolution, vz_correct=vz_correct, hr_diff=hr_diff, balloon=balloon)

	# add run info to run_info.txt file
	pyb_io.write_run_info(add_run_info=add_run_info, run=run, descent_only=descent_only, next_point=next_point, interpolate=interpolate, drift_time=drift_time, resolution=resolution, \
		vz_correct=vz_correct, hr_diff=hr_diff, balloon=balloon)

	# highest point in main-run trajectory (bit redundant for descent only)
	idx, = np.where(trajectories['alts'] == np.max(trajectories['alts']))
	latx = trajectories['lats'][idx][0]
	lonx = trajectories['lons'][idx][0]
	altx = trajectories['alts'][idx][0]
	timex = trajectories['times'][idx][0]

	# write out file for google-earth
	other_info = [(latx, lonx, altx, 'Burst point', '%.0f minutes, %.0f meters' % (timex, altx))]
	pyb_io.save_kml(kml_fname, trajectories, other_info=other_info, params=params)

	# save figs with interpolation checks
	for i in range(len(fig_dicts)):
		for key in fig_dicts[i].keys():
			fig_dicts[i][key].savefig(fig_dir + datestr + '_' + key + '_check' + str(i) + '.png')

	############################################################################################################

	print('\nProgram finished in %.1f s' % (time.time() - time0))

if __name__ == "__main__":

	runner(datestr=sys.argv[1], utc_hour=sys.argv[2], lat0=sys.argv[3], lon0=sys.argv[4], alt0=sys.argv[5], verbose=True)

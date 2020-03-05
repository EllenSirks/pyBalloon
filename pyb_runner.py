"""
Script for running a single trajectory. Certain can be manually fed in or read from the param_file.
"""

from astropy.io import ascii
import datetime as dt
import numpy as np
import sys, os
import math
import time

import pyb_traj
import get_gfs
import pyb_aux
import pyb_io

import param_file as p

#################################################################################################################

def runner(ini_conditions=None, balloon=None, params=None, run=None, overwrite=False):
	"""
	Method to run entire simulation of flight and write out calculated trajectories and predicted endpoint

	Arguments
	=========
	ini_conditions : tuple of strings
		(Date of initial point, Initial time of trajectory, (latitude in degrees, longitude in degrees, altitude in km) of initial point
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	run : string
		String indicating which run folder the results are to be stored in
	overwrite : bool
		If True, overwrite trajectory file with same ini_conditions in folder
	"""

	############################################################################################################ <---- set trajectory parameters

	# starting location
	datestr, utc_hour, loc0 = ini_conditions
	lat0, lon0, alt0 = loc0

	# general parameters, if balloon/parachute and other parameters are not supplied, they are read from param_file.py
	descent_only, drift_time, resolution, hr_diff, check_elevation, live, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	# change altitude if it is underground (not always accurate)
	if not descent_only:
		elevation = pyb_aux.get_elevation(lon=lon0, lat=lat0)
		if elevation > alt0:
			alt0 = elevation
			loc0 = lat0, lon0, alt0

	# print out parameters to terminal
	pyb_io.print_verbose(datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon)

	############################################################################################################ <---- set/create paths

	# create paths
	base_dir, kml_dir, traj_dir, fig_dir = pyb_io.get_and_make_paths(run=run)

	############################################################################################################ <---- calculation trajectories

	# get weather files
	files = get_gfs.get_interpolation_gfs_files(datestr=datestr, utc_hour=utc_hour, resolution=resolution, hr_diff=hr_diff, live=live)

	# calculate the trajectory of the balloon
	trajectories, used_weather_files = pyb_traj.run_traj(weather_files=files, ini_conditions=ini_conditions, params=params, balloon=balloon)

	############################################################################################################ <---- create/write output 

	# write out trajectories file
	pyb_io.create_trajectory_files(traj_dir=traj_dir, data=trajectories, ini_conditions=ini_conditions, overwrite=overwrite)

	# determine error
	errs = pyb_io.determine_error(data=trajectories)

	# write out kml file for google earth
	pyb_io.save_kml(kml_dir=kml_dir, data=trajectories, ini_conditions=ini_conditions, errs=errs, overwrite=overwrite)

	# save descent rate figure
	pyb_io.make_descent_rate_plot(fig_dir=fig_dir, data=trajectories, ini_conditions=ini_conditions)

	############################################################################################################

	# write parameters of this run to file
	pyb_io.write_verbose(params_dir=base_dir, params=params, balloon=balloon)

	# add run info to run_info.txt file
	pyb_io.write_run_info(run=run, params=params, balloon=balloon)

	# save weather files used
	pyb_io.save_used_weather_files(utc_hour=utc_hour, save_dir=base_dir, used_weather_files=used_weather_files, trajectories=trajectories)

	return trajectories

#################################################################################################################

if __name__ == "__main__":

	time0 = time.time()

	print('')

	############################################################################################################

	datestr = sys.argv[1]
	utc_hour = float(sys.argv[2])
	loc0 = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
	ini_conditions = (datestr, utc_hour, loc0)

	params = None

	if len(sys.argv) > 6:
		run = sys.argv[6]
	else:
		run = None

	trajectories = runner(ini_conditions=ini_conditions, params=params, run=run)

	############################################################################################################

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write(('Program finished in ' + str(round((time.time() - time0), 1)) + ' s').ljust(60) + '\n')

#################################################################################################################
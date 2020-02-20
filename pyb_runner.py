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

def runner(datestr=None, utc_hour=None, loc0=None, balloon=None, params=None, run=None, print_verbose=False, write_verbose=False, add_run_info=True, output_figs=False, overwrite=False):
	"""
	Method to run entire simulation of flight and write out calculated trajectories and predicted endpoint

	Arguments
	=========
	datestr : string
		Date of initial point
	utc_hour : float
		Time of initial point
	loc0 : floats in tuple
		(latitude in degrees, longitude in degrees, altitude in km) of initial point
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	run : string
		String indicating which run folder the results are to be stored in
	print_verbose : bool
		If True, the parameters used will be printed to the command line
	write_verbose : bool
		If True, the parameters used are written to a file
	add_run_info : bool 
		If True, the parameters used are appended to the run_info.txt file in the Output folder
	output_figs : bool
		If True, figures showing the grib data before and after interpolation between altitude steps are created and saved
	overwrite : bool
		If True, overwrite trajectory file with same name in folder
	check_elevation : bool
		If True, check elevation to make sure the code does not predict the trajectory to go underground.
	"""

	print('')
	now = dt.datetime.now()

	############################################################################################################ <---- set trajectory parameters

	# starting location
	lat0, lon0, alt0 = loc0

	# general parameters, if balloon/parachute and other parameters are not supplied, they are read from param_file.py
	descent_only, next_point, time_interpolate, grid_interpolate, drift_time, resolution, hr_diff, check_elevation, live, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	# change altitude if it is underground (not always accurate)
	if not descent_only:
		elevation = pyb_aux.get_elevation(lon=lon0, lat=lat0)
		if elevation > alt0:
			alt0 = elevation
			loc0 = lat0, lon0, alt0

	# print out parameters to terminal
	if print_verbose:
		pyb_io.print_verbose(datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon)

	############################################################################################################ <---- set/create paths

	# determine run number (datestr + no.)
	if run == None:
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)
		files = [filename for filename in os.listdir(p.path + p.output_folder) if now_str in filename]
		run = now_str + '_' + str(len(files))

	# create path names
	base_dir = p.path + p.output_folder + run + '/'
	kml_dir = base_dir + p.kml_folder
	traj_dir = base_dir + p.traj_folder
	fig_dir = base_dir + p.fig_folder + p.check_figs_folder
	fig_dir_checks = fig_dir + 'InterpolationChecks/'

	# check if the paths exist/make them
	if not os.path.exists(kml_dir):
		os.makedirs(kml_dir)
	if not os.path.exists(traj_dir):
		os.makedirs(traj_dir)
	if not os.path.exists(fig_dir_checks):
		os.makedirs(fig_dir_checks)

	############################################################################################################ <---- calculation trajectories

	# get weather files depending on if we want interpolation
	if time_interpolate:
		files = get_gfs.get_interpolation_gfs_files(datestr=datestr, utc_hour=utc_hour, resolution=resolution, hr_diff=hr_diff, live=live)
	else:
		files = get_gfs.get_closest_gfs_file(datestr=datestr, utc_hour=utc_hour, resolution=resolution, hr_diff=hr_diff)

	# calculate the trajectory of the balloon
	trajectories, fig_dicts, used_weather_files, time_diffs = pyb_traj.run_traj(weather_files=files, datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, \
		balloon=balloon, output_figs=output_figs)

	############################################################################################################ <---- create/write output 

	# write out trajectories and create kml file name based on trajectory file
	kml_fname = pyb_io.create_trajectory_files(traj_dir=traj_dir, kml_dir=kml_dir, datestr=datestr, utc_hour=utc_hour, loc0=loc0, trajectories=trajectories, params=params, overwrite=overwrite)

	# write parameters of this run to file
	if write_verbose:
		pyb_io.write_verbose(params_dir=base_dir, params=params, balloon=balloon)

	# add run info to run_info.txt file
	pyb_io.write_run_info(add_run_info=add_run_info, run=run, params=params, balloon=balloon)

	# save weather files used
	pyb_io.save_used_weather_files(utc_hour=utc_hour, save_dir=base_dir, used_weather_files=used_weather_files, trajectories=trajectories)

	# save figs with interpolation checks
	if output_figs:
		for i in range(len(fig_dicts)-1):
			for key in fig_dicts[i].keys():
				fig_dicts[i][key].savefig(fig_dir_checks + datestr + '_' + key + '_check' + str(i+1) + '.png')

	# save descent rate figure
	pyb_io.make_descent_rate_plot(directory=fig_dir, data=trajectories, datestr=datestr, utc_hour=utc_hour, loc0=loc0)

	# determine error
	parallel_err, perp_err = pyb_io.determine_error(trajectories=trajectories, params=params)

	# highest point in main-run trajectory (bit redundant for descent only)
	if descent_only:
		idx, = np.where(trajectories['alts'] == np.max(trajectories['alts']))
	else:
		idx, = np.where(trajectories['alts'] == trajectories['alts'][0])
	latx, lonx, altx, timex = trajectories['lats'][idx][0], trajectories['lons'][idx][0], trajectories['alts'][idx][0], trajectories['times'][idx][0]

	# write out file for google-earth
	other_info = [(latx, lonx, altx, 'Burst point', '%.0f minutes, %.0f meters' % (timex, altx))]
	pyb_io.save_kml(kml_fname, trajectories, other_info=other_info, params=params, parallel_err=parallel_err, perp_err=perp_err, mean_direction=np.radians(trajectories['mean_direction']))

	return trajectories

#################################################################################################################

if __name__ == "__main__":

	time0 = time.time()

	############################################################################################################

	datestr = sys.argv[1]
	utc_hour = float(sys.argv[2])
	loc0 = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
	params = None

	if len(sys.argv) > 6:
		run = sys.argv[6]
	else:
		run = None

	trajectories = runner(datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, print_verbose=False, write_verbose=True, run=run, output_figs=False)

	############################################################################################################

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write(('Program finished in ' + str(round((time.time() - time0), 1)) + ' s').ljust(60) + '\n')

#################################################################################################################
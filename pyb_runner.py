"""Script for running a single trajectory. Certain can be manually fed in or read from the param_file."""

from astropy.io import ascii
import datetime as dt
import numpy as np
import sys, os
import time

import pyb_traj
import get_gfs
import pyb_aux
import pyb_io

import param_file as p

#################################################################################################################

# method to run entire simulation of flight
# requires the starting location & starting point (0 for highest), whether or not its descent only, and the date & time of the starting point
# if parameters are not supplied, it takes it from the param file

def runner(datestr=None, utc_hour=None, loc0=None, balloon=None, params=None, run=None, print_verbose=False, write_verbose=False, add_run_info=True, output_figs=False):

	print('')
	now = dt.datetime.now()

	############################################################################################################ <---- set trajectory parameters

	# starting location
	lat0, lon0, alt0 = loc0

	# change altitude if it is underground (not always accurate)
	if pyb_aux.get_elevation(lon0, lat0) > alt0:
		alt0 = pyb_aux.get_elevation(lon0, lat0)

	# general parameters
	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

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
	trajectories, fig_dicts = pyb_traj.run_traj(weather_files=files, datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon, output_figs=output_figs)

	############################################################################################################ <---- create/write output 

	# initialise trajectory and kml file names
	traj_file = traj_dir + 'trajectory_' + datestr + '_' + str(utc_hour) + '_' + str(loc0)
	kml_fname = kml_dir + datestr + '_' + str(utc_hour) + '_' + str(loc0)

	# determine number of files already in folder (with the same starting params)
	no = len([filename for filename in os.listdir(traj_dir) if os.path.basename(traj_file) in filename])
	if no == 0:
		traj_file += '.dat'
		kml_fname += '.kml'
	else:
		traj_file += '_' + str(no + 1) + '.dat'
		kml_fname += '_' + str(no + 1) + '.kml'

	# write out trajectory file
	ascii.write([trajectories['lats'], trajectories['lons'], trajectories['alts'], trajectories['dists'], trajectories['times'], trajectories['speeds'], trajectories['z_speeds'], \
		trajectories['omegas'], trajectories['temperatures'], trajectories['grid_spreads_u'], trajectories['grid_spreads_v']], traj_file, names=['lats', 'lons', 'alts', 'dists',\
		 'times', 'speeds', 'z_speeds', 'omegas', 'temps', 'u_spread', 'v_spread'], overwrite=True)

	# write parameters of this run to file
	if write_verbose:
		pyb_io.write_verbose(params_dir=base_dir, params=params, balloon=balloon)

	# add run info to run_info.txt file
	pyb_io.write_run_info(add_run_info=add_run_info, run=run, params=params, balloon=balloon)

	# highest point in main-run trajectory (bit redundant for descent only)
	if descent_only:
		idx, = np.where(trajectories['alts'] == np.max(trajectories['alts']))
	else:
		idx, = np.where(trajectories['alts'] == trajectories['alts'][0])

	latx = trajectories['lats'][idx][0]
	lonx = trajectories['lons'][idx][0]
	altx = trajectories['alts'][idx][0]
	timex = trajectories['times'][idx][0]

	# write out file for google-earth
	other_info = [(latx, lonx, altx, 'Burst point', '%.0f minutes, %.0f meters' % (timex, altx))]
	pyb_io.save_kml(kml_fname, trajectories, other_info=other_info, params=params, radius=5)

	# save figs with interpolation checks
	if output_figs:
		for i in range(len(fig_dicts)-1):
			for key in fig_dicts[i].keys():
				fig_dicts[i][key].savefig(fig_dir + datestr + '_' + key + '_check' + str(i) + '.png')

#################################################################################################################

if __name__ == "__main__":

	time0 = time.time()

	datestr = sys.argv[1]
	utc_hour = sys.argv[2]
	loc0 = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])

	runner(datestr=datestr, utc_hour=utc_hour, loc0=loc0, print_verbose=True, write_verbose=True)

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write(('Program finished in ' + str(round((time.time() - time0), 1)) + ' s').ljust(60) + '\n')

#################################################################################################################
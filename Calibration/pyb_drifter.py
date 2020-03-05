"""
Method to determine trajectories for one starting position with different drift times.
"""

from astropy.io import ascii
import datetime as dt
import numpy as np
import sys, os
import time

import pyb_runner
import pyb_io

import param_file as p

#################################################################################################################

def drifter(ini_conditions=None, balloon=None, params=None, run=None, overwrite=False, drift_times=None):
	"""
	Method to calculate trajectories starting at the same location and time, but have different drift times

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
		If True, overwrite trajectory file with same name in folder
	drift_time: array
		Array of several drift times to check
	"""

	datestr, utc_hour, loc0 = ini_conditions
	lat0, lon0, alt0 = loc0

	out_dir = p.path + p.output_folder

	if run == None:

		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

		files = [filename for filename in os.listdir(out_dir) if now_str in filename]
		run = now_str + '_' + str(len(files))

	descent_only, drift_time, resolution, hr_diff, check_elevation, live, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	params[1] = drift_times

	pyb_io.print_verbose(datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon)

	params_dir = out_dir + run + '/'
	if not os.path.exists(params_dir):
		os.makedirs(params_dir)
	pyb_io.write_verbose(params_dir=params_dir, params=params, balloon=balloon)

	for drift_time in drift_times:

		params[1] = drift_time
		print('Running drift time: ' + str(drift_time) + ' minutes')
		pyb_runner.basic_runner(ini_conditions=ini_conditions, params=params, balloon=balloon, run=run, overwrite=overwrite)
		print('----------\n')

	pyb_io.merge_kml(datestr=datestr, run=run, params=params, balloon=balloon, drift_times=drift_times)

#################################################################################################################

if __name__ == '__main__':

	time0 = time.time()

	datestr = sys.argv[1]
	utc_hour = float(sys.argv[2])
	loc0 = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
	ini_conditions = (datestr, utc_hour, loc0)

	if len(sys.argv) > 6:
		run = sys.argv[6]
	else:
		run = None

	drift_times = np.arange(0., 40., 10.)

	drifter(ini_conditions=ini_conditions, run=run, drift_times=drift_times)

	print('Total time elapsed: %.1f s' % (time.time() - time0))

#################################################################################################################
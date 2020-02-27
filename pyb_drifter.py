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

def drifter(datestr=None, utc_hour=None, loc0=None, balloon=None, params=None, run=None, print_verbose=False, write_verbose=True, add_run_info=False, overwrite=False):
	"""
	Method to calculate trajectories starting at the same location and time, but have different drift times

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
	overwrite : bool
		If True, overwrite trajectory file with same name in folder
	"""

	out_dir = p.path + p.output_folder

	if run == None:

		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

		files = [filename for filename in os.listdir(out_dir) if now_str in filename]
		run = now_str + '_' + str(len(files))

	descent_only, drift_time, resolution, hr_diff, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	drift_times = np.arange(0., 20., 5.)
	params[4] = drift_times

	if print_verbose:
		pyb_io.print_verbose(datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon)

	if write_verbose:

		params_dir = out_dir + run + '/'
		if not os.path.exists(params_dir):
			os.makedirs(params_dir)

		pyb_io.write_verbose(params_dir=params_dir, params=params, balloon=balloon)

	for drift_time in drift_times:

		params[4] = drift_time
		print('Running drift time: ' + str(drift_time) + ' minutes')
		pyb_runner.runner(datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon, run=run, add_run_info=add_run_info, overwrite=overwrite)
		print('----------\n')

	pyb_io.merge_kml(datestr=datestr, run=run, params=params, balloon=balloon, drift_times=drift_times)

#################################################################################################################

if __name__ == '__main__':

	time0 = time.time()

	datestr = sys.argv[1]
	utc_hour = float(sys.argv[2])
	loc0 = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])

	drifter(datestr=datestr, utc_hour=utc_hour, loc0=loc0, print_verbose=True)

	print('Total time elapsed: %.1f s' % (time.time() - time0))

#################################################################################################################
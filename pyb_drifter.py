"""Method to determine trajectories for one starting position with different times"""

from astropy.io import ascii
import datetime as dt
import numpy as np
import sys, os
import time

import pyb_runner
import pyb_io

import param_file as p

#################################################################################################################

def drifter(datestr=None, utc_hour=None, loc0=None, params=None, run=None, balloon=None, print_verbose=False, write_verbose=True):

	out_dir = p.path + p.output_folder

	if loc0 != None:
		lat0, lon0, alt0 = loc0
		lat0, lon0, alt0 = float(lat0), float(lon0), float(alt0)
		loc0 = lat0, lon0, alt0

	if run == None:

		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

		files = [filename for filename in os.listdir(out_dir) if now_str in filename]
		run = now_str + '_' + str(len(files))

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, params, balloon = pyb_io.set_params(params=params, balloon=balloon)
	drift_times = np.arange(0., 60., 30.)
	params[3] = drift_times

	if utc_hour == None or loc0 == None:
		utc_hour, loc0 = pyb_io.get_ini(datestr=datestr, descent_only=descent_only, next_point=next_point)

	if print_verbose:
		pyb_io.print_verbose(datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon)

	if write_verbose:

		params_dir = out_dir + run + '/'
		if not os.path.exists(params_dir):
			os.makedirs(params_dir)

		pyb_io.write_verbose(params_dir=params_dir, params=params, balloon=balloon)

	for drift_time in drift_times:

		params[3] = drift_time
		print('Running drift time: ' + str(drift_time) + ' minutes')
		pyb_runner.runner(datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon, run=run, add_run_info=False)
		print('\n----------\n')

	pyb_io.merge_kml(datestr=datestr, run=run, params=params, balloon=balloon, drift_times=drift_times)

#################################################################################################################

if __name__ == '__main__':

	time0 = time.time()

	drifter(datestr='20181215', print_verbose=True)

	print('Total time elapsed: %.1f s' % (time.time() - time0))

#################################################################################################################
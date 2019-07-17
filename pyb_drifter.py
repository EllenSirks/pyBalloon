"""Method to determine trajectories for one starting position with different times"""

from astropy.io import ascii
import datetime as dt
import numpy as np
import sys, os
import time

import pyb_info_searcher as inf
import pyb_runner
import pyb_io

import param_file as p

def drifter(datestr=None, utc_hour=None, lat0=None, lon0=None, alt0=None, params=None, run=None, balloon=None, print_verbose=False, write_verbose=True):

	time0 = time.time()

	out_dir = '/home/ellen/Desktop/SuperBIT/Output/' 

	if run == None:

		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

		files = [filename for filename in os.listdir(out_dir) if now_str in filename]
		run = now_str + '_' + str(len(files))

	if params == None:

		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		else:
			next_point = '0'
		interpolate = p.interpolate
		resolution = p.resolution
		vz_correct = p.vz_correct
		hr_diff = p.hr_diff

	else:

		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		else:
			next_point = '0'
		interpolate = bool(params[2])
		resolution = float(params[3])
		vz_correct = bool(params[4])
		hr_diff = int(params[5])

	if balloon == None:
		balloon = p.balloon

	if utc_hour == None or lat0 == None or lon0 == None or alt0 == None:
		utc_hour, lat0, lon0, alt0 = inf.get_ini(datestr=datestr, descent_only=descent_only, next_point=next_point)

	drift_times = np.arange(0., 60., 30.)

	if print_verbose:
		pyb_io.print_verbose(datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, descent_only=descent_only, next_point=next_point, interpolate=interpolate,\
		 drift_time=drift_times, resolution=resolution, vz_correct=vz_correct, hr_diff=hr_diff, balloon=balloon)

	if write_verbose:

		params_dir = out_dir + run + '/'
		if not os.path.exists(params_dir):
			os.makedirs(params_dir)

		pyb_io.write_verbose(params_dir=params_dir, datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, descent_only=descent_only, next_point=next_point, interpolate=interpolate,\
			drift_time=drift_times, resolution=resolution, vz_correct=vz_correct, hr_diff=hr_diff, balloon=balloon)

	for drift_time in drift_times:

		params = [descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff]
		print('Running drift time: ' + str(drift_time) + ' minutes')
		pyb_runner.runner(datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, params=params, balloon=balloon, run=run, add_run_info=False)
		print('\n----------\n')

	pyb_io.merge_kml(datestr=datestr, run=run, params=params, balloon=balloon, drift_times=drift_times)

	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	drifter(datestr='20181215', print_verbose=True)
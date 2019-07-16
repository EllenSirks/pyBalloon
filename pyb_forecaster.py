from astropy.io import ascii
import datetime as dt
import numpy as np
import time
import os

import pyb_info_searcher as inf
import param_file as p
import pyb_runner
import pyb_io

def forecaster(datestr=None, utc_hour=None, loc0=None, params=None, balloon=None, run=None, print_verbose=False, write_verbose=True):

	time0 = time.time()

	lat0, lon0, alt0 = loc0

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
		drift_time = p.drift_time
		resolution = p.resolution
		vz_correct = p.vz_correct

	else:
		descent_only = params[0]
		if descent_only:
			next_point = params[1]
		else:
			next_point = '0'
		interpolate = params[2]
		drift_time = params[3]
		resolution = params[4]
		vz_correct = params[5]

	if balloon == None:
		balloon = p.balloon

	if utc_hour == None or loc0 == None:
		utc_hour, lat0, lon0, alt0 = inf.get_ini(datestr=datestr, descent_only=descent_only, next_point=next_point)

	hr_diffs = np.arange(0, 24, 6)

	if print_verbose:
		pyb_io.print_verbose(datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, descent_only=descent_only, next_point=next_point, interpolate=interpolate,\
		 drift_time=drift_time, resolution=resolution, vz_correct=vz_correct, hr_diff=hr_diffs, balloon=balloon)

	if write_verbose:

		params_dir = out_dir + run + '/'
		if not os.path.exists(params_dir):
			os.makedirs(params_dir)

		pyb_io.write_verbose(params_dir=params_dir, datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, descent_only=descent_only, next_point=next_point, interpolate=interpolate,\
			drift_time=drift_time, resolution=resolution, vz_correct=vz_correct, hr_diff=hr_diffs, balloon=balloon)

	for hr_diff in hr_diffs:

		params = [descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff]
		print('Running hr_diff: ' + str(hr_diff) + ' hours')
		pyb_runner.runner(datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, params=params, balloon=balloon, run=run, add_run_info=False)
		print('\n----------\n')

	print('Total time elapsed: %.1f s' % (time.time() - time0))

# def forecast_tester_looper(datestr=None, utc_hour=None, loc0=None, params=None, balloon=None, run=None, print_verbose=False):

# 	time0 = time.time()

# 	day = int(datestr[6:])
# 	for i in range(0, 30):
# 		datestr = datestr[:6] + str(day + i).zfill(2)
# 		forecast_tester(datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon, run=run, print_verbose=print_verbose)

# 	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	forecaster(datestr='20180901', utc_hour=14, loc0=(48.5, -81.4, 28000), print_verbose=True)
	# forecast_tester_looper(datestr='20190901', verbose=True)
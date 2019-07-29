"""Function for looping over several dates/locations and calculating their trajectories"""

from astropy.io import ascii
import datetime as dt
import time

import pyb_runner

import param_file as p

def looper(params=None, balloon=None, run=None, print_verbose=False, output_figs=False):

	time0 = time.time()

	if run == None:

		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

		files = [filename for filename in os.listdir(out_dir) if now_str in filename]
		run = now_str + '_' + str(len(files))

	dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

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
		hr_diff = p.hr_diff

		params = [descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff]

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
		hr_diff = params[6]

	if balloon == None:
		balloon = p.balloon

	if descent_only:
		fname = 'descent_only_start' + next_point + '.txt'
	else:
		print('No file with specified flight data!')

	if print_verbose:
		pyb_io.print_verbose(datestr=datestr, utc_hour=utc_hour, lat0=lat0, lon0=lon0, alt0=alt0, descent_only=descent_only, next_point=next_point, interpolate=interpolate,\
		 drift_time=drift_time, resolution=resolution, vz_correct=vz_correct, hr_diff=hr_diff, balloon=balloon)

	print('\nRunning file: ' + fname)

	lines = [line.rstrip('\n').split(' ') for line in open(dir + fname)]
	for i in range(len(lines)):

		print('\n----------')
		print('\nRunning date: ' + lines[i][0])
		print('Starting point: ' + str(lines[i][2]) + ' lat., ' + str(lines[i][3]) + ' lon., ' + str(lines[i][4]) + ' m')

		try: 
			pyb_runner.runner(datestr=lines[i][0], utc_hour=lines[i][1], lat0=lines[i][2], lon0=lines[i][3], alt0=lines[i][4], params=params, balloon=balloon, run=run, write_verbose=True, output_figs=output_figs)
		except Exception as e: 
			print(e)
			continue

	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	looper(print_verbose=True)
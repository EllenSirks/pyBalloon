"""Pipeline for performing all analysis in one go"""

from astropy.io import ascii
import datetime as dt
import time, os

import pyb_plotter
import pyb_runner
import pyb_io

import param_file as p

def looper(params=None, balloon=None, run=None, print_verbose=False, output_figs=False):

	if run == None:

		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

		files = [filename for filename in os.listdir(out_dir) if now_str in filename]
		run = now_str + '_' + str(len(files))

	flights_dir = p.path + 'Flight_data/'

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	if descent_only:
		fname = 'descent_only_start' + next_point + '.txt'
		if not os.path.isfile(flights_dir + fname):
			print('No file with specified flight data!')
	else:
		print('No file with specified flight data!')

	if print_verbose:
		pyb_io.print_verbose(params=params, balloon=balloon)

	print('\nRunning file: ' + fname)

	lines = [line.rstrip('\n').split(' ') for line in open(flights_dir + fname)]
	for i in range(len(lines)):

		print('\n----------')
		print('\nRunning date: ' + lines[i][0])
		print('Starting point: ' + str(lines[i][2]) + ' lat., ' + str(lines[i][3]) + ' lon., ' + str(lines[i][4]) + ' m')

		try: 
			pyb_runner.runner(datestr=lines[i][0], utc_hour=lines[i][1], lat0=lines[i][2], lon0=lines[i][3], alt0=lines[i][4], params=params, balloon=balloon, run=run, write_verbose=True, output_figs=output_figs)
		except Exception as e: 
			print(e)
			continue

def pipeline(params=None, balloon=None, print_verbose=False, output_figs=False):

	time0 = time.time()

	out_dir = p.path + 'Output/'

	now = dt.datetime.now()
	now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

	files = [filename for filename in os.listdir(out_dir) if now_str in filename]
	run = now_str + '_' + str(len(files))

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	if print_verbose:
		print_verbose(params=params, balloon=balloon)

	looper(params=params, balloon=balloon, run=run, output_figs=output_figs)
	pyb_plotter.plot_rates(params=params, run=run, all_plots=False)
	pyb_plotter.plot_results(params=params, run=run)

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write(('Program finished in ' + str(round((time.time() - time0), 1)) + ' s').ljust(60) + '\n')

if __name__ == '__main__':

	pipeline(print_verbose=True, output_figs=True)
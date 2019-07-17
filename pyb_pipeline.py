"""Pipeline for performing all analysis in one go"""

import datetime as dt
import time
import os

import pyb_plotter
import pyb_looper

import pyb_io

import param_file as p

def pipeline(params=None, balloon=None, print_verbose=False):

	time0 = time.time()

	out_dir = '/home/ellen/Desktop/SuperBIT/Output/'

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
		hr_diff = p.hr_diff

		params = [descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff]

	else:
		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		else:
			next_point = '0'
		interpolate = bool(params[2])
		drift_time = float(params[3])
		resolution = float(params[4])
		vz_correct = bool(params[5])
		hr_diff = int(params[6])

	if balloon == None:
		balloon = p.balloon

	if print_verbose:
		print('\nGeneral Parameters')
		print('----------')
		print('descent_only: ' + str(descent_only))
		if descent_only:
			print('starting point: ' + next_point)
		print('interpolate: ' + str(interpolate))
		print('drift time: ' + str(drift_time) + ' minutes')
		print('resolution of forecasts: ' + str(resolution) + ' degrees')
		print('correct for vertical winds: ' + str(vz_correct))
		print('difference in hrs for forecasts: ' + str(hr_diff) + ' hours')
		print('----------')
		print('\nBalloon/Parachute Parameters')
		print('----------')
		print('altitude step: ' + str(balloon['altitude_step']) + ' m')
		print('equipment mass: ' + str(balloon['equip_mass']) + ' kg')
		print('parachute Cd: ' + str(balloon['Cd_parachute']))
		print('parachute area: ' + str(round(balloon['parachute_areas'][0], 2)) + ' m^2')
		print('----------')

	pyb_looper.looper(params=params, balloon=balloon, run=run)
	pyb_plotter.plot_rates(params=params, run=run, all_plots=False)
	pyb_plotter.plot_results(params=params, run=run)

	print('Program finished. Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	pipeline(print_verbose=True)
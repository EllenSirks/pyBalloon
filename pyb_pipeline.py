""" pipeline for performing all analysis in one go """

import datetime as dt
import time
import os

import pyb_plotter
import pyb_looper

import param_file as p

def pipeline(params=None, balloon=None):

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

		params = [descent_only, next_point, interpolate, drift_time]

	else:
		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		else:
			next_point = '0'
		interpolate = bool(params[2])
		drift_time = float(params[3])

	if balloon == None:
		balloon = p.balloon

	print('\nGeneral Parameters')
	print('----------')
	print('descent_only: ' + str(descent_only))
	if descent_only:
		print('starting point: ' + next_point)
	print('interpolate: ' + str(interpolate))
	print('drift time: ' + str(drift_time) + ' minutes')
	print('----------')
	print('\nBalloon/Parachute Parameters')
	print('----------')
	print('altitude step: ' + str(balloon['altitude_step']) + ' m')
	print('equipment mass: ' + str(balloon['equip_mass']) + ' kg')
	print('parachute Cd: ' + str(balloon['Cd_parachute']))
	print('parachute area: ' + str(round(balloon['parachute_areas'][0], 2)) + ' m^2')
	print('----------')

	pyb_looper.looper(params=params, balloon=balloon, run=run)
	pyb_plotter.plot_rates(params=params, run=run)
	pyb_plotter.plot_results(params=params, run=run)

	print('Program finished. Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	pipeline()
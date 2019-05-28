""" pipeline for performing all analysis in one go """

import time

import pyb_property_plotter
import pyb_results_plotter
import param_file as p
import pyb_looper

def pipeline(params=None, balloon=None):

	time0 = time.time()

	if params == None:

		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		interpolate = p.interpolate
		drift_time = p.drift_time

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

	pyb_looper.looper(params=params, balloon=balloon)
	pyb_property_plotter.plot_rates(params=params)
	pyb_results_plotter.plot_results(params=params)

	print('Program finished. Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	pipeline()
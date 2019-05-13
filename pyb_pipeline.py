import time

import pyb_property_plotter
import pyb_results_plotter
import param_file as p
import pyb_looper

def pipeline(params=None):

	time0 = time.time()

	if params == None:

		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		interpolate = p.interpolate
		drift_time = p.drift_time

	print('\nParameters')
	print('----------')
	print('descent_only: ' + str())
	if descent_only:
		print('starting point: ' + next_point)
	print('interpolate: ' + str(interpolate))
	print('drift time: ' + str(drift_time) + ' minutes')
	print('----------')

	pyb_looper.looper(params=params)
	pyb_property_plotter.plot_rates()
	pyb_results_plotter.plot_results()

	print('Program finished. Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	pipeline()
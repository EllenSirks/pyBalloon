from astropy.io import ascii
import numpy as np
import sys

import param_file as p

def search_info(run=None, verbose=True):

	data = ascii.read('/home/ellen/Desktop/SuperBIT/Output/runs_info.txt')
	runs = data['run']

	index = np.where(runs == run)[0][0]

	if verbose:

		print('\nGeneral Parameters')
		print('----------')
		print('descent_only: ' + str(data['descent_only'][index]))
		if data['descent_only'][index]:
			print('starting point: ' + str(data['next_point'][index]))
		print('interpolate: ' + str(data['interpolate'][index]))
		print('drift time: ' + str(data['drift_time'][index]) + ' minutes')
		print('resolution of forecasts: ' + str(data['resolution'][index]) + ' degrees')
		print('----------')
		print('\nBalloon/Parachute Parameters')
		print('----------')
		print('altitude step: ' + str(data['altitude_step'][index]) + ' m')
		print('equipment mass: ' + str(data['equip_mass'][index]) + ' kg')
		print('parachute Cd: ' + str(data['Cd_parachute'][index]))
		print('parachute radius: ' + str(data['parachute_radius'][index]) + ' m')
		print('parachute area: ' + str(round(np.pi*data['parachute_radius'][index]**2, 2)) + ' m^2')
		print('----------\n')

def seach_info_drift_times(drift_times=None, params=None, balloon=None):

	data = ascii.read('/home/ellen/Desktop/SuperBIT/Output/runs_info.txt')
	runs = data['run']
	times = data['drift_time']
	descent_only = data['descent_only']
	next_point = np.array(data['next_point'])
	interpolate  = data['interpolate']
	resolution = data['resolution']
	Cd_parachute = data['Cd_parachute']
	parachute_radius = data['parachute_radius']
	altitude_step = data['altitude_step']
	equip_mass = data['equip_mass']

	if params == None:

		p_descent_only = str(p.descent_only)
		if p.descent_only:
			p_next_point = int(p.next_point)
		else:
			p_next_point = 0
		p_interpolate = str(p.interpolate)
		p_resolution = str(p.resolution)

	else:

		p_descent_only = str(params[0])
		if p_descent_only:
			p_next_point = int(params[1])
		else:
			p_next_point = 0
		p_interpolate = str(params[2])
		p_resolution = str(params[3])

	if balloon == None:

		p_Cd_parachute = float(p.balloon['Cd_parachute'])
		p_parachute_radius = float(np.sqrt(p.balloon['parachute_areas'][0]/np.pi))
		p_altitude_step = float(p.balloon['altitude_step'])
		p_equip_mass = float(p.balloon['equip_mass'])

	else:

		p_Cd_parachute = float(balloon['Cd_parachute'])
		p_parachute_radius = float(np.sqrt(balloon['parachute_areas'][0]/np.pi))
		p_altitude_step = float(balloon['altitude_step'])
		p_equip_mass = float(balloon['equip_mass'])

	folders = {}

	for drift_time in drift_times:

		index = np.where((times == drift_time) & (descent_only == p_descent_only) & (next_point == p_next_point) & (interpolate == p_interpolate) & (resolution == p_resolution) & (Cd_parachute == p_Cd_parachute) & (parachute_radius == p_parachute_radius) & \
		 (altitude_step == p_altitude_step) &  (equip_mass == p_equip_mass))[0][-1]

		folders[drift_time] = runs[index]

	return folders

if __name__ == '__main__':

	run = str(sys.argv[1])
	search_info(run=run, verbose=True)

	# drift_times = np.arange(0., 120., 30.)
	# folders = seach_info_drift_times(drift_times=drift_times)





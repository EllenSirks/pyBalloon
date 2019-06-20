from astropy.io import ascii
import numpy as np
import sys

def search_info(run=None):

	data = ascii.read('/home/ellen/Desktop/SuperBIT/Output/runs_info.txt')
	runs = data['run']

	index = np.where(runs == run)[0][0]

	print('\nGeneral Parameters')
	print('----------')
	print('descent_only: ' + str(data['descent_only'][index]))
	if data['descent_only'][index]:
		print('starting point: ' + str(data['next_point'][index]))
	print('interpolate: ' + str(data['interpolate'][index]))
	print('drift time: ' + str(data['drift_time'][index]) + ' minutes')
	print('----------')
	print('\nBalloon/Parachute Parameters')
	print('----------')
	print('altitude step: ' + str(data['altitude_step'][index]) + ' m')
	print('equipment mass: ' + str(data['equip_mass'][index]) + ' kg')
	print('parachute Cd: ' + str(data['Cd_parachute'][index]))
	print('parachute radius: ' + str(data['parachute_radius'][index]) + ' m')
	print('parachute area: ' + str(round(np.pi*data['parachute_radius'][index]**2, 2)) + ' m^2')
	print('----------\n')

if __name__ == '__main__':

	run = str(sys.argv[1])
	search_info(run=run)

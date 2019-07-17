"""Functions for getting information/parameters from specific run or date"""

from astropy.io import ascii
import numpy as np
import sys

def get_ini(datestr=None, descent_only=False, next_point='1'):

	dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

	if descent_only:
		fname = 'descent_only_start' + next_point + '.txt'
	else:
		print('No file with specified flight data!')
		sys.exit()

	lines = [line.rstrip('\n').split(' ') for line in open(dir + fname)]
	for i in range(len(lines)):
		if lines[i][0] == datestr:
			break

	return lines[i][1], lines[i][2], lines[i][3], lines[i][4]

def search_info(run=None, print_verbose=True):

	data = ascii.read('/home/ellen/Desktop/SuperBIT/Output/runs_info.txt')
	runs = data['run']

	index = np.where(runs == run)[0]

	if len(index) < 1:

		print('Cannot find the info for this run! Check the folder.')
		sys.exit()

	else:

		index = int(index)

		if print_verbose:

			print('\nGeneral Parameters')
			print('----------')
			print('descent_only: ' + str(data['descent_only'][index]))
			if data['descent_only'][index]:
				print('starting point: ' + str(data['next_point'][index]))
			print('interpolate: ' + str(data['interpolate'][index]))
			print('drift time: ' + str(data['drift_time'][index]) + ' minutes')
			print('resolution of forecasts: ' + str(data['resolution'][index]) + ' degrees')
			print('correct for vertical winds: ' + str(data['vz_correct'][index]))
			print('difference in hrs for forecasts: ' + str(data['hr_diff'][index]) + ' hours')
			print('----------')
			print('\nBalloon/Parachute Parameters')
			print('----------')
			print('altitude step: ' + str(data['altitude_step'][index]) + ' m')
			print('equipment mass: ' + str(data['equip_mass'][index]) + ' kg')
			print('parachute Cd: ' + str(data['Cd_parachute'][index]))
			print('parachute area: ' + str(data['parachute_area'][index]) + ' m^2')
			print('----------\n')

		return index

if __name__ == '__main__':

	run = str(sys.argv[1])
	index = search_info(run=run)
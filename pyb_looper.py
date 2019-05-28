from astropy.io import ascii
import time

from pyb_runner import run
import param_file as p

def looper(params=None, balloon=None, verbose=False):

	time0 = time.time()

	dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

	if params == None:

		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		interpolate = p.interpolate
		drift_time = p.drift_time

	else:

		descent_only = params[0]
		if descent_only:
			next_point = params[1]
		interpolate = params[-2]
		drift_time = params[-1]

	if descent_only:
		fname = 'descent_only_start' + next_point + '.txt'

	else:
		print('No file with specified flight data!')

	if balloon == None:
		balloon = p.balloon

	if verbose:

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

	print('\nRunning file: ' + fname)

	lines = [line.rstrip('\n').split(' ') for line in open(dir + fname)]
	for i in range(len(lines)):
		print('\n----------')
		print('\nRunning date: ' + lines[i][0])
		print('Starting point: ' + str(lines[i][2]) + ' lat., ' + str(lines[i][2]) + ' lon., ' + str(lines[i][4]) + ' m')
		try: 
			run(datestr=lines[i][0], utc_hour=lines[i][1], lat0=lines[i][2], lon0=lines[i][3], alt0=lines[i][4], params=params, balloon=balloon)
		except Exception as e: 
			print(e)
			continue

	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	looper(verbose=True)
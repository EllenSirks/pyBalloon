from astropy.io import ascii
import datetime as dt
import time

from pyb_runner import runner
import param_file as p

def looper(params=None, balloon=None, run=None, verbose=False):

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

		params = [descent_only, next_point, interpolate, drift_time, resolution]

	else:
		descent_only = params[0]
		if descent_only:
			next_point = params[1]
		else:
			next_point = '0'
		interpolate = params[2]
		drift_time = params[3]
		resolution = params[4]

	if balloon == None:
		balloon = p.balloon

	if descent_only:
		fname = 'descent_only_start' + next_point + '.txt'
	else:
		print('No file with specified flight data!')

	if verbose:

		print('\nGeneral Parameters')
		print('----------')
		print('descent_only: ' + str(descent_only))
		if descent_only:
			print('starting point: ' + next_point)
		print('interpolate: ' + str(interpolate))
		print('drift time: ' + str(drift_time) + ' minutes')
		print('resolution of forecasts: ' + str(resolution) + ' degrees')
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
		print('Starting point: ' + str(lines[i][2]) + ' lat., ' + str(lines[i][3]) + ' lon., ' + str(lines[i][4]) + ' m')

		try: 
			runner(datestr=lines[i][0], utc_hour=lines[i][1], lat0=lines[i][2], lon0=lines[i][3], alt0=lines[i][4], params=params, balloon=balloon, run=run)
		except Exception as e: 
			print(e)
			continue

	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	looper(verbose=True)
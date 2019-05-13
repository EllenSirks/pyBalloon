from astropy.io import ascii
import time

from pyb_runner import run
import param_file as p

def looper(params=None):

	time0 = time.time()

	dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

	if params == None:

		descent_only == p.descent_only:	
		next_point = p.next_point

	else:

		descent_only = params[0]
		if descent_only:
			next_point = params[1]

	if descent_only:
		fname = 'descent_only_start' + next_point + '.txt'

	else:
		print('No file with specified flight data!')

	print('\nRunning file: ' + fname)

	lines = [line.rstrip('\n').split(' ') for line in open(dir + fname)]
	for i in range(len(lines)):
		try: 
			run(datestr=lines[i][0], utc_hour=lines[i][1], lat0=lines[i][2], lon0=lines[i][3], alt0=lines[i][4], params=params)
		except Exception as e: 
			print(e)
			continue

	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	looper()
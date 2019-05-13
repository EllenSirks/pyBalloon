from astropy.io import ascii
import time

from pyb_runner import run
import param_file as p

time0 = time.time()

dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

if p.descent_only:

	next_point = p.next_point

	fname = 'descent_only_start' + next_point + '.txt'

	print('\nRunning file: ' + fname)

else:

	print('There is no file with ascent start positions!')

lines = [line.rstrip('\n').split(' ') for line in open(dir + fname)]
for i in range(len(lines)):
	try: 
		run(datestr=lines[i][0], utc_hour=lines[i][1], lat0=lines[i][2], lon0=lines[i][3], alt0=lines[i][4])
	except Exception as e: 
		print(e)
		continue

print('Total time elapsed: %.1f s' % (time.time() - time0))
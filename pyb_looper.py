from run_pyBalloon import run
from astropy.io import ascii
import time

interpolate = False
drift_time = 60.
descent_only = True

print('Options:\n')

dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

fname1 = "flights.txt"
fname2 = "flights_next.txt"
fname3 = "flights_next2.txt"

print('Option 1: ' + fname1)
print('Option 2: ' + fname2)
print('Option 3: ' + fname3 + '\n')

option = input('Please choose an option: ')

time0 = time.time()

if option == 1:
	fname = fname1
elif option == 2:
	fname = fname2
else:
	fname = fname3

print('\nRunning file: ' + fname)

if 'next' in fname:
	if 'next2' in fname:
		next_point = '2'
	else:
		next_point = '1'
else:
	next_point = '0'

lines = [line.rstrip('\n').split(' ') for line in open(dir + fname)]
for i in range(len(lines)):
	try: 
		run(datestr=lines[i][0], utc_hour=lines[i][1], lat0=lines[i][2], lon0=lines[i][3], alt0=lines[i][4], descent_only=descent_only, next_point=next_point, interpolate=interpolate, drift_time=drift_time)
	except Exception as e: 
		print(e)
		continue

print('Total time elapsed: %.1f s' % (time.time() - time0))
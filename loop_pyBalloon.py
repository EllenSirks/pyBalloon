from retrospective import retrospective
from astropy.io import ascii
import time

print('Options:\n')

dir = '/home/ellen/Desktop/SuperBIT/Flight_data/'

fname1 = "flights.txt"
fname2 = "flights_precise.txt"
fname3 = "flights+weather_files.txt"
fname4 = "flights_next+weather_files.txt"
fname5 = "flights_next2+weather_files.txt"

print('Option 1: ' + fname1)
print('Option 2: ' + fname2)
print('Option 3: ' + fname3)
print('Option 4: ' + fname4)
print('Option 5: ' + fname5 + '\n')

option = input('Please choose an option: ')

time0 = time.time()

if option == 1:
	fname = fname1
elif option == 2:
	fname = fname2
elif option == 3:
	fname = fname3
elif option == 4:
	fname = fname4
else:
	fname = fname5

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
		if option < 3:
			retrospective(datestr=lines[i][0], utc_hour=lines[i][1], lat0=lines[i][2], lon0=lines[i][3], alt0=lines[i][4], descent_only='True', next_point=next_point)
		else:
			retrospective(weather_file=lines[i][0], lat0=lines[i][1], lon0=lines[i][2], alt0=lines[i][3], descent_only='True', next_point=next_point) 
	except Exception as e: 
		print(e)
		continue

print('Total time elapsed: %.1f s' % (time.time() - time0))
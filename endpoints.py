from astropy.io import ascii
import elevation as el
import numpy as np
import time
import csv
import sys
import os

def get_file(date, fext='all_points/'):

	dir_pred = '/home/ellen/Desktop/SuperBIT/Weather_data/Trajectories/' + fext

	no = 1

	files = [file for file in os.listdir(dir_pred) if file.startswith(date)]

	if len(files) > 1:
		no = int(input('Flight no.: '))

	t1 = np.array([int(int(file[9:13])/100 + int(file[14:17])) for file in files])
	t2 = np.array([int(int(file[9:13])/100 + int(file[14:17])) for file in files])

	t1.sort()

	ind = int(np.where(t2 == t1[no - 1])[0])
	filename = files[ind]

	return filename, fext

def get_endpoint(filename, fext='all_points/'):

	if len(filename) == 8:
		filename, fext = get_file(date=filename, fext=fext)

	dir_pred = '/home/ellen/Desktop/SuperBIT/Weather_data/Trajectories/' + fext

	data  = ascii.read(dir_pred + filename)

	lats = data['lats']
	lons = data['lons']
	alts = data['alts']
	dists = data['dists']

	elevations = []

	print('Checking altitudes & elevations...')

	for i in range(1, len(lats)):
		elevation = el.get_elevation(lats[-i], lons[-i])
		elevations.append(elevation)
		if elevation < alts[-i]:
			break

	inds = (-i, -i + 1)

	dlon = lons[inds[0]] - lons[inds[1]]
	dlat = lats[inds[0]] - lats[inds[1]]

	x1 = dists[inds[1]]

	y1 = alts[inds[0]] - alts[inds[1]]
	y2 = elevations[-2] - alts[inds[1]]

	dx = x1*(y2/y1)
	x2 = x1 - dx
	f = dx/x1

	newlon = lons[inds[1]] + f*dlon
	newlat = lats[inds[1]] + f*dlat

	print('Old end-point: (' + str(lats[-1]) + ', ' + str(lons[-1]) + ')' )
	print('New end-point: (' + str(newlat) + ', ' + str(newlon) + ')' )

	f = open('/home/ellen/Desktop/SuperBIT/Weather_data/Endpoints/Precise/' + fext + filename[:-15] +'_endpoint_precise.dat','w+')
	f.write('lat lon\n')
	f.write(str(lats[-1]) + ' ' + str(lons[-1]) + '\n')
	f.write(str(newlat) + ' ' + str(newlon))
	f.close()

	return (newlat, newlon)

if __name__ == '__main__':

	time0 = time.time()

	next_point = input('Start at point after max. altitude: ')

	if next_point == '2':
		fext = 'next2_points/'
	elif next_point == '1':
		fext = 'next_points/'
	else:
		fext = 'all_points/'

	new_coords = get_endpoint(filename=sys.argv[1], fext=fext)
	print('%.1f s elapsed' % (time.time() - time0))

	
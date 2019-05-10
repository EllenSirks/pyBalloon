from astropy.io import ascii
import elevation as el
import numpy as np
import time
import csv
import sys
import os

# method to get more accurate endpoint for predictions as they can go underground.
def get_endpoint(filename, next_point='0'):

	fext = 'start_point' + next_point + '/'

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

	f = open('/home/ellen/Desktop/SuperBIT/Weather_data/Endpoints/' + fext + 'Precise/' + filename[:-15] +'_endpoint_precise.dat','w+')
	f.write('lat lon\n')
	f.write(str(lats[-1]) + ' ' + str(lons[-1]) + '\n')
	f.write(str(newlat) + ' ' + str(newlon))
	f.close()

	return (newlat, newlon)

if __name__ == '__main__':

	time0 = time.time()

	next_point = input('Start at point after max. altitude: ')
	new_coords = get_endpoint(filename=sys.argv[1], next_point=next_point)

	print('%.1f s elapsed' % (time.time() - time0))

	
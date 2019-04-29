from astropy.io import ascii
import numpy as np
import gdal
import time
import sys
import os

srtm_dir = '/home/ellen/Desktop/SuperBIT/SRTM_data/'

def find_srtm_file(lat, lon):

	data = ascii.read(srtm_dir + 'srtm_data_limits.txt')
	files, tlat, blat, llon, rlon = data['file'], data['latN'], data['latS'], data['lonW'], data['lonE']

	file = None

	for i in range(len(files)):

		if lat < tlat[i] and lat > blat[i]:
			if lon < rlon[i] and lon > llon[i]:
				file = files[i]

	if file != None:
		return str(file)
	else:
		print('Correct SRTM file data is not here!')
		return


def find_loc(lat, lon, srtm_file):

	file = srtm_dir + srtm_file

	ds = gdal.Open(file)

	band = ds.GetRasterBand(1)
	elevations = band.ReadAsArray()
	nrows, ncols = elevations.shape
	x0, dx, dxdy, y0, dydx, dy = ds.GetGeoTransform()

	lons = np.arange(x0, x0 + dx*(ncols-1) + dx, dx) 
	lon_arr = np.array([lons for i in range(nrows)])
	lats = np.arange(y0, y0 + dy *nrows, dy)                                                   
	lat_arr = np.array([np.array([l for i in range(ncols)]) for l in lats])

	diff = np.sqrt((lon_arr - lon)**2 + (lat_arr - lat)**2)
	min_diff = np.min(diff)
	i0, i1 = np.where(diff == min_diff)

	elevation = elevations[int(i0)][int(i1)]

	# elevations, diff, lons, lats = None, None, None, None

	return min_diff, elevation


def get_elevation(lat, lon, srtm_file=None):

	# print('Getting elevation at: (' + str(lat) + ', ' + str(lon) + ')')

	lat, lon = float(lat), float(lon)

	if lon < 0:
		lon = lon + 360 % 360

	if srtm_file == None:
		srtm_file = find_srtm_file(lat, lon)

	min_diff, elevation0 = find_loc(lat, lon, srtm_file)

	# print('Difference in location: ' + str(min_diff) + ' degrees')
	# print('Elevation: ' + str(elevation0) + ' m')

	return elevation0

if __name__ == "__main__":

	time0 = time.time()

	try:
		if len(sys.argv) == 4:
			print('Getting elevation at: (' + sys.argv[1] + ', ' + sys.argv[2] + ')')
			elevation = get_elevation(lat=sys.argv[1], lon=sys.argv[2], srtm_file=sys.argv[3])
			print('Elevation: ' + str(elevation) + ' m')

		elif len(sys.argv) == 3:
			print('Getting elevation at: (' + sys.argv[1] + ', ' + sys.argv[2] + ')')
			elevation = get_elevation(lat=sys.argv[1], lon=sys.argv[2])
			print('Elevation: ' + str(elevation) + ' m')

	except:

		print('Not enough information to run...')
		sys.exit()

	print('%.1f s elapsed' % (time.time() - time0))
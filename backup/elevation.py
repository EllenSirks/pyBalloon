from astropy.io import ascii
import numpy as np
import gdal
import time
import sys
import os

# methods to find elevation at given lon, lat


# method to find the srtm file with the correct lon/lat limits
def find_srtm_file(lat, lon):

	srtm_dir = '/home/ellen/Desktop/SuperBIT/SRTM_data/'

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

# find the elevation at the location (the loc closest in the grid)
def get_elevation(lat, lon, srtm_file):

	srtm_dir = '/home/ellen/Desktop/SuperBIT/SRTM_data/'

	if lon < 0:
		lon = lon + 360 % 360

	if srtm_file == None:
		srtm_file = find_srtm_file(lat, lon)

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

	elevations, diff, lons, lats = None, None, None, None

	return elevation
"""
Input and output functions used by pyBalloon
"""

from shapely.geometry import Point
from shapely.ops import transform
from functools import partial
from scipy import interpolate
from astropy.io import ascii
import pygrib as pg
import numpy as np
import sys, os
import pyproj
import math
import time

import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "serif"

import pyb_aux
import pyb_io
import pyb_traj

import param_file as p

#################################################################################################################
#################################################################################################################

def read_gfs_file(fname, area=None, alt0=0, t_0=None, descent_only=False):

	"""
	Method to read GFS file and collect relevant data (winds, temperatures, altitudes)

	Arguments
	=========
	fname : string
		The name of the grib file to be read
	area : (float, float, float, float)
		The (most northern latitude, most western longitude, most southern latitude, most eastern longitude) in degrees, indicating the area on the globe from which the data should be collected
	alt0 : float
		Initial altitude in km
	t_0 : float
		Initial temperature.
	descent_only : bool
		Option to start the trajectory at its highest point. If True, the trajectory only has a descent phase
	"""

	if area is not None:
		tlat, llon, blat, rlon = area
	else:
		print('Do you really wish to search the entire planet?')
		tlat, llon, blat, rlon = 90., 0., -90., 360.

	############################################################################################################

	grib = pg.open(fname)
	grib.seek(0)
	u_msgs = grib.select(name='U component of wind')
	v_msgs = grib.select(name='V component of wind')
	g_msgs = grib.select(name='Geopotential Height')
	t_msgs = grib.select(name='Temperature')

	lats, lons, = u_msgs[0].latlons() # lats: -90 -> 90, lons: 0 -> 360
	lats2, lons2 = u_msgs[0].latlons()

	for i in range(len(lons2)):
		for j in range(len(lons2[i])):
			if lons2[i][j] > 180:
				lons2[i][j] -= 360

	# Find closest pixel location
	locs = pyb_aux.all_and([lats <= tlat, lats >= blat, lons2 <= rlon, lons2 >= llon])

	row_idx, col_idx = np.where(locs)
	lats = lats[row_idx, col_idx]
	lons = lons[row_idx, col_idx]

	if len(lats) == 0: print('Warning! lats is empty!')
	if len(lons) == 0: print('Warning! lons is empty!')

	############################################################################################################

	# Collect U component of wind data
	u_wind = {}
	for msg in u_msgs:
		if msg.typeOfLevel == 'isobaricInhPa':
			u_wind[msg.level] = msg.values[row_idx, col_idx]
	
	# Collect V component of wind data
	v_wind = {}
	for msg in v_msgs:
		if msg.typeOfLevel == 'isobaricInhPa':
			v_wind[msg.level] = msg.values[row_idx, col_idx]

	# Collect temperatures
	temperature = {}
	for msg in t_msgs:
		if msg.typeOfLevel == 'isobaricInhPa':
			temperature[msg.level] = msg.values[row_idx, col_idx]

		# Add msg.typeOfLevel == 'surface', save to another variable for later use
		if msg.typeOfLevel == 'surface':
			t_surface = msg.values[row_idx, col_idx]

	# Collect Geopotential heights
	altitude = {}
	for msg in g_msgs:
		if msg.typeOfLevel == 'isobaricInhPa':
			altitude[msg.level] = msg.values[row_idx, col_idx]

	############################################################################################################

	pressures = list(u_wind.keys())

	# Collect data to correct altitude order. Set "surface" values before real data.
	# Use given surface temperature if available, otherwise use the model value
	if not descent_only:

		u_winds = [np.zeros(lats.shape)]
		v_winds = [np.zeros(lats.shape)]
		altitudes = [alt0*np.ones(lats.shape)]

		if t_0 is None:
			temperatures = [t_surface]
		else:
			temperatures = [t_0*np.ones(lats.shape)]

		pressures.append(max(pressures))

		ind = 0

	else:

		u_winds = []
		v_winds = []
		temperatures = []
		altitudes = []

		alts_min = np.array([np.min(alt) for alt in altitude.values()])

		diff = np.abs(alts_min - alt0)
		grid_i1, = np.where(diff == diff.min())
		grid_i1 = grid_i1[0]

		p_or = list(altitude.keys())[grid_i1]
		pressures.append(p_or)

	# Put pressures in altitude order and use them as keys
	pressures.sort()
	pressures.reverse()

	if descent_only:
		ind = np.where(np.array(pressures) == p_or)[0][0]

	############################################################################################################	

	i = 0
	for key in pressures:
		if i != ind:
			uwnd, vwnd, temp, alt, omg = [], [], [], [], []
			uwnd.append(u_wind[key])
			vwnd.append(v_wind[key])
			temp.append(temperature[key])
			alt.append(altitude[key])

			u_winds.append(np.hstack(uwnd))
			v_winds.append(np.hstack(vwnd))
			temperatures.append(np.hstack(temp))
			altitudes.append(np.hstack(alt))

			i+=1
		else:
			i+=1

	############################################################################################################

	if descent_only:

		alts_mean = np.array([np.mean(alt) for alt in altitudes])
		alts_min = np.array([np.min(alt) for alt in altitudes])

		diff = np.abs(alts_min - alt0)
		grid_i2, = np.where(diff == diff.min())
		grid_i2 = grid_i2[0]

		if alts_min[grid_i2] > alt0:
			
			deltaA1 = altitudes[grid_i2] - altitudes[grid_i2 - 1]
			deltaA2 = alt0 - altitudes[grid_i2 - 1]

			f1 = deltaA2/deltaA1

			deltaT = temperatures[grid_i2] - temperatures[grid_i2 - 1]
			deltaU = u_winds[grid_i2] - u_winds[grid_i2 - 1]
			deltaV = v_winds[grid_i2] - v_winds[grid_i2 - 1]

			altitudes.insert(grid_i2, alt0*np.ones(lats.shape))
			u_winds.insert(grid_i2, u_winds[grid_i2 - 1] + deltaU*f1)
			v_winds.insert(grid_i2, v_winds[grid_i2 - 1] + deltaV*f1)
			temperatures.insert(grid_i2, temperatures[grid_i2 - 1] + deltaT*f1)
			index = grid_i2

		else:

			deltaA1 = altitudes[grid_i2 + 1] - altitudes[grid_i2]
			deltaA2 = alt0 - altitudes[grid_i2]

			f1 = deltaA2/deltaA1

			deltaT = temperatures[grid_i2 + 1] - temperatures[grid_i2]
			deltaU = u_winds[grid_i2 + 1] - u_winds[grid_i2]
			deltaV = v_winds[grid_i2 + 1] - v_winds[grid_i2]

			altitudes.insert(grid_i2 + 1, alt0*np.ones(lats.shape))
			u_winds.insert(grid_i2 + 1, u_winds[grid_i2] + deltaU*f1)
			v_winds.insert(grid_i2 + 1, v_winds[grid_i2] + deltaV*f1)
			temperatures.insert(grid_i2 + 1, temperatures[grid_i2] + deltaT*f1)
			index = grid_i2 + 1

	############################################################################################################

	# Convert data in lists to Numpy arrays and add them to a dictionary that is returned
	data = {}
	data['lats'] = np.array(lats)
	data['lons'] = np.array(lons)
	data['u_winds'] = np.array(u_winds)
	data['v_winds'] = np.array(v_winds)
	data['temperatures'] = np.array(temperatures)
	data['altitudes'] = np.array(altitudes)
	all_pressures = []

	for dat in data['lats']:
		all_pressures.append(100*np.array(pressures)) # Convert hPa to Pa

	data['pressures'] = np.array(all_pressures).transpose()

	if descent_only: # h2 = alts_min[grid_i2], h1 = alt0
		if alts_min[grid_i2] > alt0:
			data['pressures'][grid_i2] = data['pressures'][grid_i2]/np.exp(-((p.g_0*p.M_air)/(p.T0*p.R0))*(alts_min[grid_i2] - alt0))
		else:
			data['pressures'][grid_i2 + 1] = data['pressures'][grid_i2]/np.exp(-((p.g_0*p.M_air)/(p.T0*p.R0))*(alts_min[grid_i2] - alt0))
			
	return data

#################################################################################################################

def read_gfs_single(directory=None, area=None, alt0=None, descent_only=False):
	"""
	Method to read a single GFS file

	Arguments
	=========
  	directory : string
		The path plus name of the grib file to be read (without .grb2 at the end)
	area : (float, float, float, float)
		The (most northern latitude, most western longitude, most southern latitude, most eastern longitude) in degrees, indicating the area on the globe from which the data should be collected
	alt0 : float
		Initial temperature.
	descent_only : bool
		Option to start the trajectory at its highest point. If True, the trajectory only has a descent phase
	"""

	all_data = []

	fname = os.path.join(directory, (directory + '.grb2'))
	main_run_data = read_gfs_file(fname=fname, area=area, alt0=alt0, descent_only=descent_only)
	all_data.append(main_run_data)

	return all_data

#################################################################################################################

def read_gefs_file(fname=None, area=None, alt0=0, descent_only=False):
	"""
	Method to read am ensemble GFS file (GEFS)

	Arguments
	=========
	fname : string
		The name of the grib file to be read
	area : (float, float, float, float)
		The (most northern latitude, most western longitude, most southern latitude, most eastern longitude) in degrees, indicating the area on the globe from which the data should be collected
	alt0 : float
		Initial temperature.
	descent_only : bool
		Option to start the trajectory at its highest point. If True, the trajectory only has a descent phase
	"""

	import param_file as p

	indir = p.path + p.weather_data_folder + p.GEFS_folder

	if area is not None:
		tlat, llon, blat, rlon = area
	else:
		print('Do you really wish to search the entire planet?')
		tlat, llon, blat, rlon = (90.0, 0.0, -90.0, 360.0)

	############################################################################################################

	grib = pg.open(indir + fname)
	grib.seek(0)
	u_msgs = grib.select(name='U component of wind')
	v_msgs = grib.select(name='V component of wind')
	g_msgs = grib.select(name='Geopotential Height')

	lats, lons = u_msgs[0].latlons()
	lats2, lons2 = u_msgs[0].latlons()

	for i in range(len(lons2)):
		for j in range(len(lons2[i])):
			if lons2[i][j] > 180:
				lons2[i][j] -= 360

	locs = pyb_aux.all_and([lats <= tlat, lats >= blat, lons2 <= rlon, lons2 >= llon])
	row_idx, col_idx = np.where(locs)
	lats = lats[(row_idx, col_idx)]
	lons = lons[(row_idx, col_idx)]

	if len(lats) == 0:
		print('Warning! lats is empty!')
	if len(lons) == 0:
		print('Warning! lons is empty!')

	############################################################################################################

	u_wind, v_wind, altitude = {}, {}, {}
	for msg in u_msgs:
		if msg.typeOfLevel == 'isobaricInhPa':
			u_wind[msg.level] = msg.values[(row_idx, col_idx)]

	for msg in v_msgs:
		if msg.typeOfLevel == 'isobaricInhPa':
			v_wind[msg.level] = msg.values[(row_idx, col_idx)]

	for msg in g_msgs:
		if msg.typeOfLevel == 'isobaricInhPa':
			altitude[msg.level] = msg.values[(row_idx, col_idx)]

	############################################################################################################

	pressures = list(u_wind.keys())
	pressures.sort()
	pressures.reverse()

	alt_keys, u_winds, v_winds, altitudes, alt_interp = ([], [], [], [], [])

	for key in pressures:

		uwnd, vwnd, alt = [], [], []

		uwnd.append(u_wind[key])
		vwnd.append(v_wind[key])
		u_winds.append(np.hstack(uwnd))
		v_winds.append(np.hstack(vwnd))

		if key in altitude.keys():

			alt_keys.append(key)
			alt.append(altitude[key])
			altitudes.append(np.hstack(alt))

	############################################################################################################

	p_interp_vals = list(set(pressures).symmetric_difference(set(alt_keys)))

	alt_interps = []

	for p in p_interp_vals:

		alt_interps.append([float(interpolate.interp1d(alt_keys, np.array(altitudes)[:, i])(p)) for i in range(len(lats))])

	for i in range(len(alt_interps)):

		interp_mean = np.mean(alt_interps)
		means = [np.mean(altitudes[k]) for k in range(len(altitudes))]

		index = 0
		for j in range(len(means)):
			index = j
			if interp_mean < means[j]:
				break
				
		altitudes.insert(index , np.array(alt_interps[i]))

	############################################################################################################

	data = {}
	data['lats'] = np.array(lats)
	data['lons'] = np.array(lons)
	data['u_winds'] = np.array(u_winds)
	data['v_winds'] = np.array(v_winds)
	data['altitudes'] = np.array(altitudes)

	all_pressures = []

	for dat in data['lats']:
		all_pressures.append(100 * np.array(pressures))

	data['pressures'] = np.array(all_pressures).transpose()

	return data

#################################################################################################################
#################################################################################################################

def save_kml(fname, data, other_info=None, params=None, balloon=None, shape='circle', err=0, radius=5, mean_direction=None):
	"""
	Method to save given trajectory as KML file

	Arguments
	=========
	fname : string
		Name of file containing the trajectory
	data : dict
		Dictionary containing the trajectory data
	other_info : list
		List containing the highest point in the trajectory, i.e. either the starting location if descent_only or the burst location of the balloon if not descent_only
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	shape : string
		Indicates the shape of the sigma contours drawn around the predicted landing point
	radius : float
		Radius of circle, if circle is chosen contour shape
	mean_direction : float
		Angle between starting point and landing point. Used to orientate ellips if ellips is the chosen contour shape
	"""

	descent_only, next_point, time_interpolate, grid_interpolation, drift_time, resolution, hr_diff, check_sigmas, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	kml_str = '<?xml version="1.0" encoding="UTF-8"?>\n'
	kml_str += '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">\n'
	kml_str += '<Document id="feat_2">\n'
	kml_str += '<name>pyballoon trajectory</name>\n'

	ind = fname.find(')') + 1
	if fname[61] != '/':
		ind0 = 61
	else:
		ind0 = 62

	description = str(fname[ind0:ind])
	description = 'Initial condition:\n' + description.replace('_', ', ')
	kml_str += '<description>' + description + '</description>\n'

	kml_str += '<Style id="stylesel_362">\n'
	kml_str += '<LineStyle id="substyle_363">\n'
	kml_str += '<color>BF0000DF</color>\n'
	kml_str += '<colorMode>normal</colorMode>\n'
	kml_str += '<width>3</width>\n'
	kml_str += '</LineStyle>\n'
	kml_str += '<PolyStyle id="substyle_364">\n'
	kml_str += '<color>BF0000DF</color>\n'
	kml_str += '<colorMode>normal</colorMode>\n'
	kml_str += '<fill>1</fill>\n'
	kml_str += '<outline>1</outline>\n'
	kml_str += '</PolyStyle>\n'
	kml_str += '</Style>\n'
	kml_str += '<Placemark id="feat_91">\n'
	kml_str += '<styleUrl>#stylesel_362</styleUrl>\n'
	kml_str += '<LineString id="geom_86">\n'

	num = 0

	if num == 0:

		kml_str += '<coordinates>\n'

		t_prev = -2.0
		for i in range(0, len(data['lats'])):
			kml_str += '%f,%f,%f\n' % (data['lons'][i], data['lats'][i], data['alts'][i])

		kml_str += '</coordinates>\n'
		kml_str += '<extrude>1</extrude>\n'
		kml_str += '<tessellate>1</tessellate>\n'
		kml_str += '<altitudeMode>absolute</altitudeMode>\n'
		kml_str += '</LineString>\n'
		kml_str += '</Placemark>\n'
		
	num += 1

	# Add placemarks for the trajectory end-points
	num = 0
	kml_str += '<Placemark>\n'
	if num == 0:
		kml_str += '<name>Landing point</name>\n'
		kml_str += '<description>End-point based on GFS main run</description>\n'
	else:
		kml_str += '<name># %d</name>\n' % (num-1)
		kml_str += '<description>End-point based on GFS ensemble member %d</description>\n' % (num-1)

	kml_str += '<Point>\n'
	kml_str += '<altitudeMode>absolute</altitudeMode>\n'
	kml_str += '<coordinates>%f,%f,%f</coordinates>\n' % (data['lons'][-1], data['lats'][-1], data['alts'][-1])
	kml_str += '</Point>\n'
	kml_str += '</Placemark>\n'

	num += 1

	# Add "other_info" places
	if other_info is not None:
		for dat in other_info:
			kml_str += '<Placemark>\n'
			kml_str += '<name>initial condition</name>\n' ## change this
			kml_str += '<description>initial condition</description>\n'
			kml_str += '<Point>\n'
			kml_str += '<altitudeMode>absolute</altitudeMode>\n'
			kml_str += '<coordinates>%f,%f,%f</coordinates>\n' % \
				(dat[1], dat[0], dat[2])
			kml_str += '</Point>\n'
			kml_str += '</Placemark>\n'

	end_lat, end_lon, end_alt = data['lats'][-1], data['lons'][-1], data['alts'][-1]

	if shape == 'ellips':
		sigma_T = p.sigma_T
		sigma_X = p.sigma_X
	if shape == 'circle':
		sigma_r = p.sigma_r

	if shape == 'circle':
		kml_str_add1 = create_circle(lon = end_lon, lat = end_lat, radius=np.sqrt(err**2 + sigma_r**2), color='ff0000ff') # 68th percentile
		kml_str_add2 = create_circle(lon = end_lon, lat = end_lat, radius=np.sqrt(err**2 + (2*sigma_r)**2), color='ff0080ff') # 95th percentile
		kml_str_add3 = create_circle(lon = end_lon, lat = end_lat, radius=np.sqrt(err**2 + (3*sigma_r)**2), color='ffff00ff') # 99th percentile
	elif shape == 'ellips':
		kml_str_add1 = create_ellips(lon = end_lon, lat = end_lat, add=err, times=1, x_extent=sigma_T, y_extent=sigma_X, theta=mean_direction, color='ff0000ff') # 68th percentile
		kml_str_add2 = create_ellips(lon = end_lon, lat = end_lat, add=err, times=2, x_extent=sigma_T, y_extent=sigma_X, theta=mean_direction, color='ff0080ff') # 95th percentile
		kml_str_add3 = create_ellips(lon = end_lon, lat = end_lat, add=err, times=3, x_extent=sigma_T, y_extent=sigma_X, theta=mean_direction, color='ffff00ff') # 99th percentile

	kml_str += kml_str_add1
	kml_str += kml_str_add2
	kml_str += kml_str_add3

	kml_str += '</Document>\n'
	kml_str += '</kml>\n'

	fid = open(fname, 'w')
	fid.write(kml_str)
	fid.close()

	del kml_str, fid

#################################################################################################################

def write_out_polygon(coords=None, color='BF0000DF'):
	"""
	Method to create string with given coordiinates of given shape to be added to a kml file

	Arguments
	=========
	coords : list
		List containing the lats and lons of the shape to be drawn
	color : string
		String indicating the color the shape should be drawn in (Hex Color Code)
	"""

	kml_str = '<Style id="stylesel_362">'
	kml_str += '<LineStyle id="substyle_363">'
	kml_str += '<color>'+color+'</color>'
	kml_str += '<colorMode>normal</colorMode>'
	kml_str += '<width>3</width>'
	kml_str += '</LineStyle>'
	kml_str += '<PolyStyle id="substyle_364">'
	kml_str += '<color>'+color+'</color>'
	kml_str += '<colorMode>normal</colorMode>'
	kml_str += '<fill>1</fill>'
	kml_str += '<outline>1</outline>'
	kml_str += '</PolyStyle>'
	kml_str += '</Style>'
	kml_str += '<Placemark id="feat_91">'
	kml_str += '<styleUrl>#stylesel_362</styleUrl>'
	kml_str += '<LineString id="geom_86">'
	kml_str += '<coordinates>'

	for i in range(len(coords)):
		kml_str += str(coords[i][0]) + ',' + str(coords[i][1]) + ',' +  '0' +'\n'

	kml_str += '</coordinates>\n'
	kml_str += '<extrude>1</extrude>\n'
	kml_str += '<tessellate>1</tessellate>\n'
	kml_str += '<altitudeMode>clampToGround</altitudeMode>\n'
	kml_str += '</LineString>\n'
	kml_str += '</Placemark>\n'

	return kml_str

#################################################################################################################

def merge_kml(datestr=None, run=None, params=None, balloon=None, drift_times=None):
	"""
	Method to merge the kml files output by the pyb_drifter script

	Arguments
	=========
	datestr : string
		Date of initial point
	run : string
		String indicating the folder in which the drifted trajectories are stored
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	drift_times : array
		Array of driftimes used to created the trajectories
	"""

	descent_only, next_point, time_interpolate, grid_interpolation, drift_time, resolution, hr_diff, check_sigmas, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	dir_base = p.path + p.output_folder + str(run) + '/'

	kml_str1 = '<?xml version="1.0" encoding="UTF-8"?>\n'
	kml_str1 += '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">\n'
	kml_str1 += '<Document id="feat_2">\n'

	dir1 = dir_base + p.kml_folder
	dir2 = dir_base + p.traj_folder

	endpoints = {}

	fnames1 = [f for f in os.listdir(dir1) if datestr in f]
	fnames2 = [f for f in os.listdir(dir2) if datestr in f]

	fnames1.sort()
	fnames2.sort()

	for f in range(len(fnames1)):

		info = [line.rstrip('\n').split(' ') for line in open(dir1 + fnames1[f])]
		data_pred = ascii.read(dir2 + fnames2[f])

		end_lat_pred = data_pred['lats'][-1]
		end_lon_pred = data_pred['lons'][-1]
		end_alt_pred = data_pred['alts'][-1]

		endpoints[f] = [end_lat_pred, end_lon_pred, end_alt_pred]

		for i in range(len(info)):

			if info[i] == ['<name>pyballoon', 'trajectory</name>']:
				ind1 = i
			elif info[i] ==  ['</Document>']:
				ind2 = i - 1
			elif info[i] == ['<name>initial', 'condition</name>']:
				ind3 = i - 1
				ind4 = ind3 + 7

		for i in range(len(info)):
			if i >= ind1 and i <= ind2:
				if f != 0:
					if i < ind3 or i > ind4:
						for j in range(len(info[i])):
							if j == len(info[i]) - 1:
								kml_str1 += str(info[i][j]) + '\n'
							else:
								kml_str1 += str(info[i][j]) + ' '
				else:
					for j in range(len(info[i])):
						if j == len(info[i]) - 1:
							kml_str1 += str(info[i][j]) + '\n'
						else:
							kml_str1 += str(info[i][j]) + ' '
									
	kml_str1 += '</Document>\n'
	kml_str1 += '</kml>\n'

	kml_fname1 = 'merged_' + datestr + '.kml'

	output_dir = dir1

	fid = open(output_dir + kml_fname1, 'w')
	fid.write(kml_str1)
	fid.close()

	kml_fname2 = 'endpoints_merged_' + datestr + ' .kml'

	kml_str2 = '<?xml version="1.0" encoding="UTF-8"?>\n'
	kml_str2 += '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">\n'
	kml_str2 += '<Document id="feat_2">\n'
	kml_str2 += '<name>pyballoon trajectory</name>\n'
	kml_str2 += '<Style id="stylesel_362">\n'
	kml_str2 += '<LineStyle id="substyle_363">\n'
	kml_str2 += '<color>BF0000DF</color>\n'
	kml_str2 += '<colorMode>normal</colorMode>\n'
	kml_str2 += '<width>3</width>\n'
	kml_str2 += '</LineStyle>\n'
	kml_str2 += '<PolyStyle id="substyle_364">\n'
	kml_str2 += '<color>BF0000DF</color>\n'
	kml_str2 += '<colorMode>normal</colorMode>\n'
	kml_str2 += '<fill>1</fill>\n'
	kml_str2 += '<outline>1</outline>\n'
	kml_str2 += '</PolyStyle>\n'
	kml_str2 += '</Style>\n'
	kml_str2 += '<Placemark id="feat_91">\n'
	kml_str2 += '<styleUrl>#stylesel_362</styleUrl>\n'
	kml_str2 += '<LineString id="geom_86">\n'
	kml_str2 += '<coordinates>\n'

	keys = list(endpoints.keys())
	keys.sort()

	for f in range(len(fnames1)):
		kml_str2 += str(endpoints[f][1]) + ',' + str(endpoints[f][0])  + ',' + str(endpoints[f][2]) + '\n'

	kml_str2 += '</coordinates>\n'
	kml_str2 += '<extrude>1</extrude>\n'
	kml_str2 += '<tessellate>1</tessellate>\n'
	kml_str2 += '<altitudeMode>clampToGround</altitudeMode>\n'
	kml_str2 += '</LineString>\n'
	kml_str2 += '</Placemark>\n'

	for f in range(len(fnames1)):
		kml_str2 += '<Placemark>\n'
		kml_str2 += '<name>End, drift: '+ str(drift_times[f]) + ' min.</name>\n'
		kml_str2 += '<Point>\n'

		if endpoints[f][1] >= 180.:
			endpoints[f][1] -= 360

		kml_str2 += '<coordinates>' + str(endpoints[f][1]) + ',' + str(endpoints[f][0])  + ',' + str(endpoints[f][2]) + '</coordinates>' + '\n'
		kml_str2 += '</Point>\n'
		kml_str2 += '</Placemark>\n'

	kml_str2 += '</Document>\n'
	kml_str2 += '</kml>'

	fid = open(output_dir + kml_fname2, 'w')
	fid.write(kml_str2)
	fid.close()

#################################################################################################################
#################################################################################################################

def print_verbose(datestr=None, utc_hour=None, loc0=None, params=None, balloon=None):
	"""
	Method to print out the parameters being used to the terminal

	Arguments
	=========	
	datestr : string
		Date of initial point
	utc_hour : float
		Time of initial point
	loc0 : floats in tuple
	  	(latitude in degrees, longitude in degrees, altitude in km) of initial point
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	"""

	descent_only, next_point, time_interpolate, grid_interpolate, drift_time, resolution, hr_diff, check_sigmas, params, balloon = set_params(params=params, balloon=balloon)

	print('General Parameters')
	print('----------')
	print('descent_only: ' + str(descent_only))
	if descent_only:
		print('starting point: ' + next_point)
	print('time_interpolate: ' + str(time_interpolate))
	print('grid_interpolate: ' + str(grid_interpolate))
	if drift_time is not None:
		print('drift time: ' + str(drift_time) + ' minutes')
	print('resolution of forecasts: ' + str(resolution) + ' degrees')
	if hr_diff is not None:
		print('difference in hrs for forecasts: ' + str(hr_diff) + ' hours')
	print('check sigmas: ' + str(check_sigmas))
	print('\nBalloon/Parachute Parameters')
	print('----------')
	print('altitude step: ' + str(balloon['altitude_step']) + ' m')
	print('equipment mass: ' + str(balloon['equip_mass']) + ' kg')
	print('parachute Cd: ' + str(balloon['Cd_parachute']))
	print('parachute area: ' + str(round(balloon['parachute_areas'][0], 2)) + ' m^2')

	if datestr != None and utc_hour !=  None:
		print('----------\n')
		print('Running date/time: ' + datestr + ', ' + str(utc_hour) + ' hr')
	if loc0 != None:
		print('Starting point: ' + str(loc0[0]) + ' lat., ' + str(loc0[1]) + ' lon., ' + str(loc0[2]) + ' m\n')
	
	print('----------\n')

#################################################################################################################

def write_verbose(params_dir=None, params=None, balloon=None):
	"""
	Method to write the parameters used of run to file

	Arguments
	=========	
	params_dir : string
		Directory in which to save the params.txt file
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	"""

	descent_only, next_point, time_interpolate, grid_interpolate, drift_time, resolution, hr_diff, check_sigmas, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	f = open(params_dir + 'params.txt', 'w+')
	f.write('General parameters\n')
	f.write('----------------------\n')
	f.write('descent_only: ' + str(descent_only) + '\n')
	if descent_only:
		f.write('starting point: ' + str(next_point) + '\n')
	f.write('time interpolation: ' + str(time_interpolate) + '\n')
	f.write('grid_interpolation: ' + str(grid_interpolate) + '\n')
	f.write('drift time: ' + str(drift_time) + ' min\n')
	f.write('resolution of forecasts: ' + str(resolution) + ' degrees\n')
	f.write('difference in hrs for forecasts: ' + str(hr_diff) + ' hours\n')
	f.write('check sigmas: ' + str(check_sigmas) + '\n')
	f.write('----------------------\n')
	f.write('\n')
	f.write('Balloon/Parachute parameters\n')
	f.write('----------------------\n')
	f.write('altitude step: ' + str(balloon['altitude_step']) + ' m\n')
	f.write('equipment mass: ' + str(balloon['equip_mass']) + ' kg\n')
	f.write('parachute Cd: ' + str(balloon['Cd_parachute']) + '\n')
	f.write('parachute area: ' + str(round(balloon['parachute_areas'][0], 2)) + ' m^2\n')

	f.close()

#################################################################################################################

def set_params(params=None, balloon=None):
	"""
	Method to initialise the parameters being used

	Arguments
	=========	
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	"""

	# balloon parameters
	if balloon == None:
		balloon = p.balloon

	next_point = '0'

	if params == None:

		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		time_interpolate = p.time_interpolate
		grid_interpolate = p.grid_interpolate
		drift_time = float(p.drift_time)
		resolution = float(p.resolution)
		hr_diff = int(p.hr_diff)
		check_sigmas = bool(p.check_sigmas)

		params = [descent_only, next_point, time_interpolate, grid_interpolate, drift_time, resolution, hr_diff, check_sigmas]
		
	else:

		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		time_interpolate = bool(params[-6])
		grid_interpolate = bool(params[-5])

		if type(params[-4]) != type(np.array([])):
			drift_time = float(params[-4])
		else:
			drift_time = params[-4]

		resolution = float(params[-3])

		if type(params[-2]) != type(np.array([])):
			hr_diff = int(params[-2])
		else:
			hr_diff = params[-2]

		check_sigmas = bool(params[-1])

	return descent_only, next_point, time_interpolate, grid_interpolate, drift_time, resolution, hr_diff, check_sigmas, params, balloon

#################################################################################################################

def write_run_info(add_run_info=True, run=None, params=None, balloon=None):
	"""
	Method to write the parameters being used to the runs_info.txt file which contains the info for all runs run

	Arguments
	=========	
	params_dir : string
		Directory in which to save the params.txt file
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	"""

	descent_only, next_point, time_interpolate, grid_interpolate, drift_time, resolution, hr_diff, check_sigmas = params

	if add_run_info:

		run_info_file = p.path + p.output_folder + 'runs_info.txt'

		if not os.path.isfile(run_info_file):

			f = open(run_info_file, 'w+')
			f.write('run descent_only next_point time_interpolate grid_interpolate drift_time resolution hr_diff check_sigmas Cd_parachute\
			 parachute_area altitude_step equip_mass balloon_mass fill_radius radius_empty burst_radius thickness_empty Cd_balloon simple_ascent_rate parachute_change_altitude')

		lines = [line.rstrip('\n').split(' ') for line in open(run_info_file)]
		runs = [lines[i][0] for i in range(len(lines)) if i != 0]

		f = open(run_info_file, 'a+')
		if run not in runs:
			f.write('\n' + str(run) + ' ' + str(descent_only) + ' ' + str(next_point) + ' ' + str(time_interpolate) + ' ' + str(grid_interpolate) + ' ' + str(drift_time) + ' ' + \
				str(resolution) + ' '  + str(hr_diff) + ' ' + str(check_sigmas)  + ' ' + str(balloon['Cd_parachute']) + ' ' + str(balloon['parachute_areas'][0]) \
				+ ' ' + str(balloon['altitude_step'])  + ' ' + str(balloon['equip_mass']) + ' ' + str(balloon['balloon_mass']) + ' ' + str(balloon['fill_radius']) \
				 + ' ' + str(balloon['radius_empty']) + ' ' + str(balloon['burst_radius']) + ' ' + str(balloon['thickness_empty']) + ' ' + str(balloon['Cd_balloon']) \
				 + ' ' + str(balloon['simple_ascent_rate']) + ' ' + str(balloon['parachute_change_altitude']))
		f.close()

#################################################################################################################

def search_info(run=None, print_verbose=True):
	"""
	Method to find the info for a run in the runs_info.txt file

	Arguments
	=========	
	run : string
		String indicating which run we wish to know the parameters of
	print_verbose : bool
		If True, the parameters will be printed to the command line
	"""

	data = ascii.read(p.path + p.output_folder + 'runs_info.txt')
	runs = data['run']

	index = np.where(runs == run)[0]

	if len(index) < 1:

		print('Cannot find the info for this run! Check the folder.')
		sys.exit()

	else:

		index = int(index)

		if print_verbose:

			print('\nGeneral Parameters')
			print('----------')
			print('descent_only: ' + str(data['descent_only'][index]))
			if data['descent_only'][index]:
				print('starting point: ' + str(data['next_point'][index]))
			print('time_interpolate: ' + str(data['time_interpolate'][index]))
			print('grid_interpolate: ' + str(data['grid_interpolate'][index]))
			print('drift time: ' + str(data['drift_time'][index]) + ' minutes')
			print('resolution of forecasts: ' + str(data['resolution'][index]) + ' degrees')
			print('difference in hrs for forecasts: ' + str(data['hr_diff'][index]) + ' hours')
			print('check sigmas: ' + str(data['check_sigmas'][index]))
			print('----------')
			print('\nBalloon/Parachute Parameters')
			print('----------')
			print('altitude step: ' + str(data['altitude_step'][index]) + ' m')
			print('equipment mass: ' + str(data['equip_mass'][index]) + ' kg')
			print('parachute Cd: ' + str(data['Cd_parachute'][index]))
			print('parachute area: ' + str(data['parachute_area'][index]) + ' m^2')
			print('----------\n')

		return index

#################################################################################################################

def save_used_weather_files(save_dir=None, used_weather_files=None, trajectories=None, utc_hour=None):
	"""
	Method to write out the weather files used for a run

	Arguments
	=========	
	save_dir : string
		Directory in which to store file with used weather files
	used_weather_files : dict
		Dictionary of used weather files, with the times at which the files were first used as keys
	trajectories : dict
		Dictionary containing trajectory data
	utc_hour : float
		Initial time of trajectory
	"""

	if not os.path.exists(save_dir + 'used_weather_files.txt'):
		f = open(save_dir + 'used_weather_files.txt', 'w')
		f.write('time file \n------------------------------------------\n')
		f.close()

	keys = list(used_weather_files.keys())
	keys.sort()

	f = open(save_dir + 'used_weather_files.txt', 'a+')
	for key in keys:
		if type(used_weather_files[key]) == type(list()):
			for file in range(len(used_weather_files[key])):
				f.write(str(key) + ' ' + str(used_weather_files[key][file]) + '\n')
		else:
			f.write(str(key) + ' ' + str(used_weather_files[key]) + '\n')
	f.write(str((utc_hour + trajectories['times'][-1]/60.) % 24) + ' -\n')
	f.write('------------------------------------------\n')
	f.close()

#################################################################################################################

def determine_error(trajectories=None, params=None):
	"""
	Method to determine the error in the trajectory from the difference in location, time to forecast and how far into the future the forecast looks

	Arguments
	=========
	trajectories : dict
		Dictionary containing trajectory data
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	"""

	if params == None:
		time_interpolate = bool(p.time_interpolate)
	else:
		time_interpolate = params[-7]

	mean_tfut = np.mean(trajectories['tfutures'])
	err = np.sqrt((mean_tfut*(0.19))**2) ## for forecast looking into the future

	if not time_interpolate:
		time_diff_err = np.mean(np.abs(trajectories['delta_ts']))*2.4 + -0.3 # for not using interpolation: i.e., difference between forecast time and model time
		err = np.sqrt(time_diff_err**2 + err**2)

	return err

#################################################################################################################

def create_trajectory_files(traj_dir=None, kml_dir=None, datestr=None, utc_hour=None, loc0=None, trajectories=None, params=None, overwrite=False):
	"""
	Method to create trajectory files and write out the trajectory data to them

	Arguments
	=========
	traj_dir : string
		Directory where the trajectory files will be saved
	kml_dir : string 
		Directory where the kml files will be saved
	datestr : string
		Date of initial point
	utc_hour : float
		Initial time of trajectory
	loc0 : floats in tuple
		(latitude in degrees, longitude in degrees, altitude in km) of initial point
	trajectories : dict
		Dictionary containing trajectory data
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	"""

	if params == None:
		time_interpolate = bool(p.time_interpolate)
		check_sigmas = bool(p.check_sigmas)
	else:
		time_interpolate = params[-6]
		check_sigmas = bool(params[-1])

	traj_file = traj_dir + 'trajectory_' + datestr + '_' + str(utc_hour) + '_' + str(loc0)
	kml_fname = kml_dir + datestr + '_' + str(utc_hour) + '_' + str(loc0)

	# determine number of files already in folder (with the same starting params)
	if overwrite:
		no = 0
	else:
		no = len([filename for filename in os.listdir(traj_dir) if os.path.basename(traj_file) in filename])

	if no == 0:
		traj_file += '.dat'
		kml_fname += '.kml'
	else:
		traj_file += '_' + str(no + 1) + '.dat'
		kml_fname += '_' + str(no + 1) + '.kml'

	# write out trajectory file
	out_list = [trajectories['lats'], trajectories['lons'], trajectories['alts'], trajectories['dists'], trajectories['times'], trajectories['tfutures'], trajectories['loc_diffs'], trajectories['speeds']]
	names = ['lats', 'lons', 'alts', 'dists', 'times', 'tfutures', 'loc_diffs', 'speeds']

	if check_sigmas:
		out_list += [trajectories['sigmas_u'], trajectories['sigmas_v']]
		names += ['sigmas_u', 'sigmas_v']

	if not time_interpolate:
		out_list = out_list[:5] + [trajectories['delta_ts']] + out_list[5:]
		names = names[:5] + ['delta_t'] + names[5:]

	ascii.write(out_list, traj_file, names=names, overwrite=True)

	return kml_fname	

#################################################################################################################

def make_descent_rate_plot(directory=None, data=None, datestr=None, utc_hour=None, loc0=None):
	"""
	Method to create plots showing the descent rates of the trajectory

	Arguments
	=========
	directory : string
		Directory where the figures will be saved
	data : dict
		Dictionary containing the trajectory data
	datestr : string
		Date of initial point
	utc_hour : float
		Time of initial point
	loc0 : floats in tuple
		(latitude in degrees, longitude in degrees, altitude in km) of initial point		
	"""

	fig_dir = directory + 'DescentRates/'

	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	alts = np.array(data['alts'])
	times = np.array(data['times'])
	descent_speeds = np.array(data['speeds'])

	fig = plt.figure()

	plt.plot(times, descent_speeds, 'o', markersize=1)

	plt.xlabel('Time since drop [min]', fontsize=15)
	plt.ylabel('Descent speed [m/s]', fontsize=15)

	plt.grid(True)
	plt.tight_layout()

	fig.savefig(fig_dir + 'descent_rate_vs_time.png')

	plt.clf()

	plt.plot(alts, descent_speeds, 'o', markersize=1)

	plt.xlabel('Altitude [m]', fontsize=15)
	plt.ylabel('Descent speed [m/s]', fontsize=15)

	plt.grid(True)
	plt.tight_layout()

	fig.savefig(fig_dir + 'descent_rate_vs_alt_' + datestr + '_' + str(utc_hour) + '_' + str(loc0) + '.png')

	plt.close()

#################################################################################################################
#################################################################################################################

def geodesic_point_buffer(lat=None, lon=None, radius=5):
	"""
	Method to find all latslons to create a circle of a given radius (km) around a lat/lon point

	lon : float
		Longitude in degrees of centre of circle (i.e. lon of landing point)
	lat : float
		Latitude in degrees of centre of circle (i.e. lat of landing point)
	radius : float
		Radius in km of circle
	"""

	proj_wgs84 = pyproj.Proj(init='epsg:4326')

	# Azimuthal equidistant projection
	aeqd_proj = '+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0'
	project = partial(
		pyproj.transform,
		pyproj.Proj(aeqd_proj.format(lat=lat, lon=lon)),
		proj_wgs84)
	buf = Point(0, 0).buffer(radius * 1000)  # distance in metres
	return transform(project, buf).exterior.coords[:]

#################################################################################################################

def make_ellipse(theta_num=100, phi=0, x_cent=0, y_cent=0, semimaj=1.9, semimin=1.4):
	"""
	Method to create an ellips and get x/y coordinates given ellips parameters

	Arguments
	=========
	theta_num : int
		Number of coordinates in the ellips
	phi : float
		Angle in radians indicating the rotation of the ellips
	x_cent : float
		X coordinate of the centre of the ellips
	y_cent : float
		Y coordinate of the centre of the ellips
	semimaj : float
		Semimajor-axis of ellips
	semimin : float
		Semiminor-axis of ellips
	"""

	theta = np.linspace(0, 2*np.pi, theta_num)
	r = 1/np.sqrt((np.cos(theta))**2 + (np.sin(theta))**2)
	x = r*np.cos(theta)
	y = r*np.sin(theta)

	data = np.array([x, y])

	S = np.array([[semimaj, 0], [0, semimin]])
	R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
	T = np.dot(R, S)

	data = np.dot(T, data)

	data[0] += x_cent
	data[1] += y_cent

	return data

#################################################################################################################

def get_ellips_coords(lat=None, lon=None, x_extent=1.87, y_extent=1.47, theta=0):
	"""
	#Method to find all lats/lons to create a ellips of a given minor/major axis (km) around a lat/lon point

	Arguments
	=========
	lon : float
		Longitude in degrees of centre of ellips (i.e. lon of landing point)
	lat : float
		Latitude in degrees of centre of ellips (i.e. lat of landing point)
	x_extent : float
		Extent in km in the x (longitude) direction, taken from percentile calculations
	y_extent : float
		Extent in km in the y (latitude) direction, taken from percentile calculations
	theta : float
		Angle in radians indicating the rotation of the ellipse. Calculated from the mean direction of the trajectory
	"""

	minor_axis = min([x_extent, y_extent])
	major_axis = max([x_extent, y_extent])

	if major_axis == y_extent:
		theta += np.radians(90.)

	coords = np.array(geodesic_point_buffer(lat, lon, 0.5))
	lons = coords[:, 0]
	lats = coords[:, 1]

	min_lon, max_lon = min(lons), max(lons)
	min_lat, max_lat = min(lats), max(lats)

	lon_deg_per_km =  max_lon - min_lon
	lat_deg_per_km =  max_lat - min_lat

	cart_coords = make_ellipse(theta_num=100, phi=theta, x_cent=0, y_cent=0, semimaj=major_axis, semimin=minor_axis)

	x_coords = cart_coords[0]
	y_coords = cart_coords[1]

	new_coords = []
	for i in range(len(x_coords)):

		new_lon = lon + lon_deg_per_km*x_coords[i]
		new_lat = lat + lat_deg_per_km*y_coords[i]

		new_coords.append((new_lon, new_lat))

	return new_coords

#################################################################################################################

def create_circle(lon=None, lat=None, radius=5, color='BF0000DF'):
	"""
	Method to create a circle of given radius around the landing point in the kml files.

	Arguments
	=========
	lon : float
		Longitude in degrees of centre of circle (i.e. lon of landing point)
	lat : float
		Latitude in degrees of centre of circle (i.e. lat of landing point)
	radius : float
		Radius in km of circle
	color : string
		Color of the circle to be drawn (Hex Color Code)
	"""

	coords = np.array(geodesic_point_buffer(lat, lon, radius))
	kml_str = write_out_polygon(coords=coords, color=color)
	return kml_str

#################################################################################################################

def create_ellips(lon=None, lat=None, x_extent=1.858, y_extent=1.474, theta=0, add=0, times=1, color='BF0000DF'):
	"""
	Method to create an ellips around the landing point in the kml files.

	Arguments
	=========
	lon : float
		Longitude in degrees of centre of ellips (i.e. lon of landing point)
	lat : float
		Latitude in degrees of centre of ellips (i.e. lat of landing point)
	x_extent : float
		Extent in km in the x (longitude) direction, taken from percentile calculations
	y_extent : float
		Extent in km in the y (latitude) direction, taken from percentile calculations
	theta : float
		Angle in radians indicating the rotation of the ellipse. Calculated from the mean direction of the trajectory
	add : float
		Error from trajectory to be added to the contours
	times : int
		Times indicates which sigma limit the contour represents (i.e. 1, 2, 3)
	color : string
		Color of the ellips to be drawn (Hex Color Code)
	"""

	x_extent *= times
	y_extent *= times
	x_extent = np.sqrt(add**2 + x_extent**2)
	y_extent = np.sqrt(add**2 + y_extent**2)

	coords = np.array(get_ellips_coords(lat=lat, lon=lon, x_extent=x_extent, y_extent=y_extent, theta=theta))
	kml_str = write_out_polygon(coords=coords, color=color)

	return kml_str

#################################################################################################################

if __name__ == '__main__':

	read_gfs_file('/home/ellen/Desktop/SuperBIT_DRS/Weather_data/GFS/gfs_4_20190914_1200_003.grb2', alt0=26091, descent_only=True)

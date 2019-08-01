"""Input and output functions used by pyBalloon"""

from shapely.geometry import Point
from shapely.ops import transform
from functools import partial
from scipy import interpolate
from astropy.io import ascii
import pygrib as pg
import numpy as np
import sys, os
import pyproj

import pyb_aux
import pyb_io

import param_file as p

#################################################################################################################

# method to collect relevant data from GFS file
def read_gfs_file(fname, area=None, alt0=0, t_0=None, extra_data=None, descent_only=False, step=100):

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
	omega_msgs = grib.select(name='Vertical velocity')

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

	if len(lats) == 0: print( 'Warning! lats is empty!')
	if len(lons) == 0: print( 'Warning! lons is empty!')

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

	omega = {}
	for msg in omega_msgs:
		if msg.typeOfLevel == 'isobaricInhPa':
			omega[msg.level] = msg.values[row_idx, col_idx]

	############################################################################################################

	pressures = list(u_wind.keys())

	# Collect data to correct altitude order. Set "surface" values before real data.
	# Use given surface temperature if available, otherwise use the model value
	if not descent_only:

		u_winds = [np.zeros(lats.shape)]
		v_winds = [np.zeros(lats.shape)]
		altitudes = [alt0*np.ones(lats.shape)]
		omegas = [np.zeros(lats.shape)]

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
		omegas = []

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

			if key in list(omega.keys()):
				omg.append(omega[key])

			# Add extra data to complement the currently read data, ie. 1, 2, 3, 5 and 7 hPa levels from GFS main run to ensembles. 
			# Data are expected to be in the same format as returned by this function.

			j = 0
			if extra_data is not None:

				for level in extra_data['pressures'][:, 0]/100:

					if level < min(pressures):

						for idx in range(0, len(extra_data['u_winds'][0, :])):

							uwnd.append(extra_data['u_winds'][j, idx])
							vwnd.append(extra_data['v_winds'][j, idx])
							temp.append(extra_data['temperatures'][j, idx])
							alt.append(extra_data['altitudes'][j, idx])
							omg.append(extra_data['omegas'][j, idx])
					j += 1

			u_winds.append(np.hstack(uwnd))
			v_winds.append(np.hstack(vwnd))
			temperatures.append(np.hstack(temp))
			altitudes.append(np.hstack(alt))

			if key in list(omega.keys()):
				omegas.append(np.hstack(omg))

			i+=1
		else:
			i+=1

	############################################################################################################

	if descent_only:
		main_keys = list(u_wind.keys())
	else:
		main_keys = list(pressures)

	main_keys.sort()
	main_keys.reverse()
	omega_keys = list(omega.keys())

	omegas = np.array(omegas)
	x, y = omegas.shape

	omega_ext = {}

	for key in main_keys:
		if key not in omega_keys:
			omega_ext[key] = np.zeros(y)
		else:
			omega_ext[key] = omega[key]

	omegas = [omega_ext[key] for key in main_keys]

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
			deltaOmg = omegas[grid_i2] - omegas[grid_i2 - 1]

			altitudes.insert(grid_i2, alt0*np.ones(lats.shape))
			u_winds.insert(grid_i2, u_winds[grid_i2 - 1] + deltaU*f1)
			v_winds.insert(grid_i2, v_winds[grid_i2 - 1] + deltaV*f1)
			temperatures.insert(grid_i2, temperatures[grid_i2 - 1] + deltaT*f1)
			omegas.insert(grid_i2, omegas[grid_i2 - 1] + deltaOmg*f1)
			index = grid_i2

		else:

			deltaA1 = altitudes[grid_i2 + 1] - altitudes[grid_i2]
			deltaA2 = alt0 - altitudes[grid_i2]

			f1 = deltaA2/deltaA1

			deltaT = temperatures[grid_i2 + 1] - temperatures[grid_i2]
			deltaU = u_winds[grid_i2 + 1] - u_winds[grid_i2]
			deltaV = v_winds[grid_i2 + 1] - v_winds[grid_i2]
			deltaOmg = omegas[grid_i2 + 1] - omegas[grid_i2]

			altitudes.insert(grid_i2 + 1, alt0*np.ones(lats.shape))
			u_winds.insert(grid_i2 + 1, u_winds[grid_i2] + deltaU*f1)
			v_winds.insert(grid_i2 + 1, v_winds[grid_i2] + deltaV*f1)
			temperatures.insert(grid_i2 + 1, temperatures[grid_i2] + deltaT*f1)
			omegas.insert(grid_i2 + 1, omegas[grid_i2] + deltaOmg*f1)
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
	data['omegas'] = np.array(omegas)
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

# method to read a single GFS file
def read_gfs_single(directory=None, area=None, alt0=None, descent_only=False, step=100.):

	all_data = []

	fname = os.path.join(directory, (directory + '.grb2'))
	# print('Reading GFS data from ' + fname[-28:])
	main_run_data = read_gfs_file(fname, area=area, alt0=alt0, descent_only=descent_only, step=step)
	all_data.append(main_run_data)

	return all_data

#################################################################################################################

# method to read a GEFS file
def read_gefs_file(fname=None, area=None, alt0=0, t_0=None, extra_data=None, descent_only=False, step=100):

	indir = p.path + p.weather_data_folder + p.GFS_folder

	if area is not None:
		tlat, llon, blat, rlon = area
	else:
		print 'Do you really wish to search the entire planet?'
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
		print 'Warning! lats is empty!'
	if len(lons) == 0:
		print 'Warning! lons is empty!'

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

# method to save given trajectories as KML files.
def save_kml(fname, data, model_start_idx=0, eps_mode='end', other_info=None, params=None, radius=5):

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	kml_str = '<?xml version="1.0" encoding="UTF-8"?>\n'
	kml_str += '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">\n'
	kml_str += '<Document id="feat_2">\n'
	kml_str += '<name>pyballoon trajectory</name>\n'
	kml_str += '<description>' + str(fname[52:69]) + '</description>\n'

	kml_str += '<Style id="stylesel_362">\n'
	kml_str += '<LineStyle id="substyle_363">\n'
	kml_str += '<color>BF0000DF</color>\n'
	kml_str += '<colorMode>normal</colorMode>\n'
	kml_str += '<width>5</width>\n'
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

	if num == 0 or eps_mode == 'full':

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
	kml_str += '<coordinates>%f,%f</coordinates>\n' % (data['lons'][-1], data['lats'][-1])
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

 	end_lat, end_lon = data['lats'][-1], data['lons'][-1]
	kml_str_add = create_circle(lon = end_lon, lat = end_lat, radius=radius)
	kml_str += kml_str_add

	kml_str += '</Document>\n'
	kml_str += '</kml>\n'

	fid = open(fname, 'w')
	fid.write(kml_str)
	fid.close()

#################################################################################################################

# method to find all latslons to create a circle of a given radius around a lat/lon point
def geodesic_point_buffer(lat, lon, km):

	proj_wgs84 = pyproj.Proj(init='epsg:4326')

	# Azimuthal equidistant projection
	aeqd_proj = '+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0'
	project = partial(
		pyproj.transform,
		pyproj.Proj(aeqd_proj.format(lat=lat, lon=lon)),
		proj_wgs84)
	buf = Point(0, 0).buffer(km * 1000)  # distance in metres
	return transform(project, buf).exterior.coords[:]

#################################################################################################################

# method to create a circle of given radius around the landing point in the kml files.
def create_circle(lon=None, lat=None, radius=5):

	coords = np.array(geodesic_point_buffer(lat, lon, radius))

	kml_str = '<name>error circle</name>'
	kml_str += '<description>error circle</description>'
	kml_str += '<Style id="stylesel_362">'
	kml_str += '<LineStyle id="substyle_363">'
	kml_str += '<color>BF0000DF</color>'
	kml_str += '<colorMode>normal</colorMode>'
	kml_str += '<width>5</width>'
	kml_str += '</LineStyle>'
	kml_str += '<PolyStyle id="substyle_364">'
	kml_str += '<color>BF0000DF</color>'
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
	kml_str += '<altitudeMode>relativeToGround</altitudeMode>\n'
	kml_str += '</LineString>\n'
	kml_str += '</Placemark>\n'

	return kml_str

#################################################################################################################

# method to merge the kml files output by the pyb_drifter script
def merge_kml(datestr=None, run=None, params=None, balloon=None, drift_times=None):

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	dir_base = p.path + p.output_folder + str(run) + '/'

	kml_str1 = '<?xml version="1.0" encoding="UTF-8"?>\n'
	kml_str1 += '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">\n'
	kml_str1 += '<Document id="feat_2">\n'

	dir1 = dir_base + p.kml_folder
	dir2 = dir_base + p.traj_folder

	endpoints = {}

	fnames1 = [f for f in os.listdir(dir1) if datestr in f]
	fnames2 = [f for f in os.listdir(dir2) if datestr in f]

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
	kml_str2 += '<width>5</width>\n'
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
	kml_str2 += '<altitudeMode>relativeToGround</altitudeMode>\n'
	kml_str2 += '</LineString>\n'
	kml_str2 += '</Placemark>\n'

	for f in range(len(fnames1)):
		kml_str2 += '<Placemark>\n'
		kml_str2 += '<name>End, drift = '+ str(drift_times[f]) + ' min.</name>\n'
		kml_str2 += '<Point>\n'
		kml_str2 += '<coordinates>' + str(endpoints[f][1]) + ',' + str(endpoints[f][0]) + '</coordinates>' + '\n'
		kml_str2 += '</Point>\n'
		kml_str2 += '</Placemark>\n'

		kml_str_add = create_circle(lon = endpoints[f][1], lat = endpoints[f][0])
		kml_str2 += kml_str_add

	kml_str2 += '</Document>\n'
	kml_str2 += '</kml>'

	fid = open(output_dir + kml_fname2, 'w')
	fid.write(kml_str2)
	fid.close()

#################################################################################################################

# method to print out the parameters being used to the terminal
def print_verbose(datestr=None, utc_hour=None, loc0=None, params=None, balloon=None):

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff = params

	print('General Parameters')
	print('----------')
	print('descent_only: ' + str(descent_only))
	if descent_only:
		print('starting point: ' + next_point)
	print('interpolate: ' + str(interpolate))
	if drift_time is not None:
		print('drift time: ' + str(drift_time) + ' minutes')
	print('resolution of forecasts: ' + str(resolution) + ' degrees')
	print('correct for vertical winds: ' + str(vz_correct))
	if hr_diff is not None:
		print('difference in hrs for forecasts: ' + str(hr_diff) + ' hours')
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
	
	print('----------')

#################################################################################################################

# method to write the parameters used of run to file
def write_verbose(params_dir, params, balloon):

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff = params

	f = open(params_dir + 'params.txt', 'w+')
	f.write('General parameters\n')
	f.write('----------------------\n')
	f.write('descent_only: ' + str(descent_only) + '\n')
	if descent_only:
		f.write('starting point: ' + str(next_point) + '\n')
	f.write('interpolate: ' + str(interpolate) + '\n')
	f.write('drift time: ' + str(drift_time) + ' min\n')
	f.write('resolution of forecasts: ' + str(resolution) + ' degrees\n')
	f.write('correct for vertical winds: ' + str(vz_correct) + '\n')
	f.write('difference in hrs for forecasts: ' + str(hr_diff) + ' hours\n')
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

# method to set the parameters being used
def set_params(params=None, balloon=None):

	# balloon parameters
	if balloon == None:
		balloon = p.balloon

	next_point = '0'

	if params == None:

		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		interpolate = p.interpolate
		drift_time = float(p.drift_time)
		resolution = float(p.resolution)
		vz_correct = bool(p.vz_correct)
		hr_diff = int(p.hr_diff)

		params = [descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff]
		
	else:

		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		interpolate = bool(params[2])
		drift_time = float(params[3])
		resolution = float(params[4])
		vz_correct = bool(params[5])
		hr_diff = int(params[6])

	return descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, params, balloon

#################################################################################################################

# method to write the parameters being used to the runs_info.txt file which contains the info for all runs run
def write_run_info(add_run_info=True, run=None, params=None, balloon=None):

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff = params

	if add_run_info:

		run_info_file = p.path + p.output_folder + 'runs_info.txt'

		if not os.path.isfile(run_info_file):

			f = open(run_info_file, 'w+')
			f.write('run descent_only next_point interpolate drift_time resolution vz_correct hr_diff Cd_parachute parachute_area altitude_step equip_mass balloon_mass fill_radius radius_empty burst_radius\
			 thickness_empty Cd_balloon simple_ascent_rate parachute_change_altitude')

		lines = [line.rstrip('\n').split(' ') for line in open(run_info_file)]
		runs = [lines[i][0] for i in range(len(lines)) if i != 0]

		f = open(run_info_file, 'a+')
		if run not in runs:
			f.write('\n' + str(run) + ' ' + str(descent_only) + ' ' + str(next_point) + ' ' +  str(interpolate) + ' ' + str(drift_time) + ' ' + str(resolution) + ' '  + str(vz_correct) \
				+ ' ' + str(hr_diff)  + ' ' + str(balloon['Cd_parachute']) + ' ' + str(balloon['parachute_areas'][0]) + ' ' + str(balloon['altitude_step'])  + ' ' + str(balloon['equip_mass']) \
				+ ' ' + str(balloon['balloon_mass']) + ' ' + str(balloon['fill_radius']) + ' ' + str(balloon['radius_empty']) + ' ' + str(balloon['burst_radius']) + ' ' + str(balloon['thickness_empty']) \
				+ ' ' + str(balloon['Cd_balloon']) + ' ' + str(balloon['simple_ascent_rate']) + ' ' + str(balloon['parachute_change_altitude']))
		f.close()

#################################################################################################################

# method to find the info for a run in the runs_info.txt file
def search_info(run=None, print_verbose=True):

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
			print('interpolate: ' + str(data['interpolate'][index]))
			print('drift time: ' + str(data['drift_time'][index]) + ' minutes')
			print('resolution of forecasts: ' + str(data['resolution'][index]) + ' degrees')
			print('correct for vertical winds: ' + str(data['vz_correct'][index]))
			print('difference in hrs for forecasts: ' + str(data['hr_diff'][index]) + ' hours')
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
"""Input and output functions used by pyBalloon"""

from shapely.geometry import Point, mapping
from shapely.ops import transform
from functools import partial
from astropy.io import ascii
from glob import glob
import pygrib as pg
import numpy as np
import requests
import pyproj
import json
import csv
import sys
import os

import param_file as p
import pyb_aux

def read_gfs_file(fname, area=None, alt0=0, t_0=None, extra_data=None, descent_only=False, step=100):

	"""

	Collect relevant information from GFS model.

	Required arguments:
		- fname -- File to read

	Optional arguments:
		- area -- tuple of NW and SE limits of the area to be read,
		eg.  (62., 22., 59., 24.). Default: None (read all data)
		NOTE top, left, bottom, right order!
		- alt0 -- Starting altitude above sea level, default: 0.0 m
		- t_0 -- Temperature at ground level, default: None
		- extra_data -- Add highest levels from extra_data to GFS
		ensembles which are limited to 10 hPa.

	Return the following data for the closest pixel location as Numpy
	arrays in a dictionary:
		- u_winds -- U component of wind [m/s] - E-W component,
		positive *towards* East
		- v_winds -- V component of wind [m/s] - N-S component,
		positive *towards* North
		- temperatures -- Temperature [K]
		- altitudes -- Geopotential height [m]
		
	"""

	g_0 = 9.80665 # m/s surface acc.
	R0 = 8.3144621 # Ideal gas constant, J/(mol*K)
	M_air = 0.0289644 # molar mass of air [kg/mol], altitude dependence
	T0 = 288.15 # K
	
	if area is not None:
		tlat, llon, blat, rlon = area
	else:
		print('Do you really wish to search the entire planet?')
		tlat, llon, blat, rlon = 90., 0., -90., 360.

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

	# Find closest pixel location ### FIXXXX
	locs = pyb_aux.all_and([lats <= tlat, lats >= blat, lons2 <= rlon, lons2 >= llon])

	row_idx, col_idx = np.where(locs)
	lats = lats[row_idx, col_idx]
	lons = lons[row_idx, col_idx]

	if len(lats) == 0: print( 'Warning! lats is empty!')
	if len(lons) == 0: print( 'Warning! lons is empty!')

	# Collect U component of wind data
	u_wind = {}
	for msg in u_msgs:
		if msg.typeOfLevel == 'isobaricInhPa':
			# print(msg.level, np.std(msg.values), msg['standardDeviation'])
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

	# Collect data to correct altitude order. Set "surface" values before real data.
	if not descent_only: 
		u_winds = [np.zeros(lats.shape)]
		v_winds = [np.zeros(lats.shape)]
	else:
		u_winds = []
		v_winds = []

	# Use given surface temperature if available, otherwise use the model value
	if not descent_only:
		if t_0 is None:
			temperatures = [t_surface]
		else:
			temperatures = [t_0*np.ones(lats.shape)]
	else:
		temperatures = []

	if not descent_only:
		altitudes = [alt0*np.ones(lats.shape)]
	else:
		altitudes = []

	property_array = np.array([u_wind, v_wind, temperature, altitude])
	key_lengths = np.array([len(u_wind.keys()), len(v_wind.keys()), len(temperature.keys()), len(altitude.keys())])
	index = np.where(key_lengths == min(key_lengths))[0][0]

	pressures = list(property_array[index].keys())

	if not descent_only:
		pressures.append(max(pressures))
	else:
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
	else:
		ind = 0

	i = 0
	for key in pressures:
		if i != ind:
			uwnd, vwnd, temp, alt = [], [], [], []
			uwnd.append(u_wind[key])
			vwnd.append(v_wind[key])
			temp.append(temperature[key])
			alt.append(altitude[key])

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
					j += 1

			u_winds.append(np.hstack(uwnd))
			v_winds.append(np.hstack(vwnd))
			temperatures.append(np.hstack(temp))
			altitudes.append(np.hstack(alt))

			i+=1
		else:
			i+=1

	if descent_only:

		alts_mean = np.array([np.mean(alt) for alt in altitudes])
		alts_min = np.array([np.min(alt) for alt in altitudes])

		diff = np.abs(alts_min - alt0)
		grid_i2, = np.where(diff == diff.min())
		grid_i2 = grid_i2[0]

		if alts_min[grid_i2] > alt0:
			altitudes.insert(grid_i2, alt0*np.ones(lats.shape))
			u_winds.insert(grid_i2, np.mean([u_winds[grid_i2], u_winds[grid_i2 - 1]])*np.ones(lats.shape))
			v_winds.insert(grid_i2, np.mean([v_winds[grid_i2], v_winds[grid_i2 - 1]])*np.ones(lats.shape))
			temperatures.insert(grid_i2, np.mean([temperatures[grid_i2], temperatures[grid_i2 - 1]])*np.ones(lats.shape))
		else:
			altitudes.insert(grid_i2 + 1, alt0*np.ones(lats.shape))
			u_winds.insert(grid_i2 + 1, np.mean([u_winds[grid_i2 + 1], u_winds[grid_i2]])*np.ones(lats.shape))
			v_winds.insert(grid_i2 + 1, np.mean([v_winds[grid_i2 + 1], v_winds[grid_i2]])*np.ones(lats.shape))
			temperatures.insert(grid_i2 + 1, np.mean([temperatures[grid_i2 + 1], temperatures[grid_i2]])*np.ones(lats.shape))

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
			data['pressures'][grid_i2] = data['pressures'][grid_i2]/np.exp(-((g_0*M_air)/(T0*R0))*(alts_min[grid_i2] - alt0))
		else:
			data['pressures'][grid_i2 + 1] = data['pressures'][grid_i2]/np.exp(-((g_0*M_air)/(T0*R0))*(alts_min[grid_i2] - alt0))
			
	return data

def read_gfs_set(directory1, directory2, area=None, alt0=0, main='gfs_main.grib2',
				 ens_main='ens_main.grib2',
				 ens_member_pattern='ens_??.grib2', descent_only=True,
				 use_extra=False):

	"""Read a set of GFS data consisting of 0.5 degree main run, 1
	degree ensemble main and 1 degree ensemble members. The 1 degree
	data are extended to 7 - 1 hPa levels with data from the main run.

	Required arguments:
		- directory -- Directory containing the data

	Optional arguments:
		- area -- tuple of NW and SE limits of the area to be read,
		eg.  (62., 22., 59., 24.). Default: None (read all data)
		- alt0 -- Starting altitude of the balloon, meters from the WGS84
		reference ellipsoid. Default: 0.0.
		- main -- Filename of the GFS main run. Default: gfs_main.grib2.
		- ens_main -- Filename of the GFS ensemble main run. Default:
		ens_main.grib2.
		- ens_member_pattern -- Filename pattern of the ensemble runs.
		Default: ens_??.grib2

	Return:
		- List of dictionaries containing u_winds, v_winds, temperatures,
		pressures and altitudes.
	"""

	all_data = []

	fname = os.path.join(directory1, main)
	print( "Reading GFS operational run from", fname)
	main_run_data = read_gfs_file(fname, area=area, alt0=alt0, descent_only=descent_only)
	all_data.append(main_run_data)

	if ens_main is not None:
		fname = os.path.join(directory2, ens_main)
		print( "Reading GFS ensemble main run from", fname)
		if use_extra:
			ens_main_data = read_gfs_file(fname, 
										  area=area,
										  alt0=alt0, 
										  extra_data=main_run_data, descent_only=descent_only)
		else:
			ens_main_data = read_gfs_file(fname, 
										  area=area,
										  alt0=alt0, 
										  extra_data=None, descent_only=descent_only)
		all_data.append(ens_main_data)

	if ens_member_pattern is not None:
		ens_files = glob(os.path.join(directory2, ens_member_pattern))
		ens_files.sort()
	
		for fname in ens_files:

			if '_00.' not in fname:
				print( "Reading GFS ensemble member from", fname)

				if use_extra:
					all_data.append(read_gfs_file(fname, 
												  area=area, 
												  alt0=alt0, 
												  extra_data=main_run_data, descent_only=descent_only))
				else:
					all_data.append(read_gfs_file(fname, 
												  area=area, 
												  alt0=alt0, 
												  extra_data=None, descent_only=descent_only))

	return all_data


def read_gfs_single(directory=None, area=None, alt0=0, descent_only=False, step=100.):
	"""Read a set of 0.5 degree GFS data.

	Required arguments:
		- directory -- Directory containing the data

	Optional arguments:
		- area -- tuple of NW and SE limits of the area to be read,
		eg.  (62., 22., 59., 24.). Default: None (read all data).
		Note: top, left, bottom, right order.
		- alt0 -- Starting altitude of the balloon, meters from the WGS84
		reference ellipsoid. Default: 0.0.

	Return:
		- List of dictionaries containing u_winds, v_winds, temperatures,
		pressures and altitudes.

	"""

	all_data = []

	fname = os.path.join(directory, (directory + '.grb2'))
	print('Reading GFS data from ' + fname[-28:])
	main_run_data = read_gfs_file(fname, area=area, alt0=alt0, descent_only=descent_only, step=step)
	all_data.append(main_run_data)

	return all_data


def get_sounding(station_id, date, utc_hour):
	"""Get sounding data using station ID, date and time.

	Data are retrived from Universio of Wyoming. Valid times depend on
	the station, normal sounding times are 00 and 12 UTC, or 06 and
	18 UTC.

	Required arguments:
		- station_id -- ID of the station: number or air-port code
		- date -- Date as 3-tuple (yyyy, mm, dd), eg. (2013, 3, 18)
		- utc_hour -- UTC hour of the requested sounding eg. 6

	Returns:
		- Dictionary containing station coordinates, 'u_winds',
		  'v_wind', 'temperatures', 'altitudes' and 'pressures' as
		  Numpy arrays
	"""

	knots2ms = 0.514

	year, month, day = date[0], date[1], date[2]

	url = 'http://weather.uwyo.edu/cgi-bin/sounding?' + \
		'TYPE=TEXT%3ALIST&YEAR=' + str(year) + '&MONTH=' + \
		'%02d' % month + '&FROM=' + '%02d%02d' % (day, utc_hour) + \
		'&TO=' + '%02d%02d' % (day, utc_hour) + '&STNM=' + str(station_id)

	req = requests.get(url)
	text  = req.text
	station_lat = float(text.split('Station latitude: ')[-1].split('\n')[0])
	station_lon = float(text.split('Station longitude: ')[-1].split('\n')[0])
	data = text.split('<PRE>')[1].split('</PRE>')[0].split('\n')
	data = data[5:] # only the numerical rows
	
	pressures = []
	altitudes = []
	temperatures = []
	u_winds = []
	v_winds = []

	for row in data:
		nums = row.split()
		if len(nums) == 11:
			pressures.append(float(nums[0]))
			altitudes.append(float(nums[1]))
			temperatures.append(273.15+float(nums[2]))
		
			wdir = np.radians(float(nums[6]))
			wspd = float(nums[7])*knots2ms
		
			# Towards North and East are positive
			u_winds.append(-1*wspd*np.sin(wdir))
			v_winds.append(-1*wspd*np.cos(wdir))

	data = {}
	data['lats'] = np.array(station_lat)
	data['lons'] = np.array(station_lon)
	data['u_winds'] = np.array(u_winds)
	data['v_winds'] = np.array(v_winds)
	data['temperatures'] = np.array(temperatures)
	data['pressures'] = 100*np.array(pressures)
	data['altitudes'] = np.array(altitudes)

	return data


def read_live_data(fname, delimiter=',', temp_conv=(1, 0),
				   pressure_conv=(1, 0)):
	"""Read data from live feed (via a file) from a balloon (or
	simulator). The data are assumed to be in CSV format with the
	following fields:
		- latitude, longitude, altitude
	or
		- latitude, longitude, altitude, pressure
	or
		- latitude, longitude, altitude, pressure, temperature

	Location data are compulsory (first 3 fields). Pressure and
	temperature fields can be empty, missing or in use.

	Required arguments:
		- fname -- Filename of the live data.

	Optional arguments:
		- delimiter -- CSV field delimiter. Default: ','
		- temp_conv -- Conversion factors (multiplier and offset)
		between temperature data and temperature in Kelvins. Default
		(Kelvin -> Kelvin): (1, 0).
		- pressure_conv -- Conversion factors (multiplier and offset)
		between
		- pressure data and pressure in Pascals. Default: (Pa -> Pa):
		(1, 0)

	Return:
		- None or dictionary of Numpy arrays containing the data, if
		any available.
	"""

	try:
		fid = open(fname)
	except:
		return None

	reader = csv.reader(fid, delimiter=delimiter)

	lats = []
	lons = []
	altitudes = []
	pressures = []
	temperatures = []
	num = 0
	for row in reader:
		if len(row) == 3:
			lat, lon, alt = row
			pres, temp = None, None
		else:
			if len(row) == 4:
				lat, lon, alt, pres = row
				temp = None
			else:
				lat, lon, alt = row[0], row[1], row[2]
				pres, temp = row[3], row[4]
		
		lats.append(float(lat))
		lons.append(float(lon))
		altitudes.append(float(alt))
		pressures.append(float(pres))
		temperatures.append(float(temp))
		num += 1

	if num == 0:
		return None

	live_data = {}
	live_data['lats'] = np.array(lats)
	live_data['lats'] = np.array(lons)
	live_data['altitudes'] = np.array(altitudes)
	live_data['pressures'] = np.array(pressures) * pressure_conv[0] + pressure_conv[1]
	live_data['temperatures'] = np.array(temperatures) * temp_conv[0] + temp_conv[1]

	return live_data

def save_kml(fname, data, model_start_idx=0, eps_mode='end', other_info=None, params=None, radius=5000):
	"""Save given trajectories as KML. The first trajectory is assumed
	to be from GFS main run, the second from ensemble main and the
	rest from other ensemble members.

	Required arguments:
		- fname -- File to save the KML data to
		- data -- List of dictionaries containing 'lats', 'lons',
		'altitudes' and 'times.

	Optional arguments
		- model_start_idx -- Vector index to show where model
		trajectory starts. Trajectory before this is from live
		data. Default: 0
		- eps_mode -- How to present ensemble predictions of
		trajectories. Possible modes are 'end' and 'full'. Default:
		'end'
		- other_info -- Additional information to save into the KML
		path: eg.  other_info = [(lat, lon, alt, 'marker name', 
								  'marker description')]

	Todo:
		- modelled current location of the balloon
		- latest live data location
		- important times / altitudes / ...
		- maximum modelled altitude / maximum gained altitude (balloon
		burst at T+xx:xx:xx, xx.x km)
		- different colors for model / live data parts of the trajectory?
			- change color of the trajectory at model_start_idx
			- or different paths?
	"""

	if params == None:
		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		interpolate = p.interpolate
		drift_time = p.drift_time

	else:
		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		interpolate = bool(params[-2])
		drift_time = float(params[-1])

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
			# if data['times'][i] >= t_prev + 0.5:
			#	 t_prev = data['times'][i]
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
			# kml_str += '<name>'+dat[3]+'</name>\n' ## change this
			# kml_str += '<description>'+dat[4]+'</description>\n'
			kml_str += '<name>initial condition</name>\n' ## change this
			kml_str += '<description>initial condition</description>\n'
			kml_str += '<Point>\n'
			kml_str += '<altitudeMode>absolute</altitudeMode>\n'
			kml_str += '<coordinates>%f,%f,%f</coordinates>\n' % \
				(dat[1], dat[0], dat[2])
			kml_str += '</Point>\n'
			kml_str += '</Placemark>\n'

	if descent_only:

		end_dir = fname[:42] + 'Endpoints/' + fname[52:78]
		efile = 'endpoint_' + fname[78:-4] + '.dat'

	else:

		end_dir = fname[:42] + 'Endpoints/' + fname[52:67]
		efile = 'endpoint_' + fname[67:-4] + '.dat'

	data = ascii.read(end_dir + efile)
	end_lat = data['lat'][-1]
	end_lon = data['lon'][-1]
	kml_str_add = create_circle(lon = end_lon, lat = end_lat, radius=radius)
	kml_str += kml_str_add

	kml_str += '</Document>\n'
	kml_str += '</kml>\n'

	fid = open(fname, 'w')
	fid.write(kml_str)
	fid.close()

def create_circle(lon=None, lat=None, radius=5000):

	point = Point(lon, lat)

	local_azimuthal_projection = '+proj=aeqd +R=6371000 +units=m +lat_0={point.y} +lon_0={point.x}'

	wgs84_to_aeqd = partial(pyproj.transform, pyproj.Proj('+proj=longlat +datum=WGS84 +no_defs'), pyproj.Proj(local_azimuthal_projection), )
	aeqd_to_wgs84 = partial(pyproj.transform, pyproj.Proj(local_azimuthal_projection), pyproj.Proj('+proj=longlat +datum=WGS84 +no_defs'), )

	point_transformed = transform(wgs84_to_aeqd, point)

	buffer = point_transformed.buffer(radius)

	buffer_wgs84 = transform(aeqd_to_wgs84, buffer)
	coords = json.dumps(mapping(buffer_wgs84))

	from io import StringIO
	io = StringIO(unicode(coords))
	coords = json.load(io)['coordinates'][0]

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

def merge_kml(datestr, params=None, drift_times=None):

	if params == None:
		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		interpolate = p.interpolate
		drift_time = p.drift_time

	else:

		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		interpolate = bool(params[-2])
		drift_time = float(params[-1])

	ext_str = '.'

	if datestr == '20180406_1':
		ext_str = '_8.'
	elif datestr == '20180406_2':
		ext_str = '_18.'

	dir_base = '/home/ellen/Desktop/SuperBIT/Weather_data/'

	if descent_only:
		ext = 'descent_only/start_point' + next_point + '/'

	else:
		ext = 'ascent+descent/'

	dir1 = dir_base + 'kml_files/' + ext

	lines = {}
	for fname in os.listdir(dir1):
		if datestr in fname and not 'merged' in fname and float(fname[-11:-7]) in drift_times:
			if interpolate:
				if 'interpolated' in fname:
					if ext_str in fname:
						print(fname)
						fname0 = fname
						lines[int(fname[-10:-7])] =  [line.rstrip('\n').split(' ') for line in open(dir1 + fname)]
			else:
				if not 'interpolated' in fname:
					if ext_str in fname:
						fname0 = fname
						lines[int(fname[-10:-7])] =  [line.rstrip('\n').split(' ') for line in open(dir1 + fname)]

			
	keys = list(lines.keys())

	kml_str1 = '<?xml version="1.0" encoding="UTF-8"?>\n'
	kml_str1 += '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">\n'
	kml_str1 += '<Document id="feat_2">\n'


	key0 = keys[0]

	for key in keys:
		info = lines[key]

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
				if key != key0:
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

	kml_fname1 = 'merged_' + fname0[:-13] + '_' + str(min(keys)) + '-' + str(max(keys)) + 'min.kml'

	fid = open(dir1 + kml_fname1, 'w')
	fid.write(kml_str1)
	fid.close()

	dir2 = dir_base + 'Endpoints/' + ext

	endpoints = {}
	for fname in os.listdir(dir2):
		if datestr in fname:
			if interpolate:
				if 'interpolated' in fname:
					data = ascii.read(dir2 + fname)
					endpoint = [data['lat'][1], data['lon'][1], data['alt'][1]]
					endpoints[int(fname[-11:-7])] = endpoint
			else:
				if not 'interpolated' in fname:
					data = ascii.read(dir2 + fname)
					endpoint = [data['lat'][1], data['lon'][1], data['alt'][1]]
					endpoints[int(fname[-11:-7])] = endpoint

	kml_fname2 = 'endpoints_merged_' + fname0[:-13] + '_' + str(min(keys)) + '-' + str(max(keys)) + 'min.kml'

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

	for key in keys:
		kml_str2 += str(endpoints[key][1]) + ',' + str(endpoints[key][0])  + ',' + str(endpoints[key][2]) + '\n'

	kml_str2 += '</coordinates>\n'
	kml_str2 += '<extrude>1</extrude>\n'
	kml_str2 += '<tessellate>1</tessellate>\n'
	kml_str2 += '<altitudeMode>relativeToGround</altitudeMode>\n'
	kml_str2 += '</LineString>\n'
	kml_str2 += '</Placemark>\n'

	for key in keys:
		kml_str2 += '<Placemark>\n'
		kml_str2 += '<name>End, drift = '+ str(key) + ' min.</name>\n'
		kml_str2 += '<Point>\n'
		kml_str2 += '<coordinates>' + str(endpoints[key][1]) + ',' + str(endpoints[key][0]) + '</coordinates>' + '\n'
		kml_str2 += '</Point>\n'
		kml_str2 += '</Placemark>\n'

		kml_str_add = create_circle(lon = endpoints[key][1], lat = endpoints[key][0])
		kml_str2 += kml_str_add

	kml_str2 += '</Document>\n'
	kml_str2 += '</kml>'

	fid = open(dir1 + kml_fname2, 'w')
	fid.write(kml_str2)
	fid.close()

# method to match prediction to gps file
def match_pred2gps(date):

	datestr = date[2:4] + '-' + date[4:6] + '-' + date[6:]

	if datestr[3] == '0':
		datestr = datestr[:3] + datestr[4:]
		if datestr[5] == '0':
			datestr = datestr[:5] + datestr[6:]
	else:
		if datestr[6] == '0':
			datestr = datestr[:6] + datestr[7:]

	return datestr

# method to match gps to prediction file
def match_gps2pred(date):

	if '_' in date:
		i0 = [m.start() for m in re.finditer('_', date)][0]
	else:
		i0 = len(date)

	inds = [m.start() for m in re.finditer('-', date)]
	i1, i2 = inds[0], inds[1]

	if i2 - i1 == 2:
		month = '0' + date[i1+1]
	else:
		month = date[i1+1:i1+3]

	if i0 - i2 == 2:
		day = '0' + date[i2+1:]
	else:
		day =  date[i2+1:]

	return '20' + date[:2] + month + day

"""
Functions used by pyBalloon to calculate balloon trajectories
"""

import numpy as np
import datetime
import sys, os
import math
import time

import get_gfs
import pyb_aux
import pyb_io

import param_file as p

#################################################################################################################

def update_files(used_weather_files=None, data=None, lat_rad=None, lon_rad=None, all_alts=None, current_time=None, balloon=None, datestr=None, utc_hour=None, loc0=None, \
	total_time=[0], index=0, params=None):
	"""
	Update weather files to newest and best available (closest in time to current time)

	Arguments
	=========
	used_weather_files : dict
		Dictionary of weather files that have been used so far. The keys represent the times at which the weather files started being used.
	data : dict
		Dictionary containing weather forecast data plus properties calculated using the calc_properties() function
	lon_rad : float
		Current longitude in radians
	lat_rad : float
		Current latitude in radians
	all_alts : list
		List containing altitudes of trajectory so far
	current_time : float
		Current time of the trajectory
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	datestr : string
		Date of initial point
	utc_hour : float
		Time of initial point
	loc0 : floats in tuple
		(latitude in degrees, longitude in degrees, altitude in km) of initial point
	total_time : list
		List of increments in time for steps in trajectory so far
	index : int
		Index corresponding to current weather file being used (only when not interpolating)
	params : list
	   	List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	Return:
		Data of new weather file(s), updated weather file keys, updated index, and updated used_weather_files.
	"""

	descent_only, drift_time, resolution, hr_diff, check_elevation, live, params, balloon = pyb_io.set_params(params=params, balloon=balloon)
	res = str(int(-4*resolution + 6))

	weather_files = list(data.keys())
	weather_files.sort()
	keys = weather_files

	time_keys = [(int(int(file[15:19])/100.) + int(file[20:23])) for file in weather_files]
	time_hhh = [int(file[20:23]) for file in weather_files]

	max_datestr, max_time, max_hhhh, max_hhh = keys[-1][6:14], time_keys[-1], keys[-1][15:19], [-1][20:23]

	############################################################################################################

	if (datestr == max_datestr and current_time > max_time and hr_diff == 0) or (current_time > max_time % 24 and hr_diff > 0):

		sys.stdout.write('\r')
		sys.stdout.flush()
		sys.stdout.write('Adding new weather file...'.ljust(60) + '\r')
		sys.stdout.flush()
		time.sleep(0.2)

		new_weather_files = get_gfs.get_interpolation_gfs_files(datestr=datestr, utc_hour=current_time, resolution=resolution, hr_diff=hr_diff, live=live)	
		new_hhhh, new_hhh1 = int(int(new_weather_files[-1][15:19])/100), int(new_weather_files[-1][21:24])
		new_hhhh0, new_hhh01 = int(int(new_weather_files[0][15:19])/100), int(new_weather_files[0][21:24])

		if new_weather_files[-2] not in list(data.keys()):
			data[new_weather_files[-2]] = prepare_data(weather_file=new_weather_files[-2], loc0=loc0, current_time=current_time, balloon=balloon, descent_only=descent_only)

		data[new_weather_files[-1]] = prepare_data(weather_file=new_weather_files[-1], loc0=loc0, current_time=current_time, balloon=balloon, descent_only=descent_only)

		used_weather_files[current_time] = new_weather_files

		keys = list(data.keys())
		keys.sort()
		
		index = np.where(np.array(keys) == new_weather_files[-1])[0][0]

	############################################################################################################

	weather_files = list(data.keys())
	weather_files.sort()

	keys.sort()

	return data, keys, index, used_weather_files

###############################################################################################################

def calc_time_frac(current_time=None, weather_files=None):
	"""
	Calculate how far the current time is from each interpolation file as fractions of difference in time over total time.

	E.g. if the current time is 10 then on the lhs we have 6 and on the rhs we have 12, 10 is then 2/3 of the way away from 6 and 1/3 away from 12.
	the fractions needed are then 1-2/3 for 6 and 1 - 1/3 for 12

	Arguments
	=========
	current_time : float
		Current time of trajectory
	weather_files : list
		List of the interpolation files currently being used
	Return:
		Two tuples containing the interpolation file name and corresponding fraction
	"""

	weather_files.sort()
	times = [(int(int(file[15:19])/100.) + int(file[20:23])) for file in weather_files]

	earlier_file = weather_files[-2]
	later_file = weather_files[-1]

	earlier_time = times[-2]
	later_time = times[-1]

	if later_time !=  24:
		dt1 = current_time - (earlier_time % 24)
		dt2 = (later_time % 24) - current_time
	else:
		dt1 = current_time - (earlier_time % 24)
		dt2 = later_time - current_time		

	dt_total = 3.0

	frac1 = 1. - dt1/dt_total
	frac2 = 1. - dt2/dt_total

	return (earlier_file, frac1), (later_file, frac2)

#################################################################################################################

def read_data(loc0=None, weather_file=None, descent_only=False):
	"""
	Read in model forecast data

	Arguments
	=========
	loc0 : floats in tuple
		(latitude in degrees, longitude in degrees, altitude in km) of initial point
	weather_file : string
		Name of weather file from data is to be read
	descent_only : bool
		Option to start the trajectory at its highest point. If True, the trajectory only has a descent phase
	Return:
		Dictionary containing the forecast data
	"""

	lat0, lon0, alt0 = loc0

	in_dir = p.path + p.weather_data_folder + p.GFS_folder
	usb_dir = p.usb_path

	tile_size = p.tile_size # degrees (read a tile this wide/high from the GFS grb2 file)
	area = (lat0 + (tile_size/2.), lon0 - (tile_size/2.), lat0 - (tile_size/2.), lon0 + (tile_size/2.)) # note top, left, bottom, right ordering for area

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Reading GFS data from ' + str(weather_file) + '.grb2...'.ljust(20))
	sys.stdout.flush()

	file_dir = in_dir
	if not os.path.isfile(in_dir + weather_file + '.grb2') and get_gfs.check_connection_device(p.usb_device_name):
		file_dir = usb_dir
		
	model_data = pyb_io.read_gfs_single(directory=file_dir + weather_file, area=area, alt0=alt0, descent_only=descent_only)[0]

	return model_data

#################################################################################################################

def movement2ll(lat_rad=None, lon_rad=None, alt=None, dx=None, dy=None):
	"""
	Calculate new lat/lon coordinates given original location and Cartesian displacements.

	Arguments
	=========
	lat_rad : float
		Current latitude in radians
	lon_rad : float
		Current longitude in radians
	alt : float
		Current altitude in meters
	dx : float
		East - west movement in meters (East is positive)
	dy : float
		North - south movement in meters (North is positive)
	Return:
		New coordinates: latitude [radians], longitude [radians], and distance traveled [km]
	"""

	# radius derived from WGS84 reference ellipsoid with altitude added to it.
	Re = pyb_aux.earth_radius(lat_rad)
	radius = Re + alt/1000. # convert alt to km
	
	# calculate distance travelled
	dist = np.sqrt(dx*dx + dy*dy)/1000. # Convert to km

	# calculate direction of movement
	theta = np.arctan2(dx, dy)

	cos_dr = np.cos(dist/radius)
	sin_dr = np.sin(dist/radius)
	sin_lat = np.sin(lat_rad)
	cos_lat = np.cos(lat_rad)
	cos_theta = np.cos(theta)
	sin_theta = np.sin(theta)

	# use haversine formula
	lat2 = np.arcsin(sin_lat * cos_dr + cos_lat * sin_dr * cos_theta)
	lon2 = lon_rad + np.arctan2(sin_theta * sin_dr * cos_lat, cos_dr - sin_lat * np.sin(lat2))

	return lat2, lon2, dist

#################################################################################################################

def calc_properties(data=None, weather_file=None, loc0=None, step=None, balloon=None, descent_only=False):
	"""
	Calculate necessary properties, e.g. air densities and descent speeds, and carry the interpolation to more altitudes than in forecast data

	Arguments
	=========
	data : dict
		Dictionary containing data from forecast model read in by read_data()
	weather_file : string
		Name of weather file from which the data is to be read (if data is None)
	loc0 : floats in tuple
		(latitude in degrees, longitude in degrees, altitude in km) of initial point
	step : float
		Size of the altitude step
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	descent_only : bool
		Option to start the trajectory at its highest point. If True, the trajectory only has a descent phase
	Return:
		Data appended with new properties
	"""

	if balloon == None:
		balloon = p.balloon

	if step == None:
		step = balloon['altitude_step']

	if data == None:
		data = read_data(loc0=loc0, weather_file=weather_file, descent_only=descent_only)

	lat0, lon0, alt0 = loc0

	data['air_densities'] = pyb_aux.air_density(data)

	if not descent_only:

		data['balloon_radii'], gas_mass = pyb_aux.mooney_rivlin(data, balloon['radius_empty'], balloon['fill_radius'], balloon['thickness_empty'])
		data['balloon_volumes'] = pyb_aux.balloon_volume(data)
		total_mass = balloon['equip_mass'] + balloon['balloon_mass'] + gas_mass # kg
		data['lifts'] = pyb_aux.lift(data, total_mass)
		data['ascent_speeds'] = pyb_aux.ascent_speed(data, total_mass, balloon['Cd_balloon'])
		if 'simple_ascent_rate' in balloon:
			data['ascent_speeds'] = np.ones_like(data['ascent_speeds']) * balloon['simple_ascent_rate']

	data['descent_speeds'] = pyb_aux.descent_speed(data, balloon['equip_mass'], balloon['Cd_parachute'], balloon['parachute_areas'], balloon['altitude_step'], balloon['parachute_change_altitude'])

	data = pyb_aux.data_interpolation(data=data, alt0=alt0, step=step, descent_only=descent_only)

	return data

#################################################################################################################

def calc_displacements(data=None, balloon=None, step=None, descent_only=False):
	"""
	Calculate the displacements in the x/y directions for all possible locations (lon/lat/alt) in weather data

	Arguments
	=========
	data : dict
		Dictionary containing weather forecast data plus properties calculated using the calc_properties() function
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	step : float
		Size of the altitude step
	descent_only : bool
		Option to start the trajectory at its highest point. If True, the trajectory only has a descent phase
	Return:
		Data appended with additional calculated properties  		
	"""

	if balloon == None:
		balloon = p.balloon

	if step == None:
		step = balloon['altitude_step']

 	# ascent properties		
	if not descent_only:
		for i in range(len(data['ascent_speeds'])):
			data['ascent_speeds'][i] = np.array(list(0.5*data['ascent_speeds'][i][:-1] + 0.5*data['ascent_speeds'][i][1:]) + [data['ascent_speeds'][i][-1]])

	if not descent_only:
		data['max_altitudes'], data['max_alt_idxs'] = pyb_aux.burst_altitude(data, balloon['burst_radius'])
		delta_t = step / data['ascent_speeds']
		data['ascent_time_steps'] = delta_t
		data['cumulative_ascent_times'] = np.cumsum(delta_t)/60.

 	# descent properties		
	for i in range(len(data['descent_speeds'])):
		data['descent_speeds'][i] = np.array(list(0.5*data['descent_speeds'][i][:-1] + 0.5*data['descent_speeds'][i][1:]) + [data['descent_speeds'][i][-1]])
		data['u_winds'][i] = np.array(list(0.5*data['u_winds'][i][:-1] + 0.5*data['u_winds'][i][1:]) + [data['u_winds'][i][-1]])
		data['v_winds'][i] = np.array(list(0.5*data['v_winds'][i][:-1] + 0.5*data['v_winds'][i][1:]) + [data['v_winds'][i][-1]])
	
	delta_t = -1*step / data['descent_speeds']

	# # for logaritmic steps!
	# steps = list(data['altitudes'][1:] - data['altitudes'][:-1])
	# steps = np.array([steps[0]] + steps)
	# delta_t = []
	# for i in range(len(data['descent_speeds'])):
	# 	delta_t.append(-1*steps/data['descent_speeds'][i])

	data['descent_time_steps'] = np.array(delta_t)
	data['cumulative_descent_times'] = np.cumsum(delta_t)/60

	return data

#################################################################################################################

def calc_variable(grid_i, i, lon_rad, lat_rad, data, fracs, prop, resolution):
	"""
	Calculate the value of a variable at a given grid location/altitude

	Arguments
	=========
	grid_i : int
		index of list representing current lon/lat
	i : int
		index of list representing current altitude
	lon_rad : float
		Current longitude in radians
	lat_rad : float
		Current latitude in radians
	data : dict
		Dictionary containing weather forecast data plus properties calculated using the calc_properties() function
	fracs : list
		List containing the fractions used for the interpolation between weather files
	prop : string
		Name of variable to be evaluated
	resolution : float
		Resolution of the weather forecasts
	Return:
		Value of given variable  		
	"""

	t1, f1, t2, f2 = fracs
	keys = list(data.keys())
	data_lats, data_lons = np.radians(data[keys[0]]['lats']), np.radians(data[keys[0]]['lons'])

	x1, y1, x2, y2, low_left, up_left, low_right, up_right = pyb_aux.find_bilinear_points(grid_i, i, lon_rad[-1], lat_rad[-1], data_lons, data_lats, resolution)
	coords, inds = [x1, y1, x2, y2], [low_left, up_left, low_right, up_right]

	var = f1*pyb_aux.bilinear_interpolation(i=i, lon_rad=lon_rad[-1], lat_rad=lat_rad[-1], prop=data[t1][prop], coords=coords, inds=inds)
	var += f2*pyb_aux.bilinear_interpolation(i=i, lon_rad=lon_rad[-1], lat_rad=lat_rad[-1], prop=data[t2][prop], coords=coords, inds=inds)

	return var

#################################################################################################################

def calc_movements(data=None, used_weather_files=None, ini_conditions=None, params=None, balloon=None):
	"""
	Calculate the trajectory of the balloon/parachute given a start position & other input parameters

	Arguments
	=========
	used_weather_files : dict
		Dictionary of initial weather files used. The keys represent the times at which the weather files started being used (i.e. here the starting time)	
	ini_conditions : tuple of strings
		(Date of initial point, Initial time of trajectory, (latitude in degrees, longitude in degrees, altitude in km) of initial point
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	Return:
		Output trajectory data, and dictionary of used_weather_files
	"""

	############################################################################################################

	# set general parameters and initial conditions
	datestr, utc_hour, loc0 = ini_conditions
	lat0, lon0, alt0 = loc0

	descent_only, drift_time, resolution, hr_diff, check_elevation, live, params, balloon = pyb_io.set_params(params=params, balloon=balloon)
	lat_rad, lon_rad, all_alts = [np.radians(lat0)], [np.radians(lon0)], [alt0]

	if lon_rad[0] < 0:
		lon_rad[0] += 2*np.pi
	if lat_rad[0] < 0:
		lat_rad[0] += 2*np.pi

	current_time = utc_hour
	
	keys = list(data.keys())
	time_key0 = (int(int(keys[0][15:19])/100.) + int(keys[0][20:23]))

	alts = data[keys[0]]['altitudes'] # alts, lats and lons are the same for all weather files (if we dont give it different areas)
	data_lats, data_lons  = np.radians(data[keys[0]]['lats']), np.radians(data[keys[0]]['lons'])
	size = int(np.sqrt(len(data_lats)))

	x_prop, y_prop = 'u_winds', 'v_winds'
	t_props = ['ascent_time_steps', 'descent_time_steps']
	speed_props = ['ascent_speeds', 'descent_speeds']

	tfutures, speeds = [], []
	total_time, dists, dists_u, dists_v = [0], [0], [0], [0]
	stage, index, props_index, i, max_i, timer = 1, 0, 0, 0, 0, 0

	if drift_time == 0:
		stage_update = 2
		out_str = 'Calculating descent...'
	else:
		stage_update = 1
		out_str = 'Calculating drift trajectory...'

	if not descent_only:

		sys.stdout.write('\r')
		sys.stdout.flush()
		sys.stdout.write('Calculating ascent...'.ljust(60) + '\r')
		sys.stdout.flush()
		time.sleep(0.2)

	# calc trajectory
	while True:

		if current_time + total_time[-1]/3600. >= 24.:

			year, month, day = int(datestr[:4]), int(datestr[4:6]), int(datestr[6:])
			date = datetime.datetime(year, month, day) + datetime.timedelta(days=1)
			datestr = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)

		current_time = (float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600) % 24

		# update weather files
		data, keys, index, used_weather_files = update_files(used_weather_files=used_weather_files, data=data, lat_rad=lat_rad, lon_rad=lon_rad, all_alts=all_alts, \
				balloon=balloon, datestr=datestr, utc_hour=utc_hour, loc0=loc0, total_time=total_time, current_time=current_time, index=index, params=params)

		if not (stage == 1 and descent_only):

			# determine fractions for the interpolation between forecasts
			(t1, f1), (t2, f2) = calc_time_frac(current_time=current_time, weather_files=keys)

			delta_t = current_time - (int(int(t1[15:19])/100.) + int(t1[20:23])) # difference between forecast and model time
			tfuture = f1*int(t1[20:23]) + f2*int(t2[20:23])

			# Find the closest grid point
			diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons - lon_rad[-1])**2)
			grid_i, = np.where(diff == diff.min())
			grid_i = grid_i[0]
			min_diff = diff[grid_i]

			dx = calc_variable(grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), x_prop, resolution)
			dy = calc_variable(grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), y_prop, resolution)

			if stage != 2:

				dt = calc_variable(grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), t_props[props_index], resolution)
				speed = calc_variable(grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), speed_props[props_index], resolution)

			else:

				dt = 5
				speed = 0
				
			dx *= dt
			dy *= dt

		if stage == 1:

			if not descent_only:

				sys.stdout.write('\r')
				sys.stdout.flush()
				sys.stdout.write(str(round(100.*(max_i - i)/max_i, 1)) + r' % done'.ljust(60) + '\r')
				sys.stdout.flush()

				if all_alts[-1] >= data[keys[index]]['max_altitudes'][grid_i]:

					sys.stdout.write('\r')
					sys.stdout.flush()
					stage += stage_update
					props_index += 1
					sys.stdout.write(out_str.ljust(60) + '\r')
					sys.stdout.flush()
					time.sleep(0.2)

					final_i = 0
					continue
			else:

				# Mimick the movement during ascent if we only want descent
				if alts[i] >= alt0:
					stage += stage_update
					props_index += 1
					final_i = i
					continue

			i += 1
			max_i = i

		elif stage == 2:

			timer += dt
			if timer == drift_time*60:
				stage += 1
				continue

		elif stage == 3:

			sys.stdout.write('\r')
			sys.stdout.flush()
			sys.stdout.write(str(round(100.*(max_i - i)/max_i, 1)) + r' % done'.ljust(60) + '\r')
			sys.stdout.flush()

			if i == 0:

				(t1, f1), (t2, f2) = calc_time_frac(current_time=current_time, weather_files=keys)

				delta_t = current_time - (int(int(t1[15:19])/100.) + int(t1[20:23])) # difference between forecast and model time
				tfuture = f1*int(t1[20:23]) + f2*int(t2[20:23])
				tfutures.append(tfuture)

				break

			i -= 1

		if stage == 2:

			lat, lon, dist = movement2ll(lat_rad=lat_rad[-1], lon_rad=lon_rad[-1], alt=all_alts[-1], dx=dx, dy=dy)
			alt = all_alts[-1]

		elif stage != 2 and not (stage == 1 and descent_only):

			lat, lon, dist = movement2ll(lat_rad=lat_rad[-1], lon_rad=lon_rad[-1], alt=alts[i], dx=dx, dy=dy)
			alt = alts[i]

		if not (descent_only and stage == 1):

			speeds.append(speed)
			lat_rad.append(lat)
			lon_rad.append(lon)
			dists.append(dist)
			total_time.append(dt)
			all_alts.append(alt)
			tfutures.append(tfuture)

	speeds.append(calc_variable(grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), speed_props[props_index], resolution))

	# get more accurate end-point based on elevation data
	if check_elevation:

		sys.stdout.write('\r')
		sys.stdout.flush()
		sys.stdout.write('Getting new endpoint based on elevation...'.ljust(60))

		elevation = pyb_aux.get_elevation(lon=np.degrees(lon_rad[-1]), lat=np.degrees(lat_rad[-1]))
		if np.abs(elevation - all_alts[-1]) > balloon['altitude_step']/10.:

			new_end_point, new_alt = pyb_aux.get_endpoint(data=(np.degrees(np.array(lat_rad)), np.degrees(np.array(lon_rad)), np.array(all_alts), np.array(dists)))

			diff = np.sqrt((np.degrees(np.array(lat_rad)) - new_end_point[0])**2 + (np.degrees(np.array(lon_rad)) - new_end_point[1])**2)
			index = np.where(diff == min(diff))[0][-1]

			lat_rad, lon_rad, all_alts = lat_rad[:index], lon_rad[:index], all_alts[:index], 
			dists, speeds, total_time = dists[:index+1], speeds[:index+1], total_time[:index+1]

			tfutures = tfutures[:index+1]

			lat_rad.append(np.radians(new_end_point[0]))
			lon_rad.append(np.radians(new_end_point[1]))
			all_alts.append(new_alt)

	############################################################################################################

	# Convert the result array lists to Numpy 2D-arrays
	output = {}
	output['lats'] = np.degrees(np.array(lat_rad)) # to decimal degrees
	output['lons'] = np.degrees(np.array(lon_rad)) # to decimal degrees
	output['alts'] = np.array(all_alts)

	output['dists'] = np.array(dists)
	output['times'] = np.cumsum(np.array(total_time))/60 # to minutes
	output['distance'] = np.sum(np.array(dists))
	output['speeds'] = np.array(speeds)

	output['tfutures'] = np.array(tfutures)
	output['mean_direction'] = pyb_aux.calc_mean_travel_direction(lon0=lon0, lat0=lat0, end_lon=float(output['lons'][-1]), end_lat=float(output['lats'][-1]))

	# print out relevant quantities
	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Trajectories calculated'.ljust(60))
	sys.stdout.write('\n')

	print('Maximum altitude: ' + str(np.max(all_alts)) + ' m')
	print('Landing location: (%.6f, %.6f)' % (output['lats'][-1], output['lons'][-1]))
	print('Mean direction of travel: ' + str(round(output['mean_direction'], 3)) + ' degrees')
	print('Flight time: %.3f min' % (output['times'][-1]) + ', distance travelled: %.1f' % output['distance'] + ' km')
	print('')

	return output, used_weather_files

#################################################################################################################

def prepare_data(weather_file=None, loc0=None, current_time=None, balloon=None, descent_only=False):
	"""
	Prepare the data from a weather file for calculating the trajectory
	
	Arguments
	=========
	weather_file : string
		Name of weather file to be used and read
	loc0 : floats in tuple
		(latitude in degrees, longitude in degrees, altitude in km) of initial point
	current_time : float
		Current time of the trajectory	
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
  	descent_only : bool
		Option to start the trajectory at its highest point. If True, the trajectory only has a descent phase
	Return:
		Dictionary containing data ready for starting the trajectory calculations
	"""

	model_data1 = read_data(loc0=loc0, weather_file=weather_file, descent_only=descent_only)
	model_data2 = calc_properties(data=model_data1, weather_file=weather_file, loc0=loc0, balloon=balloon, descent_only=descent_only)
	model_data3 = calc_displacements(data=model_data2, balloon=balloon, descent_only=descent_only)

	return model_data3

#################################################################################################################

def run_traj(weather_files=None, ini_conditions=None, params=None, balloon=None):
	"""
	Run all functions to calculate the trajectory
	
	Arguments
	=========
	weather_files : list
		List of weather files used at the start
	ini_conditions : tuple of strings
		(Date of initial point, Initial time of trajectory, (latitude in degrees, longitude in degrees, altitude in km) of initial point
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	Return:
		The calculated trajectories, list of fig dicionaries for each weather file, dictionary of used_weather_files, and list of time differences between trajectiry times and weather forecast files
	"""

	datestr, utc_hour, loc0 = ini_conditions

	descent_only, drift_time, resolution, hr_diff, check_elevation, live, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	data_array, used_weather_files = {}, {utc_hour : weather_files}
	for weather_file in weather_files:

		model_data= prepare_data(weather_file=weather_file, loc0=loc0, current_time=utc_hour, balloon=balloon, descent_only=descent_only)
		data_array[weather_file] = model_data
		del model_data

	trajectories, used_weather_files = calc_movements(data=data_array, used_weather_files=used_weather_files, ini_conditions=ini_conditions, \
		params=params, balloon=balloon)

	return trajectories, used_weather_files

#################################################################################################################
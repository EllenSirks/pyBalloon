"""
Functions used by pyBalloon to calculate balloon trajectories
"""

import numpy as np
import datetime
import time
import sys

import get_gfs
import pyb_aux
import pyb_io

import param_file as p

#################################################################################################################

def update_files(figs=None, used_weather_files=None, data=None, lat_rad=None, lon_rad=None, all_alts=None, current_time=None, balloon=None, datestr=None, utc_hour=None, loc0=None, \
	total_time=[0], time_diffs=[], index=0, params=None, output_figs=False):
	"""
	Update weather files to newest and best available (closest in time to current time)

	Arguments
	=========
	figs : list
		List of dictionaries containing figures showing the grib data before and after interpolation between altitude steps. Each dict corresponds to one weather file.
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
	time_diffs : list
		List containing difference in time between trajectory time and weather forecasts used
	index : int
		Index corresponding to current weather file being used (only when not interpolating)
	params : list
	   	List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	output_figs : bool
		If True, figures showing the grib data before and after interpolation between altitude steps are created and saved
	Return:
		Data of new weather file(s), updated weather file keys, updated index, updated figs, updated used_weather_files and updated time_diffs.
	"""

	descent_only, next_point, time_interpolate, grid_interpolate, drift_time, resolution, hr_diff, check_sigmas, params, balloon = pyb_io.set_params(params=params, balloon=balloon)
	res = str(int(-4*resolution + 6))

	weather_files = list(data.keys())
	weather_files.sort()
	keys = weather_files
	time_keys = [(int(int(file[15:19])/100.) + int(file[20:23])) for file in weather_files]
	time_hhh = [int(file[20:23]) for file in weather_files]

	max_datestr, max_time, max_hhhh, max_hhh = keys[-1][6:14], time_keys[-1], keys[-1][15:19], [-1][20:23]
	max_date = datetime.datetime(int(max_datestr[:4]), int(max_datestr[4:6]), int(max_datestr[6:]))
	diff_days = (datetime.datetime(int(datestr[:4]), int(datestr[4:6]), int(datestr[6:])) - max_date ).days

	############################################################################################################

	if (datestr == max_datestr and current_time > max_time and time_interpolate) or (not time_interpolate and datestr == max_datestr and current_time > max_time and (current_time - max_time) > 1.5)\
	 or (diff_days > 0 and hr_diff == 0):

		sys.stdout.write('\r')
		sys.stdout.flush()
		sys.stdout.write('Adding new weather file...'.ljust(60) + '\r')
		sys.stdout.flush()
		time.sleep(0.2)

		if not time_interpolate:

			hrs = get_gfs.get_closest_hr(datestr=datestr, utc_hour=current_time, hr_diff=hr_diff)
			new_hhhh, new_hhh1, new_hhh2, datestr = hrs[0], hrs[1], hrs[2], hrs[3]
			new_weather_file = 'gfs_' + res + '_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh1).zfill(3) + '.grb2'
			new_weather_files = get_gfs.get_gfs_files(weather_files=[new_weather_file])
			loc = (np.degrees(lat_rad[-1]), np.degrees(lon_rad[-1]), all_alts[-1])
			time_diffs[current_time] = [(new_hhhh + new_hhh1) % 24 - current_time]

		else:

			new_weather_files = get_gfs.get_interpolation_gfs_files(datestr=datestr, utc_hour=current_time, resolution=resolution, hr_diff=hr_diff)
			new_hhhh, new_hhh1 = int(int(new_weather_files[-1][15:19])/100), int(new_weather_files[-1][21:24])
			new_hhhh0, new_hhh01 = int(int(new_weather_files[0][15:19])/100), int(new_weather_files[0][21:24])
			time_diffs[current_time] = [(new_hhhh0 + new_hhh01) % 24 - current_time, (new_hhhh + new_hhh1) % 24 - current_time]

		data[new_weather_files[-1]], figs_dict = prepare_data(weather_file=new_weather_files[-1], loc0=loc0, current_time=current_time, balloon=balloon, descent_only=descent_only,\
		 check_sigmas=check_sigmas, output_figs=output_figs)

		figs.append(figs_dict)
		used_weather_files[current_time] = new_weather_files
		keys = list(data.keys())
		keys.sort()
		index = np.where(np.array(keys) == new_weather_files[-1])[0][0]

	############################################################################################################

	weather_files = list(data.keys())
	weather_files.sort()

	# update from e.g. 0600_006 to 1200_000
	if not time_interpolate:
		if len(total_time) != 1:
			prev_time = (utc_hour + np.cumsum(np.array(total_time)[:-1])[-1]/3600) % 24
		else:
			prev_time = utc_hour

		if prev_time > current_time:
			year, month, day = int(datestr[:4]), int(datestr[4:6]), int(datestr[6:])
			date = datetime.datetime(year, month, day) - datetime.timedelta(days=1)
			prev_datestr = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
		else:
			prev_datestr = datestr

		old_hhhh, old_hhh1, old_hhh2, old_datestr = get_gfs.get_closest_hr(datestr=prev_datestr, utc_hour=prev_time, hr_diff=hr_diff)
		new_hhhh, new_hhh1, new_hhh2, new_datestr = get_gfs.get_closest_hr(datestr=datestr, utc_hour=current_time, hr_diff=hr_diff)

		if ((old_hhhh, old_hhh1, old_hhh2, old_datestr) != (new_hhhh, new_hhh1, new_hhh2, new_datestr)) and \
		weather_files[-1] != 'gfs_' + res + '_' + new_datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh1).zfill(3):

			sys.stdout.write('\r')
			sys.stdout.flush()
			sys.stdout.write('Updating current weather file...'.ljust(60) + '\r')
			sys.stdout.flush()
			time.sleep(0.2)

			new_weather_file_name = 'gfs_' + res + '_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh1).zfill(3) + '.grb2'
			new_weather_file = get_gfs.get_gfs_files(weather_files=[new_weather_file_name])[0]

			data[new_weather_file_name], figs_dict = prepare_data(weather_file=new_weather_file, loc0=loc0, current_time=current_time, balloon=balloon, descent_only=descent_only, \
					check_sigmas=check_sigmas, output_figs=output_figs)

			figs.append(figs_dict)
			used_weather_files[current_time] = new_weather_file
			keys = list(data.keys())
			keys.sort()
			time_diffs[current_time] = [(new_hhhh + new_hhh1) % 24 - current_time]
			index = np.where(np.array(keys) == new_weather_file_name)[0][0]

	keys.sort()

	return data, keys, index, figs, used_weather_files, time_diffs

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

	dt_total = later_time - earlier_time

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

	tile_size = p.tile_size # degrees (read a tile this wide/high from the GFS grb2 file)
	area = (lat0 + (tile_size/2.), lon0 - (tile_size/2.), lat0 - (tile_size/2.), lon0 + (tile_size/2.)) # note top, left, bottom, right ordering for area

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Reading GFS data from ' + str(weather_file) + '.grb2...'.ljust(20))
	sys.stdout.flush()

	model_data = pyb_io.read_gfs_single(directory=in_dir + weather_file, area=area, alt0=alt0, descent_only=descent_only)[0]

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

def calc_properties(data=None, weather_file=None, loc0=None, step=None, balloon=None, descent_only=False, output_figs=False):
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
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	descent_only : bool
		Option to start the trajectory at its highest point. If True, the trajectory only has a descent phase
	output_figs : bool
		If True, figures showing the grib data before and after interpolation between altitude steps are created and saved
	Return:
		Data appended with new properties and dictionary containing figures (can be empty if output_figs is False)
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

	data, figs_dict = pyb_aux.data_interpolation(data=data, alt0=alt0, step=step, descent_only=descent_only, output_figs=output_figs)

	return data, figs_dict

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

def calc_variable(grid_interpolate, grid_i, i, lon_rad, lat_rad, data, fracs, prop, resolution):
	"""
	Calculate the value of a variable at a given grid location/altitude

	Arguments
	=========
	grid_interpolate : bool
		If True, variables are interpolated between grid points
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

	if grid_interpolate:

		x1, y1, x2, y2, low_left, up_left, low_right, up_right = pyb_aux.find_bilinear_points(grid_i, i, lon_rad[-1], lat_rad[-1], data_lons, data_lats, resolution)
		coords, inds = [x1, y1, x2, y2], [low_left, up_left, low_right, up_right]

		var = f1*pyb_aux.bilinear_interpolation(i=i, lon_rad=lon_rad[-1], lat_rad=lat_rad[-1], prop=data[t1][prop], coords=coords, inds=inds)
		var += f2*pyb_aux.bilinear_interpolation(i=i, lon_rad=lon_rad[-1], lat_rad=lat_rad[-1], prop=data[t2][prop], coords=coords, inds=inds)
	else:
		var = f1*data[t1][prop][grid_i, i] + f2*data[t2][prop][grid_i, i]

	return var

#################################################################################################################

def calc_movements(data=None, used_weather_files=None, time_diffs=None, datestr=None, utc_hour=None, loc0=None, params=None, balloon=None, output_figs=False):
	"""
	Calculate the trajectory of the balloon/parachute given a start position & other input parameters

	Arguments
	=========
	used_weather_files : dict
		Dictionary of initial weather files used. The keys represent the times at which the weather files started being used (i.e. here the starting time)
	time_diffs : list
		List containing difference in time between trajectory time and weather forecasts used		
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
	output_figs : bool
		If True, figures showing the grib data before and after interpolation between altitude steps are created and saved	
	Return:
		Output trajectory data, list of fig dictionaries for each used weather files, dictionary of used_weather_files, list time_diffs
	"""

	############################################################################################################

	# set general parameters and initial conditions
	descent_only, next_point, time_interpolate, grid_interpolate, drift_time, resolution, hr_diff, check_sigmas, params, balloon = pyb_io.set_params(params=params, balloon=balloon)
	lat0, lon0, alt0 = loc0
	lat_rad, lon_rad, all_alts = [np.radians(lat0)], [np.radians(lon0)], [alt0]

	if lon_rad[0] < 0:
		lon_rad[0] += 2*np.pi
	if lat_rad[0] < 0:
		lat_rad[0] += 2*np.pi

	current_time = utc_hour
	
	keys = list(data.keys())
	time_key0 = (int(int(keys[0][15:19])/100.) + int(keys[0][20:23]))
	initial_tfut = (int(int(keys[0][15:19])/100.))
	if time_interpolate:
		initial_tfut = utc_hour - initial_tfut

	alts = data[keys[0]]['altitudes'] # alts, lats and lons are the same for all weather files (if we dont give it different areas)
	data_lats, data_lons  = np.radians(data[keys[0]]['lats']), np.radians(data[keys[0]]['lons'])
	size = int(np.sqrt(len(data_lats)))

	if check_sigmas:
		data_lats_err, data_lons_err = np.radians(data[keys[0]]['lats_err']), np.radians(data[keys[0]]['lons_err'])
		sigmas_u, sigmas_v = [], []
		sigma_u_prop, sigma_v_prop = 'u_wind_errs', 'v_wind_errs'

	x_prop, y_prop = 'u_winds', 'v_winds'
	t_props = ['ascent_time_steps', 'descent_time_steps']
	speed_props = ['ascent_speeds', 'descent_speeds']

	speeds, figs = [], []
	delta_ts, tfutures, total_time, dists, dists_u, dists_v = [utc_hour - time_key0], [initial_tfut], [0], [0], [0], [0]
	stage, index, props_index, i, max_i, timer = 1, 0, 0, 0, 0, 0

	if drift_time == 0:
		stage_update = 2
		out_str = 'Calculating descent...'
	else:
		stage_update = 1
		out_str = 'Calculating drift trajectory...'

	# inital location in grid
	diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons - lon_rad[-1])**2)
	ini_grid_i, = np.where(diff == diff.min())
	loc_diffs = [diff[ini_grid_i[0]]]

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
		data, keys, index, figs, used_weather_files, time_diffs = update_files(figs=figs, used_weather_files=used_weather_files, data=data, lat_rad=lat_rad, lon_rad=lon_rad, all_alts=all_alts, \
			balloon=balloon, datestr=datestr, utc_hour=utc_hour, loc0=loc0, total_time=total_time, current_time=current_time, time_diffs=time_diffs, index=index, params=params, output_figs=output_figs)

		if not (stage == 1 and descent_only):

			# determine fractions for the interpolation between forecasts
			if time_interpolate:
				(t1, f1), (t2, f2) = calc_time_frac(current_time=current_time, weather_files=keys)
			else:
				t1, f1, t2, f2 = keys[index], 1, keys[index], 0

			delta_t = current_time - (int(int(t1[15:19])/100.) + int(t1[20:23])) # difference between forecast and model time
			tfuture = (int(int(t1[15:19])/100.))
			if time_interpolate:
				tfuture = current_time - tfuture

			# Find the closest grid point
			diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons - lon_rad[-1])**2)
			grid_i, = np.where(diff == diff.min())
			grid_i = grid_i[0]
			min_diff = diff[grid_i]

			if check_sigmas:
				diff = np.sqrt((data_lats_err - lat_rad[-1])**2 + (data_lons_err - lon_rad[-1])**2)
				grid_i_err, = np.where(diff == diff.min())
				grid_i_err = grid_i_err[0]

			dx = calc_variable(grid_interpolate, grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), x_prop, resolution)
			dy = calc_variable(grid_interpolate, grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), y_prop, resolution)

			if stage != 2:

				dt = calc_variable(grid_interpolate, grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), t_props[props_index], resolution)
				speed = calc_variable(grid_interpolate, grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), speed_props[props_index], resolution)

			else:

				dt = 5
				speed = 0
				
			dx *= dt
			dy *= dt

			if check_sigmas:

				sigma_u = calc_variable(grid_interpolate, grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), sigma_u_prop, resolution)*dt
				sigma_v = calc_variable(grid_interpolate, grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), sigma_v_prop, resolution)*dt

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

			if all_alts[-1] < 2000.:
				elevation = pyb_aux.get_elevation(lon=np.degrees(lon_rad[-1]), lat=np.degrees(lat_rad[-1]))
				if alts[i] <= elevation:
					break

			if i == 0:
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

			if check_sigmas:
				sigmas_u.append(sigma_u)
				sigmas_v.append(sigma_v)

			if not time_interpolate:
				delta_ts.append(delta_t)

			loc_diffs.append(min_diff)

	speeds.append(calc_variable(grid_interpolate, grid_i, i, lon_rad, lat_rad, data, (t1, f1, t2, f2), speed_props[props_index], resolution))

	# get more accurate end-point based on elevation data
	if np.abs(elevation - all_alts[-1]) > balloon['altitude_step']/10.:

		new_end_point, new_alt = pyb_aux.get_endpoint(data=(np.degrees(np.array(lat_rad)), np.degrees(np.array(lon_rad)), np.array(all_alts), np.array(dists)))

		diff = np.sqrt((np.degrees(np.array(lat_rad)) - new_end_point[0])**2 + (np.degrees(np.array(lon_rad)) - new_end_point[1])**2)
		index = np.where(diff == min(diff))[0][-1]

		lat_rad, lon_rad, all_alts = lat_rad[:index], lon_rad[:index], all_alts[:index], 
		dists, speeds, total_time, tfutures, loc_diffs = dists[:index+1], speeds[:index+1], total_time[:index+1], tfutures[:index+1], loc_diffs[:index+1]

		if check_sigmas:
			sigmas_u, sigmas_v = sigmas_u[:index+1], sigmas_v[:index+1]
		if not time_interpolate:
			delta_ts = delta_ts[:index+1]

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
	output['loc_diffs'] = np.array(loc_diffs)
	output['mean_direction'] = pyb_aux.calc_mean_travel_direction(lon0=lon0, lat0=lat0, end_lon=float(output['lons'][-1]), end_lat=float(output['lats'][-1]))

	if not time_interpolate:
		output['delta_ts'] = np.array(delta_ts)

	if check_sigmas:
		output['sigmas_u'] = np.array(sigmas_u)
		output['sigmas_v'] = np.array(sigmas_v)

	# print out relevant quantities
	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Trajectories calculated'.ljust(60))
	sys.stdout.write('\n')

	print('Maximum altitude: ' + str(np.max(all_alts)) + ' m')
	print('Landing location: (%.6f, %.6f)' % (output['lats'][-1], output['lons'][-1]))
	print('Flight time: %.3f min' % (output['times'][-1]) + ', distance travelled: %.1f' % output['distance'] + ' km')
	if check_sigmas:
		total_sigma = np.sqrt(np.sqrt(np.sum(output['sigmas_u']**2))**2 + np.sqrt(np.sum(output['sigmas_v']**2))**2)/1000.
		print('Sigma from ensemble forecasts: %.3f km' % total_sigma)
	print('')

	return output, figs, used_weather_files, time_diffs

#################################################################################################################

def prepare_data(weather_file=None, loc0=None, current_time=None, balloon=None, descent_only=False, check_sigmas=False, output_figs=False):
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
	check_sigmas : bool
		If True, use ensemble forecasts to keep track of standard deviations for winds in weather forecast data
	output_figs : bool
		If True, figures showing the grib data before and after interpolation between altitude steps are created and saved	
	Return:
		Dictionary containing data ready for starting the trajectory calculations, and dictionary containing figures for interpolation checks (empty if output_figs is False)
	"""

	if check_sigmas:

		model_data1 = read_data(loc0=loc0, weather_file=weather_file, descent_only=descent_only)
		err_data = pyb_aux.calc_gefs_errs(weather_file=weather_file, loc0=loc0, current_time=current_time, descent_only=descent_only)
		model_data2, figs_dict = calc_properties(data=model_data1, weather_file=weather_file, loc0=loc0, balloon=balloon, descent_only=descent_only, output_figs=output_figs)
		model_data2 = pyb_aux.add_uv_errs(main_data=model_data2, err_data=err_data)
		model_data3 = calc_displacements(data=model_data2, balloon=balloon, descent_only=descent_only)

	else:

		model_data1 = read_data(loc0=loc0, weather_file=weather_file, descent_only=descent_only)
		model_data2, figs_dict = calc_properties(data=model_data1, weather_file=weather_file, loc0=loc0, balloon=balloon, descent_only=descent_only, output_figs=output_figs)
		model_data3 = calc_displacements(data=model_data2, balloon=balloon, descent_only=descent_only)

	return model_data3, figs_dict

#################################################################################################################

def run_traj(weather_files=None, datestr=None, utc_hour=None, loc0=None, params=None, balloon=None, output_figs=False):
	"""
	Run all functions to calculate the trajectory
	
	Arguments
	=========
	weather_files : list
		List of weather files used at the start
	datestr : string
		Date of initial point
	utc_hour : float
		Time of initial point
	loc0 : floats in tuple
		(latitude in degrees, longitude in degrees, altitude in km) of initial point
	balloon : dict
		Dictionary of balloon parameters, e.g. burtsradius, mass etc.
	params : list
		List of parameters determining how the trajectory is calculated, e.g. with interpolation, descent_only etc.
	Return:
		The calculated trajectories, list of fig dicionaries for each weather file, dictionary of used_weather_files, and list of time differences between trajectiry times and weather forecast files
	"""

	used_weather_files = {}
	used_weather_files[utc_hour] = weather_files

	time_diffs = {}
	time_diffs[utc_hour] = [(int(int(weather_file[15:19])/100.) + int(weather_file[20:23])) % 24 - utc_hour for weather_file in weather_files]

	descent_only, next_point, time_interpolate, grid_interpolate, drift_time, resolution, hr_diff, check_sigmas, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	data_array = {}
	figs1 = []
	for weather_file in weather_files:

		model_data, figs_dict = prepare_data(weather_file=weather_file, loc0=loc0, current_time=utc_hour, balloon=balloon, descent_only=descent_only, check_sigmas=check_sigmas, output_figs=output_figs)
		data_array[weather_file] = model_data
		del model_data
		figs1.append(figs_dict)
		del figs_dict

	trajectories, figs2, used_weather_files, time_diffs = calc_movements(data=data_array, used_weather_files=used_weather_files, time_diffs=time_diffs, datestr=datestr, utc_hour=utc_hour, loc0=loc0, \
		params=params, balloon=balloon, output_figs=output_figs)
	
	for fig_dict in figs2:
		figs1.append(fig_dict)
		del fig_dict

	del figs2

	return trajectories, figs1, used_weather_files, time_diffs

#################################################################################################################
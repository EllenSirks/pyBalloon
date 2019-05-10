import numpy as np
import elevation
import get_gfs
import pyb_aux
import pyb_io
import time

import match_files
import ratio

# method to read in model data
def read_data(loc0=None, weather_file=None, balloon=None, descent_only=False):

	time0 = time.time()

	lat0, lon0, alt0 = loc0

	in_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

	tile_size = 10. # degrees (read a tile this wide/high from the GFS grb2 file)
	area = (lat0 + (tile_size/2.), lon0 - (tile_size/2.), lat0 - (tile_size/2.), lon0 + (tile_size/2.)) # note top, left, bottom, right ordering for area

	model_data = pyb_io.read_gfs_single(directory=in_dir + weather_file, area=area, alt0=alt0, descent_only=descent_only, step=balloon['altitude_step'])[0]
	print('GFS data read, %.1f s elapsed' % (time.time() - time0))

	return model_data ## need [0] otherwise get a list of a dict. (prev. fixed by looping over the data)

# method to calculate new lat/lon coordinates from Cartesian displacements. 

# Required arguments: lat_rad [rad] (current latitude), lon_rad [rad] (current longitude), alt [m] (current altitude), dx (East - west movement), dy [m] (North - south movement)
# Note: for dx, east is positive and for dy, north is positive

def movement2ll(lat_rad=None, lon_rad=None, alt=None, dx=None, dy=None):

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

	# return new coordinates: latitude [radians], longitude [radians], and distance traveled [km]
	return lat2, lon2, dist

# method to calculate necessary properties and carry out the interpolation
def calc_properties(data=None, weather_file=None, loc0=None, balloon=None, descent_only=False):

	if data == None:
		data = read_data(loc0=loc0, weather_file=weather_file)

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

	time0 = time.time()

	data = pyb_aux.data_interpolation(data, alt0, balloon['altitude_step'], mode='spline', descent_only=descent_only)

	print('Interpolation done, %.3f s elapsed' % (time.time() - time0) + '\n')

	return data

def calc_displacements(data=None, balloon=None, descent_only=False):

	data['descent_speeds'] = pyb_aux.descent_speed(data, balloon['equip_mass'], balloon['Cd_parachute'], balloon['parachute_areas'], balloon['altitude_step'], balloon['parachute_change_altitude'])

	if not descent_only:

		max_altitudes, max_alt_idxs = pyb_aux.burst_altitude(data, balloon['burst_radius'])

		delta_t = balloon['altitude_step'] / data['ascent_speeds']

		data['ascent_time_steps'] = delta_t
		data['cumulative_ascent_times'] = np.cumsum(delta_t)/60.

		# calculate movement in x/y direction for each time step
		data['dxs_up'] = data['u_winds'] * data['ascent_time_steps']
		data['dys_up'] = data['v_winds'] * data['ascent_time_steps']

	delta_t = -1*balloon['altitude_step'] / data['descent_speeds']
	data['descent_time_steps'] = delta_t
	data['cumulative_descent_times'] = np.cumsum(delta_t)/60

	# calculate movement in x/y direction for each time step
	data['dxs_down'] = data['u_winds'] * data['descent_time_steps']
	data['dys_down'] = data['v_winds'] * data['descent_time_steps']

	return data

# method to calculate movement of balloon until it has reached the ground
def calc_movements(data=None, loc0=None, datestr=None, utc_hour=None, balloon=None, alt_change_fit=10, descent_only=False, interpolate=False, drift_time=0):

	time0 = time.time()
	lat0, lon0, alt0 = loc0

	keys = list(data.keys())
	checks = [0. for key in keys]
	index = 0

	# set initial conditions
	alts = data[keys[0]]['altitudes'] # alts, lats and lons are the same for all weather files (if we dont give it different areas)
	data_lats = np.radians(data[keys[0]]['lats'])
	data_lons = np.radians(data[keys[0]]['lons'])

	lat_rad = [np.radians(lat0)]
	lon_rad = [np.radians(lon0)]
	all_alts = [alt0]
	total_time = [0]
	dists = [0]
	distance_travelled = 0
	descent_speeds = []

	# Calculate the movement during ascent
	i = 0
	print('Calculating ascent...')
	if not descent_only:

		while True:

			# Find the closest grid point
			diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons - lon_rad[-1])**2)
			grid_i, = np.where(diff == diff.min())
			grid_i = grid_i[0]

			current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600

			for j in range(len(keys)):

				key = keys[j]
				if int(current_time) == key and current_time >= key and checks[j] == 0.:

					new_hhhh, new_hhh = get_gfs.get_closest_hr(utc_hour=key)
					new_weather_file = 'gfs_4_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh).zfill(3) + '.grb2'
					new_weather_file = get_gfs.get_gfs_file(weather_file=new_weather_file)[0]
					loc = (np.degrees(lat_rad[-1]), np.degrees(lon_rad[-1]), all_alts[-1])
					data[new_hhhh + new_hhh] = prepare_data(weather_file=new_weather_file, loc0=loc0, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only)
					keys = list(data.keys())
					checks[j] += 1.

			# add new weather file if current time is past 'latest' weather file
			if current_time > max(keys) + 1.5:

				new_hhhh, new_hhh = get_gfs.get_closest_hr(utc_hour=current_time)
				new_weather_file = 'gfs_4_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh).zfill(3) + '.grb2'
				new_weather_file = get_gfs.get_gfs_file(weather_file=new_weather_file)[0]
				loc = (np.degrees(lat_rad[-1]), np.degrees(lon_rad[-1]), all_alts[-1])
				data[new_hhhh + new_hhh] = prepare_data(weather_file=new_weather_file, loc0=loc0, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only)
				keys = list(data.keys())

				if new_hhh == 0.:
					checks.append(1.)
				else:
					checks.append(0.)

				index += 1

			# calculate change in latitude & longitude, and distance travelled
			if interpolate:

				(t1, f1), (t2, f2) = calc_time_frac(current_time = current_time, weather_times=keys)

				dx = f1*data[t1]['dxs_up'][grid_i, i] + f2*data[t2]['dxs_up'][grid_i, i]
				dy = f1*data[t1]['dys_up'][grid_i, i] + f2*data[t2]['dys_up'][grid_i, i]
				dt = f1*data[t1]['ascent_time_steps'][grid_i, i] + f2*data[t2]['ascent_time_steps'][grid_i, i]

			else:

				dx = data[keys[index]]['dxs_up'][grid_i, i]
				dy = data[keys[index]]['dys_up'][grid_i, i]
				dt = data[keys[index]]['ascent_time_steps'][grid_i, i]

			total_time.append(dt)

			lat, lon, dist = movement2ll(lat_rad=lat_rad[-1], lon_rad=lon_rad[-1], alt=alts[i], dx=dx, dy=dy)

			lat_rad.append(lat)
			lon_rad.append(lon)
			all_alts.append(alts[i])
			distance_travelled += dist

			if all_alts[-1] >= max_altitudes[grid_i]:
				break

			i += 1

	else: 

		# Mimick the movement during ascent if we only want descent
		while True:

			if alts[i] >= alt0:
				break

			i += 1

		i -= 1

	# add first descent speed (needs to be outside of the loop atm.)
	current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600
	

	diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons-lon_rad[-1])**2)
	grid_i, = np.where(diff == diff.min())
	grid_i = grid_i[0]

	if interpolate:
		(t1, f1), (t2, f2) = calc_time_frac(current_time = current_time, weather_times=keys)
		descent_speed = f1*data[t1]['descent_speeds'][grid_i, i+1] + f2*data[t2]['descent_speeds'][grid_i, i+1]
	else:
		descent_speed = data[keys[index]]['descent_speeds'][grid_i, i+1]
	descent_speeds.append(descent_speed)

	# option for drifting the balloon before descent
	timer, dt = 0, 1 # seconds
	while timer < 60*drift_time: # seconds (drift time is in min.)

		if timer == 0:
			print('Calculating drift trajectory...')

		diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons-lon_rad[-1])**2)
		grid_i, = np.where(diff == diff.min())
		grid_i = grid_i[0]

		current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600
		
		for j in range(len(keys)):

			key = keys[j]
			if int(current_time) == key and current_time >= key and checks[j] == 0.:

				new_hhhh, new_hhh = get_gfs.get_closest_hr(utc_hour=key)
				new_weather_file = 'gfs_4_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh).zfill(3) + '.grb2'
				new_weather_file = get_gfs.get_gfs_file(weather_file=new_weather_file)[0]
				loc = (np.degrees(lat_rad[-1]), np.degrees(lon_rad[-1]), all_alts[-1])
				data[new_hhhh + new_hhh] = prepare_data(weather_file=new_weather_file, loc0=loc0, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only)
				keys = list(data.keys())
				checks[j] += 1.

		# add new weather file if current time is past 'latest' weather file
		if current_time > max(keys) + 1.5:

			new_hhhh, new_hhh = get_gfs.get_closest_hr(utc_hour=current_time)
			new_weather_file = 'gfs_4_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh).zfill(3) + '.grb2'
			new_weather_file = get_gfs.get_gfs_file(weather_file=new_weather_file)[0]
			loc = (np.degrees(lat_rad[-1]), np.degrees(lon_rad[-1]), all_alts[-1])
			data[new_hhhh + new_hhh] = prepare_data(weather_file=new_weather_file, loc0=loc0, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only)
			keys = list(data.keys())

			if new_hhh == 0.:
				checks.append(1.)
			else:
				checks.append(0.)

			index += 1

		# calculate change in latitude & longitude, and distance travelled
		if interpolate:

			(t1, f1), (t2, f2) = calc_time_frac(current_time = current_time, weather_times=keys)

			dx = (f1*data[t1]['u_winds'][grid_i, i] + f2*data[t2]['u_winds'][grid_i, i])*dt
			dy = (f1*data[t1]['v_winds'][grid_i, i] + f2*data[t2]['v_winds'][grid_i, i])*dt
			descent_speed = f1*data[t1]['descent_speeds'][grid_i, i] + f2*data[t2]['descent_speeds'][grid_i, i]

		else:

			dx = data[keys[index]]['u_winds'][grid_i, i]*dt
			dy = data[keys[index]]['v_winds'][grid_i, i]*dt
			descent_speed = data[keys[index]]['descent_speeds'][grid_i, i]

		total_time.append(dt)
		descent_speeds.append(descent_speed)

		lat, lon, dist = movement2ll(lat_rad=lat_rad[-1], lon_rad=lon_rad[-1], alt=all_alts[-1], dx=dx, dy=dy)

		all_alts.append(all_alts[-1])
		lat_rad.append(lat)
		lon_rad.append(lon)
		dists.append(dist)
		distance_travelled += dist

		timer += dt

	print('Calculating descent...')
	# Calculate the movement during descent (same for either option)
	while i >= 0:

		# Find the closest grid point
		diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons-lon_rad[-1])**2)
		grid_i, = np.where(diff == diff.min())
		grid_i = grid_i[0]

		current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600

		# update current weather files if current time is past latest data release (i.e. 0600_006 -> 1200_000)
		for j in range(len(keys)):

			key = keys[j]
			if int(current_time) == key and current_time >= key and checks[j] == 0.:

				new_hhhh, new_hhh = get_gfs.get_closest_hr(utc_hour=key)
				new_weather_file = 'gfs_4_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh).zfill(3) + '.grb2'
				new_weather_file = get_gfs.get_gfs_file(weather_file=new_weather_file)[0]
				loc = (np.degrees(lat_rad[-1]), np.degrees(lon_rad[-1]), all_alts[-1])
				data[new_hhhh + new_hhh] = prepare_data(weather_file=new_weather_file, loc0=loc0, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only)
				keys = list(data.keys())
				checks[j] += 1.

		# add new weather file if current time is past 'latest' weather file
		if current_time > max(keys) + 1.5:

			new_hhhh, new_hhh = get_gfs.get_closest_hr(utc_hour=current_time)
			new_weather_file = 'gfs_4_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh).zfill(3) + '.grb2'
			new_weather_file = get_gfs.get_gfs_file(weather_file=new_weather_file)[0]
			loc = (np.degrees(lat_rad[-1]), np.degrees(lon_rad[-1]), all_alts[-1])
			data[new_hhhh + new_hhh] = prepare_data(weather_file=new_weather_file, loc0=loc0, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only)
			keys = list(data.keys())

			if new_hhh == 0.:
				checks.append(1.)
			else:
				checks.append(0.)

			index += 1

		# calculate change in latitude & longitude, and distance travelled
		if interpolate:

			(t1, f1), (t2, f2) = calc_time_frac(current_time = current_time, weather_times=keys)

			dx = f1*data[t1]['dxs_down'][grid_i, i] + f2*data[t2]['dxs_down'][grid_i, i]
			dy = f1*data[t1]['dys_down'][grid_i, i] + f2*data[t2]['dys_down'][grid_i, i]
			dt = f1*data[t1]['descent_time_steps'][grid_i, i] + f2*data[t2]['descent_time_steps'][grid_i, i]
			descent_speed = f1*data[t1]['descent_speeds'][grid_i, i] + f2*data[t2]['descent_speeds'][grid_i, i]

		else:

			dx = data[keys[index]]['dxs_down'][grid_i, i]
			dy = data[keys[index]]['dys_down'][grid_i, i]
			dt = data[keys[index]]['descent_time_steps'][grid_i, i]
			descent_speed = data[keys[index]]['descent_speeds'][grid_i, i]

		total_time.append(dt)
		descent_speeds.append(descent_speed)

		lat, lon, dist = movement2ll(lat_rad=lat_rad[-1], lon_rad=lon_rad[-1], alt=alts[i], dx=dx, dy=dy)

		lat_rad.append(lat)
		lon_rad.append(lon)
		all_alts.append(alts[i])
		dists.append(dist)
		distance_travelled += dist

		i -= 1

	# Convert the result array lists to Numpy 2D-arrays
	output = {}
	output['lats'] = np.degrees(np.array(lat_rad)) # to decimal degrees
	output['lons'] = np.degrees(np.array(lon_rad)) # to decimal degrees
	output['alts'] = np.array(all_alts)
	output['dists'] = np.array(dists)
	output['descent_speeds'] = np.array(descent_speeds)
	output['times'] = np.cumsum(np.array(total_time))/60 # to minutes
	output['distance'] = distance_travelled

	# print out relevant quantities
	print('\nTrajectories calculated, %.3f s elapsed' % (time.time() - time0) + '\n')
	print('Maximum altitude: ' + str(np.max(all_alts)) + ' m')
	print('Landing location: (%.6f, %.6f)' % (output['lats'][-1], output['lons'][-1]))
	print('Total flight time: %d min' % (int(output['times'][-1])) + ', total distance travelled: %.1f' % distance_travelled + ' km')

	return output

def prepare_data(weather_file=None, loc0=None, utc_hour=None, balloon=None, descent_only=False, drift_time=0):

	model_data = read_data(loc0=loc0, weather_file=weather_file, balloon=balloon, descent_only=descent_only)
	model_data = calc_properties(data=model_data, weather_file=weather_file, loc0=loc0, balloon=balloon, descent_only=descent_only)
	model_data = calc_displacements(data=model_data, balloon=balloon, descent_only=descent_only)

	return model_data

def run_traj(weather_files=None, loc0=None, datestr=None, utc_hour=None, balloon=None, descent_only=False, interpolate=False, drift_time=0):

	time0 = time.time()

	data_array = {}
	for weather_file in weather_files:

		model_data = prepare_data(weather_file=weather_file, loc0=loc0, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only)
		data_array[int(int(weather_file[15:19])/100.) + int(weather_file[20:23])] = model_data

	trajectories = calc_movements(data=data_array, loc0=loc0, datestr=datestr, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only, interpolate=interpolate, drift_time=drift_time)

	return trajectories

def calc_time_frac(current_time=None, weather_times=None, weather_files=None):

	if weather_times == None:
		weather_times = []
		for file in weather_files:
			weather_times.append(int(int(file[15:19])/100.) + int(file[20:23]))

	for time in weather_times:
		if time < current_time:
			earlier_time = time
		elif time > current_time:
			later_time = time
			break

	dt1 = current_time - earlier_time
	dt2 = later_time - current_time

	dt_total = later_time - earlier_time

	frac1 = 1. - dt1/dt_total
	frac2 = 1. - dt2/dt_total

	return (earlier_time, frac1), (later_time, frac2)

if __name__ == '__main__':

	(t1, f1), (t2, f2) = calc_time_frac(current_time=11, weather_files= ['gfs_4_20180304_0600_003.grb2', 'gfs_4_20180304_0600_006.grb2', 'gfs_4_20180304_0600_009.grb2'])
	(t1, f1), (t2, f2) = calc_time_frac(current_time=11, weather_times=[9, 12, 15])
import numpy as np
import elevation
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

	tile_size = 2. # degrees (read a tile this wide/high from the GFS grb2 file)
	area = (lat0 + (tile_size/2.), lon0 - (tile_size/2.), lat0 - (tile_size/2.), lon0 + (tile_size/2.)) # note top, left, bottom, right ordering for area

	model_data = pyb_io.read_gfs_single(in_dir + weather_file, area=area, alt0=alt0, descent_only=descent_only, step = balloon['altitude_step'])
	print('GFS data read, %.1f s elapsed' % (time.time() - time0))

	return model_data[0] ## need [0] otherwise get a list of a dict. (prev. fixed by looping over the data)

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

	# rho_pre = pyb_aux.air_density(data)
	# u_pre = data['u_winds'].copy()
	# v_pre = data['v_winds'].copy()
	# T_pre = data['temperatures'].copy()

	if not descent_only:

		data['balloon_radii'], gas_mass = pyb_aux.mooney_rivlin(data, balloon['radius_empty'], balloon['fill_radius'], balloon['thickness_empty'])
		data['balloon_volumes'] = pyb_aux.balloon_volume(data)

		total_mass = balloon['equip_mass'] + balloon['balloon_mass'] + gas_mass # kg

		data['lifts'] = pyb_aux.lift(data, total_mass)
		data['ascent_speeds'] = pyb_aux.ascent_speed(data, total_mass, balloon['Cd_balloon'])

		if 'simple_ascent_rate' in balloon:
			data['ascent_speeds'] = np.ones_like(data['ascent_speeds']) * balloon['simple_ascent_rate']

	time0 = time.time()

	# alts_pre = data['altitudes'].copy()

	data = pyb_aux.data_interpolation(data, alt0, balloon['altitude_step'], mode='spline', descent_only=descent_only)

	print('Interpolation done, %.1f s elapsed' % (time.time() - time0))

	# return data, rho_pre, u_pre, v_pre, T_pre, alts_pre
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
# def calc_movements(data=None, weather_file=None, loc0=None, balloon=None, alt_change_fit=10, descent_only=False, interpolation_files=1, name=None, ext='all_points/'):
def calc_movements(data=None, loc0=None, balloon=None, alt_change_fit=10, descent_only=False):

	time0 = time.time()
	lat0, lon0, alt0 = loc0

	if type(data) == dict: 

	# calculate movement in x/y direction for each time step
	if not descent_only:

		dxs_up = data['dxs_up']
		dys_up = data['dys_up']

	dxs_down = data['dxs_down']
	dys_down = data['dys_down']

	# set initial conditions
	alts = data['altitudes']
	data_lats = np.radians(data['lats'])
	data_lons = np.radians(data['lons'])
			
	lat_rad = [np.radians(lat0)]
	lon_rad = [np.radians(lon0)]
	all_alts = [alt0]
	total_time = [0]
	dists = [0]
	distance_travelled = 0

	i = 0
	if not descent_only:

		# Calculate the movement during ascent
		while True:

			# Find the closest grid point
			diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons - lon_rad[-1])**2)
			grid_i, = np.where(diff == diff.min())
			grid_i = grid_i[0]

			# calculate change in latitude & longitude, and distance travelled
			lat, lon, dist = movement2ll(lat_rad=lat_rad[-1], lon_rad=lon_rad[-1], alt=alts[i], dx=dxs_up[grid_i, i], dy=dys_up[grid_i, i])
			lat_rad.append(lat)
			lon_rad.append(lon)
			total_time.append(data['ascent_time_steps'][grid_i, i])
			all_alts.append(alts[i])
			distance_travelled += dist

			if all_alts[-1] >= max_altitudes[grid_i]:
				break

			i += 1

	else: 

		# Mimick the movement during ascent
		while alts[i] < alt0:
			i += 1

		i -= 1

	ind = i

	diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons-lon_rad[-1])**2)
	grid_i, = np.where(diff == diff.min())
	grid_i = grid_i[0]

	descent_speeds = [data['descent_speeds'][grid_i, ind + 1]]

	# f = open('/home/ellen/Desktop/SuperBIT/properties/' + ext + 'prop_afterinterp_' + name + '.dat', 'w+')
	# f.write('alt0 grid alt rho u v T\n')

	# Calculate the movement during descent (same for either option)

	# grids = []
	while i >= 0:

		# Find the closest grid point

		diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons-lon_rad[-1])**2)
		grid_i, = np.where(diff == diff.min())
		grid_i = grid_i[0]

		descent_speeds.append(data['descent_speeds'][grid_i, i])
		# f.write(str(alt0) + ' ' + str(grid_i) + ' '  + str(data['altitudes'][i]) + ' ' + str(data['air_densities'][grid_i, i]) + ' ' + str(data['u_winds'][grid_i, i]) + ' ' + str(data['v_winds'][grid_i, i]) + ' ' + str(data['temperatures'][grid_i, i]) + '\n')

		# calculate change in latitude & longitude, and distance travelled

		lat, lon, dist = movement2ll(lat_rad=lat_rad[-1], lon_rad=lon_rad[-1], alt=alts[i], dx=dxs_down[grid_i, i], dy=dys_down[grid_i, i])
		# grids.append(grid_i)

		lat_rad.append(lat)
		lon_rad.append(lon)
		total_time.append(data['descent_time_steps'][grid_i, i])
		all_alts.append(alts[i])
		dists.append(dist)
		distance_travelled += dist

		i -= 1

	# f.close()

	# f = open('/home/ellen/Desktop/SuperBIT/properties/' + ext + 'prop_preinterp_' + name + '.dat', 'w+')
	# f.write('alt rho u v T\n')
	# for i in range(len(rho_pre)):
	#	 f.write(str(alts_pre[i][int(np.mean(grids))]) + ' ' + str(rho_pre[i][int(np.mean(grids))]) + ' ' + str(u_pre[i][int(np.mean(grids))]) + ' ' +str(v_pre[i][int(np.mean(grids))]) + ' ' + str(T_pre[i][int(np.mean(grids))]) + '\n')
	# f.close()

	# Convert the result array lists to Numpy 2D-arrays
	output = {}
	output['lats'] = np.degrees(np.array(lat_rad)) # to decimal degrees
	output['lons'] = np.degrees(np.array(lon_rad)) # to decimal degrees
	output['alts'] = np.array(all_alts)
	output['dists'] = np.array(dists)
	output['descent_speeds'] = np.array(descent_speeds)
	output['times'] = np.cumsum(np.array(total_time))/60 # to minutes
	output['distance'] = distance_travelled

	print('Trajectories calculated, %.1f s elapsed' % (time.time() - time0) + '\n')

	print('Maximum altitude: ' + str(np.max(all_alts)) + ' m')
	print('Landing location: (%.6f, %.6f)' % (output['lats'][-1], output['lons'][-1]))
	print('Total time: %d min of flight' % (int(output['times'][-1])) + ', total distance travelled: %.1f' % distance_travelled + ' km')

	return output

def single_weatherfile(weather_file=None, loc0=None, balloon=None, descent_only=False):

	time0 = time.time()

	model_data = read_data(loc0=loc0, weather_file=weather_file, balloon=balloon, descent_only=descent_only)
	model_data = calc_properties(data=model_data, weather_file=weather_file, loc0=loc0, balloon=balloon, descent_only=descent_only)
	model_data = calc_displacements(data=model_data, balloon=balloon, descent_only=descent_only)

	trajectories = calc_movements(data=model_data, loc0=loc0, balloon=balloon, descent_only=descent_only)

	print('Trajectories calculated, %.1f s elapsed' % (time.time() - time0))

	return trajectories

def multiple_weatherfiles(weather_files=None, loc0=None, balloon=None, descent_only=False):

	time0 = time.time()

	data_array = {}
	for weather_file in weather_files:

		model_data = read_data(loc0=loc0, weather_file=weather_file, balloon=balloon, descent_only=descent_only)
		model_data = calc_properties(data=model_data, weather_file=weather_file, loc0=loc0, balloon=balloon, descent_only=descent_only)

		data_array[int(int(weather_file[15:19])/100.) + int(weather_file[20:23])] = model_data

	print('Trajectories calculated, %.1f s elapsed' % (time.time() - time0))

	# return trajectories

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

	frac1 = dt1/dt_total
	frac2 = dt2/dt_total

	return (earlier_time, frac1), (later_time, frac2)

if __name__ == '__main__':

	(t1, f1), (t2, f2) = calc_time_frac(current_time=11, weather_files= ['gfs_4_20180304_0600_003.grb2', 'gfs_4_20180304_0600_006.grb2', 'gfs_4_20180304_0600_009.grb2'])

	print(t1, f1, t2, f2)

	(t1, f1), (t2, f2) = calc_time_frac(current_time=11, weather_times=[9, 12, 15])

	print(t1, f1, t2, f2)

# method to calculate movement of balloon until it has reached the ground
def calc_movements_2(data=None, loc0=None, datestr=None, utc_hour=None, balloon=None, alt_change_fit=10, descent_only=False, interpolate=False, drift_time=0):

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
	while True:

		if not descent_only:

			# Find the closest grid point
			diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons - lon_rad[-1])**2)
			grid_i, = np.where(diff == diff.min())
			grid_i = grid_i[0]

			data, keys, checks, index = update_files(data=data, lat_rad=lat_rad, lon_rad=lon_rad, all_alts=all_alts, balloon=balloon, datestr=datestr, utc_hour=utc_hour, \
				loc0=loc0, total_time=total_time, checks=checks, index=index, descent_only=descent_only)

			# calculate change in latitude & longitude, and distance travelled
			if interpolate:

				current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600
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
			if alts[i] >= alt0:
				break
			i += 1

	i -= 1

	# add first descent speed (needs to be outside of the loop atm.)
	diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons-lon_rad[-1])**2)
	grid_i, = np.where(diff == diff.min())
	grid_i = grid_i[0]

	if interpolate:
		current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600
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

		data, keys, checks, index = update_files(data=data, lat_rad=lat_rad, lon_rad=lon_rad, all_alts=all_alts, balloon=balloon, datestr=datestr, utc_hour=utc_hour, \
			loc0=loc0, total_time=total_time, checks=checks, index=index, descent_only=descent_only)

		# calculate change in latitude & longitude, and distance travelled
		if interpolate:

			current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600
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

		data, keys, checks, index = update_files(data=data, lat_rad=lat_rad, lon_rad=lon_rad, all_alts=all_alts, balloon=balloon, datestr=datestr, utc_hour=utc_hour, \
			loc0=loc0, total_time=total_time, checks=checks, index=index, descent_only=descent_only)

		# calculate change in latitude & longitude, and distance travelled
		if interpolate:

			current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600
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
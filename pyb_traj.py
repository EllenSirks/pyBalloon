import numpy as np
import time

import get_gfs
import pyb_aux
import pyb_io
import test

# method to calculate the fractions needed for interpolation
# e.g. if time is 10 then on the lhs we have 6 and on the rhs we have 12, 10 is then 2/3 of the way away from 6 and 1/3 away from 12.
# the fractions needed are then 1-2/3 for 6 and 1 - 1/3 for 12
def calc_time_frac(current_time=None, weather_times=None, weather_files=None):

	if weather_times == None:
		weather_times = []
		for file in weather_files:
			weather_times.append(int(int(file[15:19])/100.) + int(file[20:23]))

	weather_times.sort()

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

# method to read in model data
def read_data(loc0=None, weather_file=None, balloon=None, descent_only=False):

	time0 = time.time()

	lat0, lon0, alt0 = loc0

	in_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/GFS/'

	tile_size = 10. # degrees (read a tile this wide/high from the GFS grb2 file)
	area = (lat0 + (tile_size/2.), lon0 - (tile_size/2.), lat0 - (tile_size/2.), lon0 + (tile_size/2.)) # note top, left, bottom, right ordering for area

	model_data = pyb_io.read_gfs_single(directory=in_dir + weather_file, area=area, alt0=alt0, descent_only=descent_only, step=balloon['altitude_step'])[0]
	print('GFS data read, %.1f s elapsed' % (time.time() - time0))

	return model_data

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

	# data['descent_speeds'] = pyb_aux.descent_speed(data, balloon['equip_mass']+balloon['balloon_mass'], balloon['Cd_parachute'], balloon['parachute_areas'], balloon['altitude_step'], balloon['parachute_change_altitude'])
	data['descent_speeds'] = pyb_aux.descent_speed(data, balloon['equip_mass'], balloon['Cd_parachute'], balloon['parachute_areas'], balloon['altitude_step'], balloon['parachute_change_altitude']) # check this!

	if not descent_only:

		data['max_altitudes'], data['max_alt_idxs'] = pyb_aux.burst_altitude(data, balloon['burst_radius'])
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

# method to prepare the data for the calculations of the trajectory
def prepare_data(weather_file=None, loc0=None, current_time=None, balloon=None, descent_only=False, drift_time=0):

	model_data1 = read_data(loc0=loc0, weather_file=weather_file, balloon=balloon, descent_only=descent_only)
	err_data = pyb_aux.calc_uv_errs(weather_file=weather_file, loc0=loc0, current_time=current_time, descent_only=descent_only, main_data=model_data1)
	model_data2 = calc_properties(data=model_data1, weather_file=weather_file, loc0=loc0, balloon=balloon, descent_only=descent_only)
	model_data2 = pyb_aux.add_uv_errs(main_data=model_data2, err_data=err_data)
	model_data3 = calc_displacements(data=model_data2, balloon=balloon, descent_only=descent_only)

	return model_data3

# method to pull everything together and return the trajectory
def run_traj(weather_files=None, loc0=None, datestr=None, utc_hour=None, balloon=None, descent_only=False, interpolate=False, drift_time=0):

	time0 = time.time()

	data_array = {}
	for weather_file in weather_files:

		model_data = prepare_data(weather_file=weather_file, loc0=loc0, current_time=utc_hour, balloon=balloon, descent_only=descent_only)
		data_array[int(int(weather_file[15:19])/100.) + int(weather_file[20:23])] = model_data

	trajectories = calc_movements(data=data_array, loc0=loc0, datestr=datestr, utc_hour=utc_hour, balloon=balloon, descent_only=descent_only, interpolate=interpolate, drift_time=drift_time)

	return trajectories

# method to update weather files to newest available
def update_files(data, lat_rad, lon_rad, all_alts, balloon, datestr, utc_hour, loc0, total_time, checks, index, descent_only):

	keys = list(data.keys())
	current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600
			
	for j in range(len(keys)):

		key = keys[j]
		if int(current_time) == key and current_time >= key and checks[j] == 0.:

			new_hhhh, new_hhh1, new_hhh2 = get_gfs.get_closest_hr(utc_hour=key)
			new_weather_file = 'gfs_4_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh1).zfill(3) + '.grb2'
			new_weather_file = test.get_gfs_file(weather_files=[new_weather_file])[0]
			loc = (np.degrees(lat_rad[-1]), np.degrees(lon_rad[-1]), all_alts[-1])
			data[new_hhhh + new_hhh1] = prepare_data(weather_file=new_weather_file, loc0=loc0, current_time=current_time, balloon=balloon, descent_only=descent_only)
			keys = list(data.keys())
			checks[j] += 1.

	# add new weather file if current time is past 'latest' weather file
	if current_time > max(keys) + 1.5:

		new_hhhh, new_hhh1, new_hhh2 = get_gfs.get_closest_hr(utc_hour=current_time)
		new_weather_file = 'gfs_4_' + datestr + '_' + str(new_hhhh*100).zfill(4) + '_' + str(new_hhh1).zfill(3) + '.grb2'
		new_weather_file = test.get_gfs_file(weather_files=[new_weather_file])[0]
		loc = (np.degrees(lat_rad[-1]), np.degrees(lon_rad[-1]), all_alts[-1])
		data[new_hhhh + new_hhh1] = prepare_data(weather_file=new_weather_file, loc0=loc0, current_time=current_time, balloon=balloon, descent_only=descent_only)
		keys = list(data.keys())

		if new_hhh1 == 0.:
			checks.append(1.)
		else:
			checks.append(0.)

		index += 1

	return data, keys, checks, index

# method to calculate the trajectory of the balloon/parachute given a start position & other input parameters
def calc_movements(data=None, loc0=None, datestr=None, utc_hour=None, balloon=None, alt_change_fit=10, descent_only=False, interpolate=False, drift_time=0):

	time0 = time.time()
	
	lat0, lon0, alt0 = loc0

	keys = list(data.keys())
	checks = [0. for key in keys]

	hrs = np.array([0., 6., 12., 18.])
	if not interpolate:
		if int(float(utc_hour)) in hrs:
			hr = hrs[np.where(hrs == int(float(utc_hour)))[0][0]]
			if float(utc_hour) >= hr:
				checks[0] = 1.

	index, i, timer = 0, 0, 0
	stage = 1

	# set initial conditions
	alts = data[keys[0]]['altitudes'] # alts, lats and lons are the same for all weather files (if we dont give it different areas)
	data_lats, data_lons  = np.radians(data[keys[0]]['lats']), np.radians(data[keys[0]]['lons'])
	data_lats_err, data_lons_err = np.radians(data[keys[0]]['lats_err']), np.radians(data[keys[0]]['lons_err'])

	lat_rad, lon_rad, all_alts = [np.radians(lat0)], [np.radians(lon0)], [alt0]

	if lon_rad[0] < 0:
		lon_rad[0] += 2*np.pi
	if lat_rad[0] < 0:
		lat_rad[0] += 2*np.pi

	total_time, dists, dists_u, dists_v = [0], [0], [0], [0]
	speeds, temperatures = [], []
	sigmas_u, sigmas_v = [], []

	initial_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600
	if interpolate:
		(t1_ini, f1_ini), (t2_ini, f2_ini) = calc_time_frac(current_time = initial_time, weather_times=keys)

	while True:

		data, keys, checks, index = update_files(data=data, lat_rad=lat_rad, lon_rad=lon_rad, all_alts=all_alts, balloon=balloon, datestr=datestr, utc_hour=utc_hour, \
			loc0=loc0, total_time=total_time, checks=checks, index=index, descent_only=descent_only)

		# Find the closest grid point
		diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons - lon_rad[-1])**2)
		grid_i, = np.where(diff == diff.min())
		grid_i = grid_i[0]

		diff = np.sqrt((data_lats_err - lat_rad[-1])**2 + (data_lons_err - lon_rad[-1])**2)
		grid_i_err, = np.where(diff == diff.min())
		grid_i_err = grid_i_err[0]

		if interpolate:

			current_time = float(utc_hour) + np.cumsum(np.array(total_time))[-1]/3600
			(t1, f1), (t2, f2) = calc_time_frac(current_time = current_time, weather_times=keys)

		if stage == 1:
			if i == 0:
				print('Calculating ascent...')

			if descent_only:

				# Mimick the movement during ascent if we only want descent
				if alts[i] >= alt0:
					if drift_time == 0:
						stage += 2
					else:
						stage += 1
					final_i = i
					continue

			else:

				if interpolate:

					dx = f1*data[t1]['dxs_up'][grid_i, i] + f2*data[t2]['dxs_up'][grid_i, i]
					dy = f1*data[t1]['dys_up'][grid_i, i] + f2*data[t2]['dys_up'][grid_i, i]
					dt = f1*data[t1]['ascent_time_steps'][grid_i, i] + f2*data[t2]['ascent_time_steps'][grid_i, i]
					speed = f1*data[t1]['ascent_speeds'][grid_i, i] + f2*data[t2]['ascent_speeds'][grid_i, i]
					T = f1*data[t1]['temperatures'][grid_i, i] + f2*data[t2]['temperatures'][grid_i, i]
					sigma_u = (f1*data[t1]['u_wind_errs'][grid_i_err, i] + f2*data[t2]['u_wind_errs'][grid_i_err, i])*dt
					sigma_v = (f1*data[t1]['v_wind_errs'][grid_i_err, i] + f2*data[t2]['v_wind_errs'][grid_i_err, i])*dt

				else:

					dx = data[keys[index]]['dxs_up'][grid_i, i]
					dy = data[keys[index]]['dys_up'][grid_i, i]
					dt = data[keys[index]]['ascent_time_steps'][grid_i, i]
					speed = data[keys[index]]['ascent_speeds'][grid_i, i]
					T = data[keys[index]]['temperatures'][grid_i, i]
					sigma_u = (data[keys[index]]['u_wind_errs'][grid_i_err, i])*dt
					sigma_v = (data[keys[index]]['v_wind_errs'][grid_i_err, i])*dt

				if all_alts[-1] >= data[keys[index]]['max_altitudes'][grid_i]:
					if drift_time == 0:
						stage += 2
					else:
						stage += 1
					final_i = i
					continue

			i += 1

		if stage == 2:

			if timer == 0:
				print('Calculating drift trajectory...')

			dt = 5 # seconds

			if timer >= drift_time*60:
				stage += 1
				continue

			# calculate change in latitude & longitude, and distance travelled
			if interpolate:

				dx = (f1*data[t1]['u_winds'][grid_i, i] + f2*data[t2]['u_winds'][grid_i, i])*dt
				dy = (f1*data[t1]['v_winds'][grid_i, i] + f2*data[t2]['v_winds'][grid_i, i])*dt
				T = f1*data[t1]['temperatures'][grid_i, i] + f2*data[t2]['temperatures'][grid_i, i]
				sigma_u = (f1*data[t1]['u_wind_errs'][grid_i_err, i] + f2*data[t2]['u_wind_errs'][grid_i_err, i])*dt
				sigma_v = (f1*data[t1]['v_wind_errs'][grid_i_err, i] + f2*data[t2]['v_wind_errs'][grid_i_err, i])*dt

			else:

				dx = data[keys[index]]['u_winds'][grid_i, i]*dt
				dy = data[keys[index]]['v_winds'][grid_i, i]*dt
				T = data[keys[index]]['temperatures'][grid_i, i]
				sigma_u = (data[keys[index]]['u_wind_errs'][grid_i_err, i])*dt
				sigma_v = (data[keys[index]]['v_wind_errs'][grid_i_err, i])*dt

			speed = 0
			timer += dt

		if stage == 3:

			if i == final_i:
				print('Calculating descent...\n')

			# calculate change in latitude & longitude, and distance travelled
			if interpolate:

				dx = f1*data[t1]['dxs_down'][grid_i, i] + f2*data[t2]['dxs_down'][grid_i, i]
				dy = f1*data[t1]['dys_down'][grid_i, i] + f2*data[t2]['dys_down'][grid_i, i]
				dt = f1*data[t1]['descent_time_steps'][grid_i, i] + f2*data[t2]['descent_time_steps'][grid_i, i]
				speed = f1*data[t1]['descent_speeds'][grid_i, i] + f2*data[t2]['descent_speeds'][grid_i, i]
				T = f1*data[t1]['temperatures'][grid_i, i] + f2*data[t2]['temperatures'][grid_i, i]
				sigma_u = (f1*data[t1]['u_wind_errs'][grid_i_err, i] + f2*data[t2]['u_wind_errs'][grid_i_err, i])*dt
				sigma_v = (f1*data[t1]['v_wind_errs'][grid_i_err, i] + f2*data[t2]['v_wind_errs'][grid_i_err, i])*dt

			else:

				dx = data[keys[index]]['dxs_down'][grid_i, i]
				dy = data[keys[index]]['dys_down'][grid_i, i]
				dt = data[keys[index]]['descent_time_steps'][grid_i, i]
				speed = data[keys[index]]['descent_speeds'][grid_i, i]
				T = data[keys[index]]['temperatures'][grid_i, i]
				sigma_u = (data[keys[index]]['u_wind_errs'][grid_i_err, i])*dt
				sigma_v = (data[keys[index]]['v_wind_errs'][grid_i_err, i])*dt

			elevation = pyb_aux.get_elevation(lat=np.degrees(lat_rad[-1]), lon=np.degrees(lon_rad[-1]))
			if i == 0 or alts[i] <= elevation:
				break

			i -= 1

		if not (stage == 1 and descent_only):

			total_time.append(dt)

			if stage == 2:
				lat, lon, dist = movement2ll(lat_rad=lat_rad[-1], lon_rad=lon_rad[-1], alt=all_alts[-1], dx=dx, dy=dy)
				all_alts.append(all_alts[-1])
			else:
				lat, lon, dist = movement2ll(lat_rad=lat_rad[-1], lon_rad=lon_rad[-1], alt=alts[i], dx=dx, dy=dy)
				all_alts.append(alts[i])

			sigmas_u.append(sigma_u)
			sigmas_v.append(sigma_v)
			speeds.append(speed)
			temperatures.append(T)
			lat_rad.append(lat)
			lon_rad.append(lon)
			dists.append(dist)
			dists_u.append(np.sqrt(dx*dx))
			dists_v.append(np.sqrt(dy*dy))

	if interpolate:

		if descent_only:
			speed = f1_ini*data[t1]['descent_speeds'][grid_i, final_i] + f2_ini*data[t2_ini]['descent_speeds'][grid_i, final_i]
		else:
			speed = f1_ini*data[t1]['ascent_speeds'][grid_i, final_i] + f2_ini*data[t2_ini]['ascent_speeds'][grid_i, final_i]
		T = f1_ini*data[t1_ini]['temperatures'][grid_i, final_i] + f2_ini*data[t2_ini]['temperatures'][grid_i, final_i]

	else:

		if descent_only:
			speed = data[keys[index]]['descent_speeds'][grid_i, final_i]
		else:
			speed = data[keys[index]]['ascent_speeds'][grid_i, final_i] 
		T = data[keys[index]]['temperatures'][grid_i, final_i]

	speeds = [speed] + speeds
	temperatures = [T] + temperatures

	# get more accurate end-point based on elevation data
	new_end_point, new_alt = pyb_aux.get_endpoint(data=(np.degrees(np.array(lat_rad)), np.degrees(np.array(lon_rad)), np.array(all_alts), np.array(dists)))

	diffs = np.abs(np.degrees(np.array(lat_rad)) - new_end_point[0])
	diff = np.sqrt((np.degrees(np.array(lat_rad)) - new_end_point[0])**2 + (np.degrees(np.array(lon_rad)) - new_end_point[1])**2)
	index = np.where(diff == min(diff))[0][0]

	lat_rad, lon_rad, all_alts, dists, speeds, temperatures, total_time, sigmas_u, sigmas_v = \
		lat_rad[:index], lon_rad[:index], all_alts[:index], dists[:index+1], speeds[:index+1], temperatures[:index+1], total_time[:index+1], sigmas_u[:index+1], sigmas_v[:index+1]

	lat_rad.append(np.radians(new_end_point[0]))
	lon_rad.append(np.radians(new_end_point[1]))
	all_alts.append(new_alt)

	# Convert the result array lists to Numpy 2D-arrays
	output = {}
	output['lats'] = np.degrees(np.array(lat_rad)) # to decimal degrees
	output['lons'] = np.degrees(np.array(lon_rad)) # to decimal degrees
	output['alts'] = np.array(all_alts)
	output['dists'] = np.array(dists)
	output['speeds'] = np.array(speeds)
	output['temperatures'] = np.array(temperatures)
	output['times'] = np.cumsum(np.array(total_time))/60 # to minutes
	output['distance'] = np.sum(np.array(dists))
	output['distance_u'] = np.sum(np.array(dists_u))
	output['distance_v'] = np.sum(np.array(dists_v))
	output['sigmas_u'] = np.array(sigmas_u)
	output['sigmas_v'] = np.array(sigmas_v)
	output['error_forecast'] = (1./(output['distance']*1000))*np.sqrt((output['distance_u']*np.sqrt(np.sum(output['sigmas_u']**2)))**2 \
		+ (output['distance_v']*np.sqrt(np.sum(output['sigmas_u']**2)))**2)

	# print out relevant quantities
	print('Trajectories calculated, %.3f s elapsed' % (time.time() - time0) + '\n')
	print('Maximum altitude: ' + str(np.max(all_alts)) + ' m')
	print('Landing location: (%.6f, %.6f)' % (output['lats'][-1], output['lons'][-1]))
	print('Flight time: %d min' % (int(output['times'][-1])) + ', distance travelled: %.1f' % output['distance'] + ' km, error from forecast: ' + str(round(output['error_forecast'], 1)) + ' m')

	return output
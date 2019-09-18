"""Functions used by pyBalloon to calculate balloon trajectories"""

from haversine import haversine
import numpy as np
import datetime
import time
import sys

import get_gfs
import pyb_aux
import pyb_io

import param_file as p

#################################################################################################################

# method to calculate the fractions needed for interpolation
# e.g. if time is 10 then on the lhs we have 6 and on the rhs we have 12, 10 is then 2/3 of the way away from 6 and 1/3 away from 12.
# the fractions needed are then 1-2/3 for 6 and 1 - 1/3 for 12
def calc_time_frac(current_time=None, weather_times=None, weather_files=None):

	weather_files.sort()
	times = [(int(int(file[15:19])/100.) + int(file[20:23])) for file in weather_files]

	earlier_file = weather_files[-2]
	later_file = weather_files[-1]

	earlier_time = times[-2]
	later_time = times[-1]

	if later_time % 24 == 0:
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

# method to read in model data
def read_data(loc0=None, weather_file=None, descent_only=False):

	lat0, lon0, alt0 = loc0

	in_dir = p.path + p.weather_data_folder + p.GFS_folder

	tile_size = 60. # degrees (read a tile this wide/high from the GFS grb2 file)
	area = (lat0 + (tile_size/2.), lon0 - (tile_size/2.), lat0 - (tile_size/2.), lon0 + (tile_size/2.)) # note top, left, bottom, right ordering for area

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Reading GFS data from ' + str(weather_file) + '.grb2...'.ljust(20))
	sys.stdout.flush()

	model_data = pyb_io.read_gfs_single(directory=in_dir + weather_file, area=area, alt0=alt0, descent_only=descent_only)[0]

	return model_data

#################################################################################################################

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

#################################################################################################################

# method to calculate necessary properties and carry out the interpolation
def calc_properties(data=None, weather_file=None, loc0=None, balloon=None, descent_only=False, output_figs=False):

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

	data, figs_dict = pyb_aux.data_interpolation(data, alt0, balloon['altitude_step'], mode='spline', descent_only=descent_only, output_figs=output_figs)

	return data, figs_dict

#################################################################################################################

# method to calculate the displacements in the x/y directions
def calc_displacements(data=None, balloon=None, descent_only=False, vz_correct=False):

	g = p.g_0 * (p.R_e / (p.R_e + data['altitudes']))**2
	data['z_speeds'] = -data['omegas']/(data['air_densities']*g)

	data['descent_speeds'] = pyb_aux.descent_speed(data, balloon['equip_mass'], balloon['Cd_parachute'], balloon['parachute_areas'], balloon['altitude_step'], balloon['parachute_change_altitude'])

	if vz_correct:
		data['descent_speeds'] += data['z_speeds']

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

#################################################################################################################

# method to update weather files to newest available
def update_files(figs=None, used_weather_files=None, data=None, lat_rad=None, lon_rad=None, all_alts=None, current_time=None, balloon=None, datestr=None, utc_hour=None, loc0=None, \
	total_time=[0], time_diffs=[], index=0, params=None, output_figs=False):

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, check_sigmas, params, balloon = pyb_io.set_params(params=params, balloon=balloon)
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

	if (datestr == max_datestr and current_time > max_time and interpolate) or (not interpolate and datestr == max_datestr and current_time > max_time and (current_time - max_time) > 1.5) or (diff_days > 0 and hr_diff == 0):

		sys.stdout.write('\r')
		sys.stdout.flush()
		sys.stdout.write('Adding new weather file...'.ljust(60) + '\r')
		sys.stdout.flush()
		time.sleep(0.2)

		if not interpolate:

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
		 vz_correct=vz_correct, check_sigmas=check_sigmas, output_figs=output_figs)

		figs.append(figs_dict)
		used_weather_files[current_time] = new_weather_files
		keys = list(data.keys())
		index = np.where(np.array(keys) == new_weather_files[-1])[0][0]

	############################################################################################################

	weather_files = list(data.keys())
	weather_files.sort()

	# update from e.g. 0600_006 to 1200_000
	if not interpolate:
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
					vz_correct=vz_correct, check_sigmas=check_sigmas, output_figs=output_figs)

			figs.append(figs_dict)
			used_weather_files[current_time] = new_weather_file
			keys = list(data.keys())
			time_diffs[current_time] = [(new_hhhh + new_hhh1) % 24 - current_time]
			index = np.where(np.array(keys) == new_weather_file_name)[0][0]

	return data, keys, index, figs, used_weather_files, time_diffs

###############################################################################################################

# method to calculate the trajectory of the balloon/parachute given a start position & other input parameters
def calc_movements(data=None, used_weather_files=None, time_diffs=None, datestr=None, utc_hour=None, loc0=None, params=None, balloon=None, alt_change_fit=10, output_figs=False):

	time0 = time.time()

	############################################################################################################

	# set general parameters and initial conditions
	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, check_sigmas, params, balloon = pyb_io.set_params(params=params, balloon=balloon)
	lat0, lon0, alt0 = loc0
	lat_rad, lon_rad, all_alts = [np.radians(lat0)], [np.radians(lon0)], [alt0]

	current_time = utc_hour

	if lon_rad[0] < 0:
		lon_rad[0] += 2*np.pi
	if lat_rad[0] < 0:
		lat_rad[0] += 2*np.pi

	keys = list(data.keys())
	time_key0 = (int(int(keys[0][15:19])/100.) + int(keys[0][20:23]))

	alts = data[keys[0]]['altitudes'] # alts, lats and lons are the same for all weather files (if we dont give it different areas)
	data_lats, data_lons  = np.radians(data[keys[0]]['lats']), np.radians(data[keys[0]]['lons'])
	size = int(np.sqrt(len(data_lats)))

	if check_sigmas:
		data_lats_err, data_lons_err = np.radians(data[keys[0]]['lats_err']), np.radians(data[keys[0]]['lons_err'])

	speeds, z_speeds, omegas, temperatures, sigmas_u, sigmas_v, grid_spreads_u, grid_spreads_v, figs, f1s, f2s = [], [], [], [], [], [], [], [], [], [], []

	delta_ts, total_time, dists, dists_u, dists_v = [utc_hour - time_key0], [0], [0], [0], [0]

	index, i, timer = 0, 0, 0

	diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons - lon_rad[-1])**2)
	ini_grid_i, = np.where(diff == diff.min())
	ini_diff = diff[ini_grid_i[0]]

	loc_diffs = [ini_diff]

	z_speed_prop, omega_prop, temp_prop = 'z_speeds', 'omegas', 'temperatures' 
	sigma_u_prop, sigma_v_prop = 'u_wind_errs', 'v_wind_errs'

	stage = 1

	if not descent_only:
		sys.stdout.write('\r')
		sys.stdout.flush()
		sys.stdout.write('Calculating ascent...'.ljust(60) + '\r')
		sys.stdout.flush()
		time.sleep(0.2)

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

			############################################################################################################

			if stage == 1:
				x_prop, y_prop, t_prop, speed_prop,  = 'dxs_up', 'dys_up', 'ascent_time_steps', 'ascent_speeds'
			if stage == 2:
				x_prop, y_prop = 'u_winds', 'v_winds'
			if stage == 3:
				x_prop, y_prop, t_prop, speed_prop = 'dxs_down', 'dys_down', 'descent_time_steps', 'descent_speeds'

			############################################################################################################

			if interpolate:
				(t1, f1), (t2, f2) = calc_time_frac(current_time=current_time, weather_files=keys)
			else:
				t1, f1, t2, f2 = keys[index], 1, keys[index], 0

			# Find the closest grid point
			diff = np.sqrt((data_lats - lat_rad[-1])**2 + (data_lons - lon_rad[-1])**2)
			grid_i, = np.where(diff == diff.min())
			grid_i = grid_i[0]

			min_diff = diff[grid_i]

			grids = [grid_i-1, grid_i+1, grid_i-len(set(data_lats)), grid_i+len(set(data_lats))]

			if check_sigmas:
				diff = np.sqrt((data_lats_err - lat_rad[-1])**2 + (data_lons_err - lon_rad[-1])**2)
				grid_i_err, = np.where(diff == diff.min())
				grid_i_err = grid_i_err[0]

			# bilinear interpolation or not
			# dx = f1*data[t1][x_prop][grid_i, i] + f2*data[t2][x_prop][grid_i, i]
			# dy = f1*data[t1][y_prop][grid_i, i] + f2*data[t2][y_prop][grid_i, i]
			dx = f1*pyb_aux.bilinear_interpolation(grid_i, i, lon_rad[-1], lat_rad[-1], data_lons, data_lats, data[t1][x_prop], resolution) \
			 + f2*pyb_aux.bilinear_interpolation(grid_i, i, lon_rad[-1], lat_rad[-1], data_lons, data_lats, data[t2][x_prop], resolution)
			dy = f1*pyb_aux.bilinear_interpolation(grid_i, i, lon_rad[-1], lat_rad[-1], data_lons, data_lats, data[t1][y_prop], resolution) \
			 + f2*pyb_aux.bilinear_interpolation(grid_i, i, lon_rad[-1], lat_rad[-1], data_lons, data_lats, data[t2][y_prop], resolution)

			z_speed = f1*data[t1][z_speed_prop][grid_i, i] + f2*data[t2][z_speed_prop][grid_i, i]
			omega = f1*data[t1][omega_prop][grid_i, i] + f2*data[t2][omega_prop][grid_i, i]
			T = f1*data[t1][temp_prop][grid_i, i] + f2*data[t2][temp_prop][grid_i, i]

			grid_spread_u = f1*np.std([data[t1][x_prop][grid_i-1, i], data[t1][x_prop][grid_i+1, i], data[t1][x_prop][grid_i-len(set(data_lats)), i], data[t1][x_prop][grid_i+len(set(data_lats)), i]]) + \
			f2*np.std([data[t2][x_prop][grid_i-1, i], data[t2][x_prop][grid_i+1, i], data[t2][x_prop][grid_i-len(set(data_lats)), i], data[t2][x_prop][grid_i+len(set(data_lats)), i]])

			grid_spread_v = f1*np.std([data[t1][y_prop][grid_i-1, i], data[t1][y_prop][grid_i+1, i], data[t1][y_prop][grid_i-len(set(data_lats)), i], data[t1][y_prop][grid_i+len(set(data_lats)), i]]) + \
			f2*np.std([data[t2][y_prop][grid_i-1, i], data[t2][y_prop][grid_i+1, i], data[t2][y_prop][grid_i-len(set(data_lats)), i], data[t2][y_prop][grid_i+len(set(data_lats)), i]])

			delta_t = current_time - (int(int(t1[15:19])/100.) + int(t1[20:23]))

			if stage != 2:
				# bilinear interpolation or not
				# dt = f1*data[t1][t_prop][grid_i, i] + f2*data[t2][t_prop][grid_i, i]
				dt = f1*pyb_aux.bilinear_interpolation(grid_i, i, lon_rad[-1], lat_rad[-1], data_lons, data_lats, data[t1][t_prop], resolution) \
				 + f2*pyb_aux.bilinear_interpolation(grid_i, i, lon_rad[-1], lat_rad[-1], data_lons, data_lats, data[t2][t_prop], resolution)

				speed = f1*data[t1][speed_prop][grid_i, i] + f2*data[t2][speed_prop][grid_i, i]
			else:
				dt = 5
				speed = 0
				dx *= dt
				dy *= dt
				grid_spread_u *= dt
				grid_spread_v *= dt

			if check_sigmas:
				sigma_u = (f1*data[t1][sigma_u_prop][grid_i_err, i] + f2*data[t2][sigma_u_prop][grid_i_err, i])*dt
				sigma_v = (f1*data[t1][sigma_v_prop][grid_i_err, i] + f2*data[t2][sigma_v_prop][grid_i_err, i])*dt

		if stage == 1:
			if not descent_only:
				if all_alts[-1] >= data[keys[index]]['max_altitudes'][grid_i]:
					sys.stdout.write('\r')
					sys.stdout.flush()
					if drift_time == 0:
						stage += 2
						sys.stdout.write('Calculating descent...'.ljust(60) + '\r')
					else:
						stage += 1
						sys.stdout.write('Calculating drift trajectory...'.ljust(60) + '\r')
					sys.stdout.flush()
					time.sleep(0.2)
					final_i = 0
					continue
				i += 1
			else:
				# Mimick the movement during ascent if we only want descent
				if alts[i] >= alt0:
					if drift_time == 0:
						stage += 2
					else:
						stage += 1
					final_i = i
					continue
				i += 1
		elif stage == 2:
			timer += dt
			if timer == drift_time*60:
				stage += 1
				continue
		elif stage == 3:
			elevation = pyb_aux.get_elevation(lat=np.degrees(lat_rad[-1]), lon=np.degrees(lon_rad[-1]))
			if i == 0 or alts[i] <= elevation:
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
			z_speeds.append(z_speed)
			omegas.append(omega)
			temperatures.append(T)
			lat_rad.append(lat)
			lon_rad.append(lon)
			dists.append(dist)
			total_time.append(dt)
			all_alts.append(alt)
			grid_spreads_u.append(grid_spread_u)
			grid_spreads_v.append(grid_spread_v)

			if check_sigmas:
				sigmas_u.append(sigma_u)
				sigmas_v.append(sigma_v)

			if not interpolate:
				delta_ts.append(delta_t)

			loc_diffs.append(min_diff)

			f1s.append(f1)
			f2s.append(f2)

	# get more accurate end-point based on elevation data
	new_end_point, new_alt = pyb_aux.get_endpoint(data=(np.degrees(np.array(lat_rad)), np.degrees(np.array(lon_rad)), np.array(all_alts), np.array(dists)))

	lat_rad, lon_rad, all_alts, total_time, dists, loc_diffs = lat_rad[:-1], lon_rad[:-1], all_alts[:-1], total_time[:-1], dists[:-1], loc_diffs[:-1]

	diffs = np.abs(np.degrees(np.array(lat_rad)) - new_end_point[0])
	diff = np.sqrt((np.degrees(np.array(lat_rad)) - new_end_point[0])**2 + (np.degrees(np.array(lon_rad)) - new_end_point[1])**2)
	index = np.where(diff == min(diff))[0][0]

	lat_rad, lon_rad, all_alts, dists, speeds, z_speeds, omegas, temperatures, total_time, loc_diffs, grid_spreads_u, grid_spreads_v, f1s, f2s = \
		lat_rad[:index], lon_rad[:index], all_alts[:index], dists[:index+1], speeds[:index+1], z_speeds[:index+1], omegas[:index+1], temperatures[:index+1],\
		 total_time[:index+1], loc_diffs[:index+1], grid_spreads_u[:index+1], grid_spreads_v[:index+1], f1s[:index+1], f2s[:index+1]

	if check_sigmas:
		sigmas_u, sigmas_v = sigmas_u[:index+1], sigmas_v[:index+1]

	if not interpolate:
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
	output['speeds'] = np.array(speeds)
	output['z_speeds'] = np.array(z_speeds)
	output['omegas'] = np.array(omegas)
	output['temperatures'] = np.array(temperatures)
	output['times'] = np.cumsum(np.array(total_time))/60 # to minutes
	output['distance'] = np.sum(np.array(dists))
	output['grid_spreads_u'] = np.array(grid_spreads_u)
	output['grid_spreads_v'] = np.array(grid_spreads_v)

	if check_sigmas:
		output['sigmas_u'] = np.array(sigmas_u)
		output['sigmas_v'] = np.array(sigmas_v)

	if not interpolate:
		output['delta_ts'] = np.array(delta_ts)

	output['loc_diffs'] = np.array(loc_diffs)
	output['f1s'] = np.array(f1s)
	output['f2s'] = np.array(f2s)

	start_lat, start_lon = lat0, lon0
	end_lat, end_lon = float(output['lats'][-1]), float(output['lons'][-1])

	dx_end = haversine((start_lat, start_lon), (start_lat, end_lon))
	if end_lon < start_lon:
		dx_end *= -1
	dy_end = haversine((start_lat, start_lon), (end_lat, start_lon))
	if end_lat < start_lat:
		dy_end *= -1

	mean_theta = np.degrees(np.arctan2(dy_end, dx_end))

	output['mean_direction'] = mean_theta

	# print out relevant quantities
	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Trajectories calculated'.ljust(60))
	sys.stdout.write('\n')

	print('Maximum altitude: ' + str(np.max(all_alts)) + ' m')
	print('Landing location: (%.6f, %.6f)' % (output['lats'][-1], output['lons'][-1]))
	print('Flight time: %d min' % (int(output['times'][-1])) + ', distance travelled: %.1f' % output['distance'] + ' km')
	if check_sigmas:
		total_sigma = np.sqrt(np.sqrt(np.sum(output['sigmas_u']**2))**2 + np.sqrt(np.sum(output['sigmas_v']**2))**2)/1000.
		print('Sigma from ensemble forecasts: %.3f km' % total_sigma)
	print('')

	return output, figs, used_weather_files, time_diffs

#################################################################################################################

# method to prepare the data for the calculations of the trajectory
def prepare_data(weather_file=None, loc0=None, current_time=None, balloon=None, descent_only=False, drift_time=0, vz_correct=False, check_sigmas=False, output_figs=False):

	if check_sigmas:

		model_data1 = read_data(loc0=loc0, weather_file=weather_file, descent_only=descent_only)
		err_data = pyb_aux.calc_gefs_errs(weather_file=weather_file, loc0=loc0, current_time=current_time, descent_only=descent_only)
		model_data2, figs_dict = calc_properties(data=model_data1, weather_file=weather_file, loc0=loc0, balloon=balloon, descent_only=descent_only, output_figs=output_figs)
		model_data2 = pyb_aux.add_uv_errs(main_data=model_data2, err_data=err_data)
		model_data3 = calc_displacements(data=model_data2, balloon=balloon, descent_only=descent_only, vz_correct=vz_correct)

	else:

		model_data1 = read_data(loc0=loc0, weather_file=weather_file, descent_only=descent_only)
		model_data2, figs_dict = calc_properties(data=model_data1, weather_file=weather_file, loc0=loc0, balloon=balloon, descent_only=descent_only, output_figs=output_figs)
		model_data3 = calc_displacements(data=model_data2, balloon=balloon, descent_only=descent_only, vz_correct=vz_correct)

	return model_data3, figs_dict

#################################################################################################################

# method to pull everything together and return the trajectory
def run_traj(weather_files=None, datestr=None, utc_hour=None, loc0=None, params=None, balloon=None, output_figs=False):

	used_weather_files = {}
	used_weather_files[utc_hour] = weather_files

	time_diffs = {}
	time_diffs[utc_hour] = [(int(int(weather_file[15:19])/100.) + int(weather_file[20:23])) % 24 - utc_hour for weather_file in weather_files]

	descent_only, next_point, interpolate, drift_time, resolution, vz_correct, hr_diff, check_sigmas, params, balloon = pyb_io.set_params(params=params, balloon=balloon)

	data_array = {}
	figs1 = []
	for weather_file in weather_files:

		model_data, figs_dict = prepare_data(weather_file=weather_file, loc0=loc0, current_time=utc_hour, balloon=balloon, descent_only=descent_only, vz_correct=vz_correct, check_sigmas=check_sigmas, output_figs=output_figs)
		data_array[weather_file] = model_data
		figs1.append(figs_dict)

	trajectories, figs2, used_weather_files, time_diffs = calc_movements(data=data_array, used_weather_files=used_weather_files, time_diffs=time_diffs, datestr=datestr, utc_hour=utc_hour, loc0=loc0, params=params, balloon=balloon, output_figs=output_figs)

	for fig_dict in figs2:
		figs1.append(fig_dict)

	return trajectories, figs1, used_weather_files, time_diffs

#################################################################################################################
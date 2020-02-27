"""
Auxiliary functions used in pyBalloon
"""

from haversine import haversine
from osgeo import gdal
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.io import ascii
import numpy as np
import sys, os
import math
import time

import get_gfs
import pyb_io

import param_file as p

#################################################################################################################

def all_and(data):
	"""
	Method to create logical and for a list of arrays.
	
	Arguments
	=========
	data : list
		List of Numpy boolean arrays
	"""

	result = data.pop()
	for d in data:
		result = np.logical_and(result, d)

	return result

#################################################################################################################

def earth_radius(lat_rad):
	"""
	Calculate Earth WGS84 based radius on a given latitude.

	Arguments
	=========
	lat_rad : float
		Latitude in radians
	"""
	
	# WGS84 reference ellipsoid
	a = 6378.137	 # Earth equatorial radius, km
	b = 6356.7523142 # Earth polar radius, km
	
	cos_lat = np.cos(lat_rad)
	sin_lat = np.sin(lat_rad)

	r = np.sqrt(((a*a*cos_lat)**2 + (b*b*sin_lat)**2) / ((a*cos_lat)**2 + (b*sin_lat)**2))

	return r

#################################################################################################################

def air_density(data):
	"""
	Calculate air density from pressure and temperature.
	
	Arguments
	=========
	data : dict
		Dictonary returned by read_gfs_data()
	Return
		Air densities (as Numpy array) for each altitude step in data.
	"""

	ps = data['pressures'] # Pa
	Ts = data['temperatures'] # K

	if ps.shape != Ts.shape:
		x, y = ps.shape
		rho = [np.array((ps[:, i] * p.M_air)/(p.R0 * Ts)) for i in range(0, y)]
		rho = np.array(rho).transpose()
	else:
		rho = (ps * p.M_air)/(p.R0 * Ts)

	return rho # kg m-3

#################################################################################################################

def data_interpolation(data, alt0, step, descent_only=False, output_figs=False):
	"""
	Interpolate (and extrapolate in the low end, if mode='spline' is used) vertical data from alt0 to maximum level present in the data.

	Arguments
	=========
	data : dict
		Dictionary returned by read_gfs_data()
	alt0 : float
		Starting altitude of the interpolation in meters from WGS84 reference ellipsoid
	step : float
		Interpolation altitude step
	mode : string (optional)
		Interpolation method. Supported methods are 'spline' and 'linear'. Default: 'spline'
	Return:
		Interpolated version of the data.
	"""

	altitudes = data['altitudes']

	new_data = {}

	if descent_only:
		new_data['altitudes'] = np.arange(alt0 % step, alt0 + step, step)
		# new_data['altitudes'] = np.logspace(math.log(1, 10), math.log(alt0+1000, 10), 400) # logarithmic steps
	else:
		new_data['altitudes'] = np.arange(alt0, altitudes.max() + step, step)

	new_data['lats'] = data['lats']
	new_data['lons'] = data['lons']

	checks = ['altitudes', 'lats', 'lons']

	for key in data.keys():		
		if key not in checks:

			arr = []
			d = data[key]

			try:
				x, y = d.shape

			except ValueError:

				x = 0
				y, = d.shape

			if not descent_only:
				if x > 0:
					for i in range(0, y):
						ok_idxs = altitudes[:, i] >= alt0
						tck = interpolate.splrep(altitudes[ok_idxs, i], d[ok_idxs, i])
						arr.append(np.array(interpolate.splev(new_data['altitudes'], tck)))
				else:
					tck = interpolate.splrep(altitudes, d)
					arr.append(np.array(interpolate.splev(new_data['altitudes'], tck)))

			elif descent_only:					   
				if x > 0:
					for i in range(0, y):
						ok_idxs = altitudes[:, i] <= alt0
						tck = interpolate.splrep(altitudes[ok_idxs, i], d[ok_idxs, i])
						arr.append(np.array(interpolate.splev(new_data['altitudes']-0.5*step, tck)))
				else:
					tck = interpolate.splrep(altitudes, d)
					arr.append(np.array(interpolate.splev(new_data['altitudes']-0.5*step, tck)))

			new_data[key] = np.array(arr)

	figs = {}

	if output_figs:

		ind1 = np.where(altitudes[:, 0] == alt0)[0][0]
		ind2 = len(new_data['altitudes'])-1

		for key in data.keys():
			if key not in checks:

				fig = plt.figure()

				plt.axvline(alt0, linewidth=1, linestyle='--', label='Initial alt.')
				plt.plot(altitudes[:, 0], data[key][:, 0], 'ro--', label='Before interp.', markersize=3.5)
				plt.plot(altitudes[:, 0][ind1], data[key][:, 0][ind1], 'go', markersize=5, label='Inserted at alt. 0')
				plt.plot(new_data['altitudes'][:ind2+1], new_data[key][0][:ind2+1], 'bo', markersize=0.5, label='After interp.')
				plt.ylabel(key.capitalize().replace('_', ' '), fontsize=15)
				plt.xlabel('Altitude [m]', fontsize=15)
				plt.legend(loc='best')
				plt.grid(True)
				plt.tight_layout()

				figs[key] = fig

				plt.close()

	return new_data, figs

#################################################################################################################

def lift(data, mass):
	"""
	Calculate effective lift (force, in Newtons) caused by the balloon.

	Arguments
	=========
	data : dict
		Dictionary containing 'altitudes' (meters), balloon 'volumes' (m^3) and 'air_densities' (kg/m^3)
	mass: float
		Mass of the whole balloon
	Return:
		Resultant lift force.
	"""

	h = data['altitudes']
	V_b = data['balloon_volumes']
	rho_air = data['air_densities']
	g = p.g_0 * (p.R_e / (p.R_e + h))**2 # Gravitational acceleration at height h

	F_lift = g * (V_b*rho_air - mass)

	return F_lift

#################################################################################################################

def balloon_volume(data):
	"""
	Calculate volume of a sphere.

	Arguments
	=========
	data : dict
		Dictionary containing 'balloon_radii' (in meters)
	Return:
		Volume of the balloon at each level in data
	"""

	r = data['balloon_radii']
	V = 4/3. * np.pi * r**3

	return V

#################################################################################################################

def balloon_volume_ideal_gas(data, gas_mass, gas_molar_mass=p.M_helium):
	"""
	Calculate gas (~balloon) volume based on ideal gas law: pV = nRT.

	Arguments
	=========
	data : dict
		Dictionary returned by read_gfs_data()
	gas_mass : float
		Mass of the gas in question
	gas_molar_mass : float (optional)
		Gas molar mask (g/mol). Default: 4.002602 (Helium)
	Return:
		Gas volume at each level in data
	"""

	m = gas_mass # 
	M = gas_molar_mass/1000. # default to Helium, convert to kg/mol

	Ts = data['temperatures']
	ps = data['pressures'] # pressure in Pascals

	# gas volume without the balloon
	V = m*p.R0*Ts/(M*ps) # pV = nRT, n = m/M

	return V

#################################################################################################################

def burst_altitude(data, burst_radius):
	"""
	Find the altitude where balloon radius gets greater than the given burst radius.

	Arguments
	=========
	data : dict
		Dictionary containing 'altitudes' and corresponding 'balloon_radii'
	burst_radius : float
		Balloon burst radius
	Return:
		Altitude of burst and corresponding array index
	"""

	radii = data['balloon_radii']

	alt_out = []
	i_out = []
	i = 0
	for radius in radii:
		diff = np.abs(radius-burst_radius)
		idx = np.where(diff == diff.min())	
		i_out.append(idx[0][0])
		alt_out.append(data['altitudes'][idx])
		i+=1

	return np.array(alt_out), np.array(i_out)

#################################################################################################################

def neutral_buoyancy_level(data):
	"""
	Find the level of neutral buoyancy (or more precise, the level where effective lift is closest to zero).

	Arguments
	=========
	data : dict
		Dictionary containing 'altitudes' and corresponding 'lifts'
	Return:
		Altitude of neutral buoyancy and corresponding array index
	"""

	alt_out = []
	i_out = []
	idx = 0
	for lft in data['lifts']:
		lft = np.abs(lft)
		idx2 = np.where(lft == lft.min())
		i_out.append(idx2[0][0])
		alt_out.append(data['altitudes'][idx, idx2])
		idx += 1

	return np.array(alt_out), np.array(i_out)

#################################################################################################################

def ascent_speed(data, mass, Cd=p.Cd_sphere):
	"""
	Calculate the rate of ascent (in m/s) for the inflated balloon at given levels.

	Arguments
	=========
	data : dict
		Dictionary with 'altitudes' and corresponding 'air_densities', 'balloon_radii' and 'balloon_volumes'
	mass : float
		Full balloon mass (in kg)
	Cd : float (optional)
		Coefficient of drag. Default: 0.47 (sphere)
	Return:
		Ascent speed (m/s) at every level in input data.
	"""

	m = mass
	rho = data['air_densities']
	A = np.pi*data['balloon_radii']**2
	V = data['balloon_volumes']
	h = data['altitudes']

	g = p.g_0 * (p.R_e / (p.R_e + h))**2

	Fb = V*rho*g # byoyance
	Fg = m*g	 # gravity
	F = Fb-Fg

	# Set negative buyoyancies to zero (we won't get there)
	idxs = np.where(F <= 0)
	if len(idxs[0]) > 0:
		idx = idxs[0][0]
		F[idx:] = 1e-30

	v = np.sqrt(2*F/(rho*Cd*A)) # m/s

	return v

#################################################################################################################

def descent_speed(data, mass, Cd, areas, alt_step, change_alt=None):
	"""
	Calculate the rate of descent for deflated (burst) balloon with 1 or 2 different sized parachutes with given areas, change altitude and drag-coefficient.

	Arguments
	=========
	data : dict
		Dictionary with 'altitudes', and corresponding 'air_densities'
	mass: float
		Mass of the payload + assumed remnants of the balloon
	Cd : floats in tuple
		Coefficients of drag for one or two parachutes in a tuple
	areas: floats in tuple
		Effective areas (in m^2) of one or two parachutes in a tuple
	change_alt : float
		Altitude where first parachute is changed to the second one. If None, only one parachute is used. Default: None
	Return:
		Rate of descent (m/s) for each level in input data
	"""

	m = mass
	h = data['altitudes']
	rho = data['air_densities']
	g = p.g_0 * (p.R_e / (p.R_e + h))**2 # Gravitational acceleration at height h

	# speeds = []
	# if change_alt is not None:
	# 	idxs = h < change_alt
	# 	for rho in data['air_densities']:
	# 		v = np.sqrt(2*m*g/(rho*Cd*areas[0])) # m/s
	# 		v[idxs] = np.sqrt(2*m*g[idxs]/(rho[idxs]*Cd*areas[1]))
	# 		speeds.append(v)
	# else:
	# 	for rho in data['air_densities']:
	# 		v = np.sqrt(2*m*g/(rho*Cd*areas)) # m/s
	# 		speeds.append(v)

	if change_alt is not None:
		idxs = h < change_alt
		v = np.sqrt(2*m*g/(rho*Cd*areas[0])) # m/s
		v[idxs] = np.sqrt(2*m*g[idxs]/(rho[idxs]*Cd*areas[1]))
	else:
		v = np.sqrt(2*m*g/(rho*Cd*areas)) # m/s

	speeds = v

	return -1*np.array(speeds)

#################################################################################################################	

def mooney_rivlin(data, radius_empty, radius_filled, thickness_empty, gas_molar_mass=p.M_helium):
	"""
	Calculate balloon radii for given pressures/temperatures/initial conditions using inversion of Mooney-Rivlin equation.
	See description of the equation at: http://www.zmatt.net/weather-balloon-physics/

	Arguments
	=========
	data : dict
		Dictionary containing 'pressures' and 'temperatures'
	radius_empty : float
		Radius of the empty balloon
	radius_filled : float
		Radius when filled at ground level
	thickness_empty : float
		Balloon rubber initial thickness
	gas_molar_mass : float (optional)
		Molar mass (g/mol) of the gas used to fill the balloon. Default: 4.002602 (Helium)
	Return:
		Radius of the balloon at each level in input data, and the mass of the gas.
	"""

	r0 = radius_empty # in meters
	r1 = radius_filled # in meters
	t0 = thickness_empty # in meters
	M = gas_molar_mass / 1000. # convert to kg/mol
	ps = data['pressures']
	p0 = ps[0]
	Ts = data['temperatures']
	T0 = Ts[0]

	mu = 300000. # balloon shear modulus in Pascals
	alfa = 10./11.

	# Amount of gas in moles
	n = (4/3. * np.pi * r1**3)/(p.R0*T0) * (p0 - 2*mu*(t0/r0)*((r0/r1) - (r0/r1)**7) * (1 + (1/alfa - 1) * (r1/r0)**2))
	gas_mass = n*M

	# Solve balloon-radius-roots for each height level
	
	# Constants for the polynomials
	a8 = (1/alfa-1)*2*mu*t0/(r0**2)
	a6 = 2*mu*t0
	a2 = -2*mu*t0*(1/alfa-1)*(r0**4)
	a0 = -2*mu*t0*(r0**6)
	
	all_radii = []

	try:
		x, y = ps.shape
		
		for i in range(0, x):
			radii = []
			for j in range(0, y):
				a4 = -3*n[j]*p.R0/(4*np.pi)
				# 8th degree polynomial
				poly = [a8,		# r^8
						ps[i,j],	# r^7
						a6,		# r^6
						0,		 # r^5
						a4*Ts[i,j], # r^4
						0,		 # r^3
						a2,		# r^2
						0,		 # r^1
						a0]		# r^0

				roots = np.roots(poly)
		
				for r in roots:
					if r.real > 0 and r.imag == 0:
						radii.append(r.real)
			all_radii.append(np.array(radii))


	except ValueError:		
		for i in range(0, len(p)):
			a4 = -3*n*p.R0/(4*np.pi)
			# 8th degree polynomial
			poly = [a8,		# r^8
					ps[i],	# r^7
					a6,		# r^6
					0,		 # r^5
					a4*Ts[i], # r^4
					0,		 # r^3
					a2,		# r^2
					0,		 # r^1
					a0]		# r^0

			roots = np.roots(poly)
		
			for r in roots:
				if r.real > 0 and r.imag == 0:
					all_radii.append(r.real)

	all_radii = np.array(all_radii)

	return all_radii, gas_mass

#################################################################################################################

# find the elevation at the location (the loc closest in the grid)
def get_elevation(lon, lat):

	srtm_dir = p.path + p.elevation_data_folder

	if lat > 180:
		lat -= 360
	if lon > 180:
		lon -= 360

	if np.abs(lat) < 60: # better accuracy, but slower!

		lon_box = int(math.ceil((lon + 180) / 5.))
		lat_box = 25 - int(math.ceil((lat + 60) / 5.))

		file = srtm_dir + 'tif_files/srtm_' + str(lon_box).zfill(2) + '_' + str(lat_box).zfill(2) + '.tif'

		ds = gdal.Open(file, gdal.GA_ReadOnly)
		elevations = ds.GetRasterBand(1).ReadAsArray()

		nrows, ncols = elevations.shape
		x0, dx, dxdy, y0, dydx, dy = ds.GetGeoTransform()

		lons = np.arange(x0, x0 + dx*ncols, dx)
		lats = np.arange(y0, y0 + dy*nrows, dy)                                     

		diff1 = np.abs(lons - lon)
		diff2 = np.abs(lats - lat)

		i1 = np.where((diff1 == np.min(diff1)))[0][0]
		i0 = np.where((diff2 == np.min(diff2)))[0][0]

		elevation = elevations[i0][i1]

	else:

		SAMPLES = 1201  # Change this to 3601 for SRTM1
		HGTDIR = p.path + p.elevation_data_folder + '/hgt_files/'

		if lat >= 0:
			ns = 'N'
		elif lat < 0:
			ns = 'S'

		if lon >= 0:
			ew = 'E'
		elif lon < 0:
			ew = 'W'

		hgt_file = "%(ns)s%(lat)02d%(ew)s%(lon)03d.hgt" % {'lat': abs(lat), 'lon': abs(lon), 'ns': ns, 'ew': ew}
		hgt_file_path = os.path.join(HGTDIR, hgt_file)

		if not os.path.isfile(hgt_file_path):
			print(str(hgt_file) + ' does not exist!')
			elevation = -32768

		else:
			with open(hgt_file_path, 'rb') as hgt_data:
				# HGT is 16bit signed integer(i2) - big endian(>)
				elevations = np.fromfile(hgt_data, np.dtype('>i2'), SAMPLES*SAMPLES).reshape((SAMPLES, SAMPLES))

				lat_row = int(round((lat - int(lat)) * (SAMPLES - 1), 0))
				lon_row = int(round((lon - int(lon)) * (SAMPLES - 1), 0))

				elevation = elevations[SAMPLES - 1 - lat_row, lon_row].astype(int)

	return elevation

#################################################################################################################

def get_endpoint(data=None, run=None, filename=None, params=None):
	""" 
	Calculate more accurate endpoint using data from: http://viewfinderpanoramas.org/
	PyBalloon does not know the elevation of given locations, so trajectories can continue underground untill i = 0 (see pyb_traj)

	Arguments
	=========
	data: dict
		Dictionary containing 'latitudes', 'longitudes', 'altitudes', and 'distances'. If None, read data from run and filename
	run : string
		String indicating which run folder the trajectory data is stored in, if None use last folder created
	filename : string
		Name of the trajectory file which contains the 'latitudes', 'longitudes', 'altitudes', and 'distances'.
	params : list
		List of parameters used for given trajectory. If None, use parameters from param_file.py
	Return:
	   Latitude and longitude of new endpoint and elevation at this longitude and latitude
	"""

	if params == None:
		descent_only = p.descent_only
		drift_time = p.drift_time
	else:
		descent_only = bool(params[0])
		drift_time = float(params[1])

	if data == None and run != None:

		dir_pred = p.path + p.output_folder + p.traj_folder + run + '/'
		data  = ascii.read(dir_pred + filename)

		lats = data['lats']
		lons = data['lons']
		alts = data['alts']
		dists = data['dists']

	else:

		lats, lons, alts, dists = data

	elevations = []

	for i in range(1, len(lats)):
		elevation = get_elevation(lat=lats[-i], lon=lons[-i])
		elevations.append(elevation)
		if elevation < alts[-i]:
			break

	if i == 1:
		return (lats[-i], lons[-i]), elevation

	inds = (-i, -i + 1)

	if list(set(dists)) == [0]:
		newlon, newlat = lons[inds[0]], lats[inds[0]]

	else:

		dlon = lons[inds[1]] - lons[inds[0]]
		dlat = lats[inds[1]] - lats[inds[0]]

		x1 = dists[inds[1]]

		y1 = alts[inds[0]] - alts[inds[1]]
		y2 = elevations[-2] - alts[inds[1]]

		dx = x1*(y2/y1)
		f = dx/x1

		newlon = lons[inds[1]] - f*dlon
		newlat = lats[inds[1]] - f*dlat

	return (newlat, newlon), get_elevation(lat=newlat, lon=newlon)

#################################################################################################################

def calc_gefs_errs(weather_file=None, current_time=None, loc0=None, descent_only=False):
	"""
	Calculate std and mean from the GEFS ensemble forecasts

	Arguments
	=========
	weather_file : string
		Name of main weather forecast model (GFS)
	current_time : float
		Current time of the trajectory
	loc0 : floats in tuple
		(latitude in degrees, longitude in degrees, altitude in km) of initial point
	descent_only : bool
		Option to start the trajectory at its highest point. If True, the trajectory only has a descent phase
	Return:
		Dictionary containing std and mean data from the GEFS ensemble forecasts
	"""

	import param_file as p

	indir = p.path + p.weather_data_folder + p.GEFS_folder

	tile_size = p.tile_size
	lat0, lon0, alt0 = loc0
	area = (lat0 + tile_size / 2, lon0 - tile_size / 2, lat0 - tile_size / 2, lon0 + tile_size / 2)
	tlat, llon, blat, rlon = area

	datestr = weather_file[6:14]
	hhhh, hhh = int(int(weather_file[15:19]) / 100), int(weather_file[20:23])
	hour = hhhh + hhh

	hrs = get_gfs.get_closest_hr(utc_hour=hour)
	closest_hhhh, closest_hhh1, closest_hhh2 = hrs[0], hrs[1], hrs[2]

	files = [filename for filename in os.listdir(indir)]

	count = 0
	while count < 2:

		files = [filename for filename in os.listdir(indir)]

		if 'gespr_4_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2' in files and 'geavg_4_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2' in files:

			sys.stdout.write('\r')
			sys.stdout.flush()
			sys.stdout.write('Reading GEFS file...'.ljust(60) + '\r')
			sys.stdout.flush()
			time.sleep(0.2)

			data1 = pyb_io.read_gefs_file(fname='gespr_4_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2', area=area, alt0=alt0, descent_only=descent_only)
			data2 = pyb_io.read_gefs_file(fname='geavg_4_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2', area=area, alt0=alt0, descent_only=descent_only)
			
			data = {}
			data['lats'] = data1['lats']
			data['lons'] = data1['lons']
			data['u_winds'] = data1['u_winds']
			data['v_winds'] = data1['v_winds']
			data['altitudes'] = data2['altitudes']
			data['pressures'] = data1['pressures']

			data = data_interpolation(data=data, alt0=alt0, step=100, descent_only=descent_only, output_figs=False)[0]

			return data

		if 'gens-a_3_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(closest_hhh2).zfill(3) + '_00.grb2' in files:

			u_winds, v_winds, altitudes = {}, {}, {}

			sys.stdout.write('\r')
			sys.stdout.flush()
			sys.stdout.write('Reading 20 GEFS files...'.ljust(60) + '\r')
			sys.stdout.flush()
			time.sleep(0.2)

			for no in range(0, 20):

				ens_fname = 'gens-a_3_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(closest_hhh2).zfill(3) + '_' + str(no + 1).zfill(2) + '.grb2'
				data = pyb_io.read_gefs_file(fname=ens_fname, area=area, alt0=alt0, descent_only=descent_only)
				u_winds[no], v_winds[no], altitudes[no] = data['u_winds'], data['v_winds'], data['altitudes']

			u_winds_std, v_winds_std, altitudes_mean = {}, {}, {}
			lats, lons, pressures = data['lats'], data['lons'], data['pressures']

			for p in range(len(pressures)):

				key = pressures[p][0]
				u_winds_std[key] = []
				v_winds_std[key] = []
				altitudes_mean[key] = []

				for j in range(len(u_winds[0][p])):

					u_winds_std[key].append(np.std([u_winds[i][p][j] for i in range(0, 20)]))
					v_winds_std[key].append(np.std([v_winds[i][p][j] for i in range(0, 20)]))
					altitudes_mean[key].append(np.mean([altitudes[i][p][j] for i in range(0, 20)]))

			u_winds, v_winds, altitudes = [], [], []

			for key in pressures[:, 0]:

				uwnd, vwnd, alt = [], [], []
				uwnd.append(u_winds_std[key])
				vwnd.append(v_winds_std[key])
				alt.append(altitudes_mean[key])
				u_winds.append(np.hstack(uwnd))
				v_winds.append(np.hstack(vwnd))
				altitudes.append(np.hstack(alt))

			data = {}
			data['lats'] = lats
			data['lons'] = lons
			data['u_winds'] = np.array(u_winds)
			data['v_winds'] = np.array(v_winds)
			data['altitudes'] = np.array(altitudes)
			data['pressures'] = pressures

			data = data_interpolation(data=data, alt0=alt0, step=100, descent_only=descent_only, output_figs=False)[0]

			return data

		get_gfs.get_gefs_files(datestr=datestr, utc_hour=current_time)
		
		count += 1

#################################################################################################################

def add_uv_errs(main_data=None, err_data=None):
	"""
	Method to add error data from ensemble weather forecasts to main weather data

	Arguments
	=========
	main_data : dict
		Dictionary containing the data from the main weather forecast model (GFS)
	err_data : dict
		Dictionary containing the standard deviation and mean data from the twenty ensemble weather forecast models (GEFS)
	Return
		Dictionary made out of both the main and error data
	"""

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Adding u/v wind errors...'.ljust(60) + '\r')
	sys.stdout.flush()
	time.sleep(0.2)

	main_data['u_wind_errs'] = err_data['u_winds']
	main_data['v_wind_errs'] = err_data['u_winds']
	main_data['lats_err'] = err_data['lats']
	main_data['lons_err'] = err_data['lons']

	return main_data

#################################################################################################################

def calc_mean_travel_direction(lon0, lat0, end_lon, end_lat):
	"""
	Calculate mean direction of travel from start lat/lon to end lat/lon

	Arguments
	=========
	lon0 : float
		Initial longitude
	lat0 : float
		Initial latitude
	end_lon : float
		Final (predicted) longitude
	end_lat : float
		Final (predicted) latitude
	"""

	if end_lon > 180:
		end_lon -= 360

	dx_end = haversine((lat0, lon0), (lat0, end_lon))
	dy_end = haversine((lat0, lon0), (end_lat, lon0))

	theta = np.degrees(np.arctan2(dy_end, dx_end))

	if end_lat < lat0 and end_lon > lon0:
		theta = 360. - theta
	elif end_lat > lat0 and end_lon < lon0:
		theta = 180. - theta
	elif end_lat < lat0 and end_lon < lon0:
		theta += 180.

	return theta

#################################################################################################################

def find_bilinear_points(grid_i, i, lon_rad, lat_rad, data_lons, data_lats, resolution):
	"""
	Calculate locations and values of points to be used for bilinear interpolation given a list of data and point

	Arguments
	=========
	grid_i : float
		Index corresponding to closest longitude and latitude to given point
	i : float
		Index corresponding to current altitude
	lon_rad : float
		Longitude in radians of point at which we wish to know the value of the property
	lat_rad : float
		Latitude in radians of point at which we wish to know the value of the property
	data_lons : list
		List of longitudes in weather forecast model (GFS)
	data_lats : list
		List of latitudes in weather forecast model (GFS)
	resolution : float
		Resolution of weather forecast models (GFS) used
	Return
		Four tuples containing the x and y coordinates and values of the points to be used for bilinear interpolation
	"""

	curr_lat, curr_lon = data_lats[grid_i], data_lons[grid_i]

	lat_diff = curr_lat - lat_rad
	lon_diff = curr_lon - lon_rad

	lon_add = np.radians(-resolution*np.sign(lon_diff))
	lat_add = np.radians(-resolution*np.sign(lat_diff))

	if lon_diff != 0 and lat_diff != 0:
		lons = [curr_lon, curr_lon, curr_lon + lon_add, curr_lon + lon_add]
		lats = [curr_lat, curr_lat + lat_add, curr_lat, curr_lat + lat_add]
	else:
		if lon_diff == 0 and lat_diff != 0:
			lon_add += np.radians(resolution)
			lons = [curr_lon - lon_add, curr_lon - lon_add, curr_lon + lon_add, curr_lon + lon_add]
			lats = [curr_lat, curr_lat + lat_add, curr_lat, curr_lat + lat_add]
		elif lon_diff != 0 and lat_diff == 0:
			lat_add += np.radians(resolution)
			lons = [curr_lon, curr_lon + lon_add, curr_lon, curr_lon + lon_add]
			lats = [curr_lat - lat_add, curr_lat - lat_add, curr_lat + lat_add, curr_lat + lat_add]
		else:
			lons = [curr_lon, curr_lon, curr_lon, curr_lon]
			lats = [curr_lat, curr_lat, curr_lat, curr_lat]

	bottom_lat, bottom_lon, top_lat, top_lon = min(lats), min(lons), max(lats), max(lons)

	up_left = np.where((np.isclose(data_lats, top_lat)) & (np.isclose(data_lons, bottom_lon)))[0][0]
	low_left = np.where((np.isclose(data_lats, bottom_lat)) & (np.isclose(data_lons, bottom_lon)))[0][0]
	up_right = np.where((np.isclose(data_lats, top_lat)) & (np.isclose(data_lons, top_lon)))[0][0]
	low_right = np.where((np.isclose(data_lats, bottom_lat)) & (np.isclose(data_lons, top_lon)))[0][0]

	x1, y1, x2, y2 = data_lons[low_left], data_lats[low_left], data_lons[up_right], data_lats[up_right]

	return x1, y1, x2, y2, low_left, up_left, low_right, up_right

#################################################################################################################

def bilinear_interpolation(grid_i=None, i=None, lon_rad=None, lat_rad=None, data_lons=None, data_lats=None, prop=None, resolution=None, coords=None, inds=None):
	"""
	Calculate property at given latitude/longitude using bilinear interpolation

	Arguments
	=========
	grid_i : float
		Index corresponding to closest longitude and latitude to given point
	i : float
		Index corresponding to current altitude
	lon_rad : float
		Longitude in radians of point at which we wish to know the value of the property
	lat_rad : float
		Latitude in radians of point at which we wish to know the value of the property
	data_lons : list
		List of longitudes in weather forecast model (GFS)
	data_lats : list
		List of latitudes in weather forecast model (GFS)
	prop : list
		List containing the data we wish to interpolate
	resolution : float
		Resolution of weather forecast models (GFS) used
	coords : list
		List of lon and lat coordinates of points used for the bilinear interpolation
	inds : list
		List of indices of the points in the weather forecast grid
	Return
		Value of data at given latitude and longitude based on bilinear interpolation
	"""

	x = lon_rad
	y = lat_rad

	if coords == None or inds == None:
		x1, y1, x2, y2, low_left, up_left, low_right, up_right = find_bilinear_points(grid_i, i, lon_rad, lat_rad, data_lons, data_lats, resolution)
	else:
		x1, y1, x2, y2 = coords
		low_left, up_left, low_right, up_right = inds

	q11, q12, q21, q22 = prop[low_left, i], prop[up_left, i], prop[low_right, i], prop[up_right, i]
	points = (x1, y1, q11), (x1, y2, q12), (x2, y1, q21), (x2, y2, q22)

	points = sorted(points)
	(x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

	if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
		raise ValueError('points do not form a rectangle')
	if not x1 <= x <= x2 or not y1 <= y <= y2:
		raise ValueError('(x, y) not within the rectangle')

	return (q11 * (x2 - x) * (y2 - y) +
			q21 * (x - x1) * (y2 - y) +
			q12 * (x2 - x) * (y - y1) +
			q22 * (x - x1) * (y - y1)
		   ) / ((x2 - x1) * (y2 - y1) + 0.0)

#################################################################################################################

if __name__ == '__main__':

	time0 = time.time()

	lat, lon = 67.95563298887704, 311.74684661346487

	el = get_elevation(lon, lat)

	print(time.time() - time0)


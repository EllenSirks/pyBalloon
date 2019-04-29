from astropy.io import ascii
import numpy as np
import warnings
import pyb_traj
import get_gfs
import pyb_io
import time
import sys
import os

warnings.filterwarnings("ignore")

def retrospective(weather_file=None, datestr=None, utc_hour=None, lat0=None, lon0=None, alt0=None, descent_only=True, next_point='0'):

	time0 = time.time()

	in_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

	lat0 = float(lat0)
	lon0 = float(lon0)
	alt0 = float(alt0)

	if descent_only == 'True':
		descent_only = True
	else:
		descent_only = False

	if next_point == '2':
		ext = 'next2_points/'
	elif next_point == '1':
		ext = 'next_points/'		
	else:
		ext = 'all_points/'

	if weather_file is None:
		hhhh, hhh = get_gfs.get_closest_gfs_file(datestr=datestr, utc_hour=utc_hour, verbose=True)
		file = 'gfs_4_' + datestr +'_' + str(hhhh) + '_' +str(hhh)
	else:
		hhhh, hhh = get_gfs.get_gfs_file(weather_file=weather_file, verbose=True)
		file, datestr = weather_file[:-5], weather_file[6:14]

	# example balloon:
	balloon = {}
	balloon['altitude_step'] = 100.0 # meters
	balloon['equip_mass'] = 1.608 # kg
	balloon['balloon_mass'] = 1.50 # kg
	balloon['fill_radius'] = 2.122/2 # meters
	balloon['radius_empty'] = 2.41/2 # meters (flaccid body length - neck length)
	balloon['burst_radius'] = 9.44/2 # meters
	balloon['thickness_empty'] = 0.2 * 10**-3 # mm -> meters
	balloon['Cd_balloon'] = 0.5
	balloon['Cd_parachute'] = 0.5
	radius = 1.0
	balloon['parachute_areas'] = np.pi * np.array([radius])**2 # m^2
	# balloon['parachute_areas'] = np.array([radius*np.sqrt(2)])**2 # m^2
	balloon['parachute_change_altitude'] = None # meters
	balloon['simple_ascent_rate'] = 5.0 # m/s

	tile_size = 2. # degrees (read a tile this wide/high from the GFS grb2 file)
	area = (lat0 + (tile_size/2.), lon0 - (tile_size/2.), lat0 - (tile_size/2.), lon0 + (tile_size/2.)) # note top, left, bottom, right ordering for area

	model_data = pyb_io.read_gfs_single(in_dir + file, area=area, alt0=alt0, descent_only=descent_only, step = balloon['altitude_step'])
	print('GFS data read, %.1f s elapsed' % (time.time() - time0))

	loc0 = (lat0, lon0, alt0)

	trajectories = []
	for data in model_data:
		trajectories.append(pyb_traj.calc_movements(data, loc0, balloon, descent_only=descent_only, name=file[6:] + '_'+str(loc0), ext=ext))
	print('Trajectories calculated, %.1f s elapsed' % (time.time() - time0))

	traj_file = '/home/ellen/Desktop/SuperBIT/Weather_data/Trajectories/' + ext + file[6:] + '_'+str(loc0) +'_trajectory.dat'
	ascii.write([trajectories[0]['lats'], trajectories[0]['lons'], trajectories[0]['alts'], trajectories[0]['dists'], trajectories[0]['times'], trajectories[0]['descent_speeds']], traj_file, names=['lats', 'lons', 'alts', 'dists', 'times', 'descent_speeds'], overwrite=True)

	# highest point in main-run trajectory
	idx, = np.where(trajectories[0]['alts'] == np.max(trajectories[0]['alts']))
	latx = trajectories[0]['lats'][idx]
	lonx = trajectories[0]['lons'][idx]
	altx = trajectories[0]['alts'][idx]
	timex = trajectories[0]['times'][idx]

	print(latx, lonx, altx, '%.0f minutes' % (timex))

	other_info = [(latx, lonx, altx, 'Burst point', '%.0f minutes, %.0f meters' % (timex, altx))]

	kml_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/kml_files/' + ext
	kml_fname = kml_dir + file[6:] + '_' + str(loc0) + '.kml'
	pyb_io.save_kml(kml_fname, trajectories, other_info=other_info)

	f = open('/home/ellen/Desktop/SuperBIT/Weather_data/Endpoints/' + ext + file[6:] + '_' + str(loc0) + '_endpoint.dat','w+')
	f.write('lat lon alt\n')
	f.write(str(lat0) + ' ' + str(lon0) + ' ' + str(alt0) + '\n')
	f.write(str(trajectories[0]['lats'][-1]) + ' ' + str(trajectories[0]['lons'][-1]) + ' ' + str(trajectories[0]['alts'][-1]))
	f.close()

	print('Program finished in %.1f s' % (time.time() - time0))

if __name__ == "__main__":

	if sys.argv[1].endswith('grb2'):
		retrospective(weather_file=sys.argv[1], lat0=sys.argv[2], lon0=sys.argv[3], alt0=sys.argv[4], descent_only=sys.argv[5], next_point=sys.argv[6])
	else:
		retrospective(datestr=sys.argv[1], utc_hour=sys.argv[2], lat0=sys.argv[3], lon0=sys.argv[4], alt0=sys.argv[5], descent_only=sys.argv[6], next_point=sys.argv[7])

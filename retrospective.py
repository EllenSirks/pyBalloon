"""Retrospective pyBalloon"""

# Download GFS 0.5 degree archive data in grb2 format from the NOAA NCEI:
# https://nomads.ncdc.noaa.gov/data/gfs4/

import numpy as np
import sys
import time
import pyb_io
import pyb_traj
import get_gfs
import os

if __name__ == "__main__":

	# GFS doesn't seem to like negative longitudes
	# Use 0-360 instead

	datestr = sys.argv[1]
	utc_hour = sys.argv[2]

	lat0 = float(sys.argv[3]) #69.385818 # degrees
	lon0 = float(sys.argv[4]) #309.093365 = -50.906635 # degrees
	alt0 = float(sys.argv[5]) #0. # meters

	if sys.argv[6] == 'True':
		descent_only = True
	elif sys.argv[6] == 'False':
		descent_only = False

	in_dir = '/cosma5/data/durham/dc-sirk1/SuperBIT/Weather_data/'

	#hhhh, hhh = get_gfs.get_gfs_data(datestr, utc_hour, verbose=True)
        hhhh, hhh = '1200', '000'

	file = 'gfs_4_' + datestr +'_' + str(hhhh) + '_' +str(hhh)

	# if not os.path.isfile(indir + file):
		
	tile_size = 2. # degrees (read a tile this wide/high from the GFS grb2 file)

	time0 = time.time()	

	# example balloon:
	balloon = {}
	balloon['altitude_step'] = 50.0 # meters
	balloon['equip_mass'] = 1.608 # kg
	balloon['balloon_mass'] = 1.50 # kg
	balloon['fill_radius'] = 2.122/2 # meters
	balloon['radius_empty'] = 2.41/2 # meters (flaccid body length - neck length)
	balloon['burst_radius'] = 9.44/2 # meters
	balloon['thickness_empty'] = 0.2 * 10**-3 # mm -> meters
	balloon['Cd_balloon'] = 0.5
	balloon['Cd_parachute'] = 0.5
	#balloon['parachute_areas'] = np.pi * np.array([1.07])**2 # m^2
	balloon['parachute_areas'] = np.pi * np.array([1.0])**2 # m^2
	balloon['parachute_change_altitude'] = None # meters
	balloon['simple_ascent_rate'] = 5.0 # m/s

	# note top, left, bottom, right ordering for area
	model_data = pyb_io.read_gfs_single(in_dir + file, (lat0+(tile_size/2.), lon0-(tile_size/2.), lat0-(tile_size/2.), lon0+(tile_size/2.)), alt0, descent_only=descent_only)

	print('GFS data read, %.1f s elapsed' % (time.time() - time0))

	loc0 = (lat0, lon0, alt0)

	trajectories = []

	for data in model_data:
		trajectories.append(pyb_traj.calc_movements(data, loc0, balloon, descent_only=descent_only))

	print('Trajectories calculated, %.1f s elapsed' % (time.time() - time0))
        
	# highest point in main-run trajectory
	idx, = np.where(trajectories[0]['alts'] == np.max(trajectories[0]['alts']))
	latx = trajectories[0]['lats'][idx]
	lonx = trajectories[0]['lons'][idx]
	altx = trajectories[0]['alts'][idx]
	timex = trajectories[0]['times'][idx]
	print(latx, lonx, altx, '%.0f minutes' % (timex))
	other_info = [(latx, lonx, altx, 'Burst point', '%.0f minutes, %.0f meters' % (timex, altx))]

	klm_dir = '/cosma5/data/durham/dc-sirk1/SuperBIT/Weather_data'

	kml_fname = in_dir + file + '.kml'
	pyb_io.save_kml(kml_fname, trajectories, other_info=other_info)
	
	print('Program finished in %.1f s' % (time.time() - time0))

"""File with relevant running parameters"""

import numpy as np

# general parameters
descent_only = True
next_point = '1'
interpolate = False
drift_time = 0.
resolution = 0.5
vz_correct = False
hr_diff = 0

# balloon parameters
balloon = {}
balloon['altitude_step'] = 100.0 # (meters); (~100-800)
balloon['equip_mass'] = 1.608 # kg
balloon['balloon_mass'] = 1.50 # kg
balloon['fill_radius'] = 2.122/2 # meters
balloon['radius_empty'] = 2.41/2 # meters (flaccid body length - neck length)
balloon['burst_radius'] = 9.44/2 # meters
balloon['thickness_empty'] = 0.2 * 10**-3 # mm -> meters
balloon['Cd_balloon'] = 0.5
balloon['simple_ascent_rate'] = 5.0 # m/s

# parachute parameters
balloon['Cd_parachute'] = 0.97
balloon['parachute_areas'] = [1.62]
# radius = 0.5
# balloon['parachute_areas'] = np.pi * np.array([radius])**2 # m^2
# balloon['parachute_areas'] = 2*np.array([radius])**2 # m^2
balloon['parachute_change_altitude'] = None # meters

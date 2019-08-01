"""File with relevant running parameters"""

# general path
path = '/home/ellen/Desktop/SuperBIT_DRS/'

weather_data_folder = 'Weather_data/'
GFS_folder =  'GFS/'
GEFS_folder = 'GEFS/'
elevation_data_folder = 'SRTM_data/'

output_folder = 'Output/'
kml_folder = 'Kml_files/'
traj_folder = 'Trajectories/'

fig_folder = 'Figs/'
check_figs_folder = 'Checks/'

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
balloon['altitude_step'] = 100.0 # (meters); (~50-800)
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
balloon['parachute_areas'] = [1.62] # meters
balloon['parachute_change_altitude'] = None # meters

# constants
g_0 = 9.80665 # m/s surface acc.
R0 = 8.3144621 # Ideal gas constant, J/(mol*K)
R_e = 6371009 # Mean Earth radius in meters
M_air = 0.0289644 # Molar mass of air [kg/mol], altitude dependence
T0 = 288.15 # K
M_helium = 4.002602 # molar mass of helium [kg/mol], altitude dependence
Cd_sphere = 0.47 # Drag coefficient for a sphere

# login for 0.25 degrees resolution weather data
email = 'ellen.l.sirks@durham.ac.uk'
password = 'Petten36'
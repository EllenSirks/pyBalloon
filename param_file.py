"""File with relevant running parameters"""

# paths
path = '/path/to/location/of/pyBalloon/folder/'

# folders
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
next_point = '1' # 0 if not descent_only
time_interpolate = True
grid_interpolate = True
drift_time = 0. # minutes
resolution = 0.5 # 0.25 or 0.5
vz_correct = False
hr_diff = 0 # hrs, multiples of 6
check_sigmas = False

# balloon parameters
balloon = {}
balloon['altitude_step'] = 100.0 # (meters); (~50-800)
balloon['equip_mass'] = 1.608 # kg 1
balloon['balloon_mass'] = 1.50 # kg
balloon['fill_radius'] = 2.122/2 # meters
balloon['radius_empty'] = 2.41/2 # meters (flaccid body length - neck length)
balloon['burst_radius'] = 11.2/2 # meters 9.44
balloon['thickness_empty'] = 0.2 * 10**-3 # mm -> meters
balloon['Cd_balloon'] = 0.5
balloon['simple_ascent_rate'] = 5.0 # m/s

# parachute parameters
balloon['Cd_parachute'] = (0.97*1.4)/((0.57*1.75)**2) #  1.36
balloon['parachute_areas'] = [(0.57*1.75)**2] # meters [0.57**2]
balloon['parachute_change_altitude'] = None # meters

# constants
g_0 = 9.80665 # m/s**2 surface acc.
R0 = 8.3144621 # Ideal gas constant, J/(mol*K)
R_e = 6371009 # Mean Earth radius in meters
M_air = 0.0289644 # Molar mass of air [kg/mol], altitude dependence not used here
T0 = 288.15 # K
M_helium = 4.002602 # molar mass of helium [kg/mol], altitude dependence not used here
Cd_sphere = 0.47 # Drag coefficient for a sphere

# data tile size
tile_size = 6.

# login for 0.25 degrees resolution weather data
email = ''
password = ''
"""
Script to be run to initialise paths before running pyBalloon.
"""

import sys, os

import param_file as p

print('Have you changed the path and folder names in param_file.py?')
answer = raw_input('y or n: ')

if answer == 'y':

	if not os.path.exists(p.path + p.weather_data_folder + p.GFS_folder):
		os.mkdir(p.path + p.weather_data_folder + p.GFS_folder)
	if not os.path.exists(p.path + p.elevation_data_folder):
		os.mkdir(p.path + p.elevation_data_folder)
	if not os.path.exists(p.path + p.output_folder):
		os.mkdir(p.path + p.output_folder)

else:
	print('Please go to param_file.py and change the path to where the pyBalloon folder is located.')
	sys.exit()





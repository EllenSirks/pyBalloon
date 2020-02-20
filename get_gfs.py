"""
Functions for finding and downloading GFS files.
"""

from shutil import copyfile
import datetime as dt
import numpy as np
import subprocess
import requests
import sys, os
import math
import time
import re

import param_file as p

#################################################################################################################

def check_connection_device(device_name):

	device_re = re.compile("Bus\s+(?P<bus>\d+)\s+Device\s+(?P<device>\d+).+ID\s(?P<id>\w+:\w+)\s(?P<tag>.+)$", re.I)
	df = subprocess.check_output("lsusb")
	devices = []
	for i in df.split('\n'.encode()):
		if i:
			info = device_re.match(i.decode('utf-8'))
			if info:
				dinfo = info.groupdict()
				dinfo['device'] = '/dev/bus/usb/%s/%s' % (dinfo.pop('bus'), dinfo.pop('device'))
				devices.append(dinfo)

	tags = [device['tag'] for device in devices]

	return device_name in tags

#################################################################################################################

def get_latest_gfs_file(resolution=0.5):
	"""
	Find and download latest weather file available

	Arguments
	=========
	resolution : float
		Resolution of weather forecast model (GFS)
	Return:
		Weather file name in list ([])
	"""

	res1 = 'pgrb2.0p' + str(int(resolution*100))
	res2 = str(int(-4*resolution + 6))

	now = dt.datetime.utcnow()
	now_hour = int(round(now.hour - 1 + now.minute/60.)) ### using UTC time
	now_date_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

	out_dir = p.path + p.weather_data_folder + p.GFS_folder
	url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'

	############################################################################################################

	req = requests.get(url_base)
	if req.ok is not True:
		print("Could not connect!")
		sys.exit()

	text = req.content.split('</a>')
	available_folders = [t.split('>')[-1] for t in text if 'gfs.20' in t.split('>')[-1]]
	available_folders.sort()
	url_base += available_folders[-1]

	req = requests.get(url_base)
	text = req.content.split('</a>')
	available_times = [t.split('>')[-1] for t in text if '/' in t.split('>')[-1]]
	available_times.sort()

	req = requests.get(url_base + available_times[-1])
	text = req.content.split('</a>')
	available_files = [t.split('>')[-1] for t in text if res1 in t.split('>')[-1] and 'idx' not in t.split('>')[-1] and 'anl' not in t.split('>')[-1]]

	hhhh = int(available_times[-1][-3:-1])

	if len(available_files) == 0:

		req = requests.get(url_base + available_times[-2])
		text = req.content.split('</a>')
		available_files = [t.split('>')[-1] for t in text if res1 in t.split('>')[-1] and 'idx' not in t.split('>')[-1] and 'anl' not in t.split('>')[-1]]

		hhhh = int(available_times[-2][-3:-1])

	diffs = [np.abs(hhhh + int(file[-3:]) - now_hour) for file in available_files]
	weather_file = available_files[int(np.where(diffs == min(diffs))[0])]
	weather_files = ['gfs_' + res2 + '_' + str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2) + '_' + weather_file[5:7] + '00_' + weather_file[-3:] + '.grb2']

	file = get_gfs_files(weather_files=weather_files)

	return weather_files

#################################################################################################################

def get_interpolation_gfs_files(datestr=None, utc_hour=None, resolution=0.5, hr_diff=0, live=True):
	"""
	Find and download weather forecast files to be used for interpolation

	Arguments
	=========
	datestr : string
		The current date of the trajectory (yyyymmdd)
	utc_hour : float
		Current time of trajectory
	resolution : float
		Resolution of weather forecast model (GFS)
	hr_diff : float
		Number of hours in the past we wish the forecasts to be from (e.g. when we wish to download them in the morning for later in the day) needs to be a multiple of 6.
	Return:
		List of weather_file names (without .grb2 at the end)		
	"""

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Finding interpolation files...'.ljust(60) + '\r')
	sys.stdout.flush()
	time.sleep(0.2)

	res = str(int(-4*resolution + 6))

	add = 0
	if live and np.abs(get_interval(datestr=datestr, utc_hour=utc_hour)[0] - utc_hour) < 3.5: # upload lag
		add = 6

	left_hr, left_hhh, right_hr, right_hhh, left_datestr, right_datestr = get_interval(datestr=datestr, utc_hour=utc_hour, hr_diff=hr_diff+add, live=live)
	weather_files = ['gfs_' + res + '_' + left_datestr + '_' + str(left_hr*100).zfill(4) + '_' + str(left_hhh).zfill(3) + '.grb2', \
	'gfs_' + res + '_' + right_datestr + '_' + str(right_hr*100).zfill(4) + '_' + str(right_hhh).zfill(3) + '.grb2']

	file = get_gfs_files(weather_files=weather_files)

	return [weather_files[i][:-5] for i in range(len(weather_files))]

#################################################################################################################

def get_closest_gfs_file(datestr=None, utc_hour=None, resolution=0.5, hr_diff=0):
	"""
	Find & download weather_file nearest in time to current time of trajectory

	Arguments
	=========
	datestr : string
		The current date of the trajectory (yyyymmdd)
	utc_hour : float
		Current time of trajectory
	resolution : float
		Resolution of weather forecast model (GFS)
	hr_diff : float
		Number of hours in the past we wish the forecasts to be from (e.g. when we wish to download them in the morning for later in the day) needs to be a multiple of 6.
	Return:
		Weather file name (without .grb2 at the end) in list ([])
	"""

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Finding closest gfs file...'.ljust(60) + '\r')
	sys.stdout.flush()
	time.sleep(0.1)

	res = str(int((-4*resolution + 6)))

	out_dir = p.path + p.weather_data_folder + p.GFS_folder
	now = dt.datetime.now()

	hrs = get_closest_hr(datestr=datestr, utc_hour=utc_hour, hr_diff=hr_diff)
	h1, h2, h3, datestr = hrs[0], hrs[1], hrs[2], hrs[3]

	month = datestr[0:6]
	day = datestr

	hhhh = str(h1*100).zfill(4)
	hhh = str(h2).zfill(3)

	files = ['gfs_' + res + '_' + datestr + '_' + hhhh + '_' + hhh + '.grb2']
	file = get_gfs_files(weather_files=files)

	return [files[0][:-5]]

#################################################################################################################

def get_gfs_files(weather_files=None):
	"""
	Download specified weather forecasts

	Arguments
	=========
	weather_files : list
		List containing weather file from which we wish to use the data. File names should be formatted as: gfs_x_datestr_hhhh_hhh.grb2
	Return:
		List of weather file names (without .grb2 at the end)
	"""

	out_dir = p.path + p.weather_data_folder + p.GFS_folder
	usb_path = p.usb_path

	device_connected = check_connection_device(p.usb_device_name)

	now = dt.datetime.utcnow()

	for file in weather_files:

		res1 = '0p' + str(int(100*(int(file[4]) - 6)/-4.))
		res2 = file[4]

		datestr = file[6:14]
		year = datestr[:4]
		month = datestr[:6]
		hhhh, hhh = file[-13:-9], file[-8:-5]

		if datestr in [str(now.year) + str(now.month).zfill(2) + str(now.day - i).zfill(2) for i in range(0, 10)]:
			url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
			filename = url_base + 'gfs.' + datestr + '/' + str(int(hhhh)/100).zfill(2) + '/' + 'gfs.t' + str(int(hhhh)/100).zfill(2) +  'z.pgrb2.' + res1 + '.f' + hhh

		else:
			if res2 == '4':
				url_base = 'https://nomads.ncdc.noaa.gov/data/gfs4/'
				filename = url_base + month + '/' + datestr + '/' + file
			elif res2 == '5':
				url_base = 'http://rda.ucar.edu/data/ds084.1/'
				filename = url_base + year + '/' + datestr + '/gfs.0p25.' + datestr + str(int(int(hhhh)/100)).zfill(2) + '.f' + hhh + '.grib2'
			
		second_file = 'gfs_' + res2 + '_' + datestr + '_' + hhhh + '_' + hhh + '.grb2'

		if os.path.isfile(out_dir + second_file) or (device_connected and os.path.isfile(usb_path + second_file)):
			# copyfile(usb_path + second_file, out_dir + second_file)
			continue
		else:
			download_file(path_file=filename, out_dir=out_dir)

			diff_name = os.path.basename(filename) != second_file

			if device_connected:
				copyfile(out_dir + os.path.basename(filename), usb_path + second_file)
				if diff_name:
					copyfile(out_dir + os.path.basename(filename), usb_path + second_file)
					os.remove(usb_path + os.path.basename(filename))

			if diff_name:
				copyfile(out_dir + os.path.basename(filename), out_dir + second_file)
				os.remove(out_dir + os.path.basename(filename))

	return [weather_files[0][:-5]]

#################################################################################################################

def download_file(path_file=None, out_dir=None):
	"""
	Method to download 0.25 resolution weather forecast data

	Arguments
	=========
	path_file : string
		String containing the path to the file and the file name
	out_dir : strings
		Directory file is to be saved in
	"""

	file = os.path.basename(path_file)

	if 'gfs' in path_file and ('0p25' in file or '_5_' in file):
		payload = {'email': p.email, 'passwd' : p.password, 'do': 'login'}
		ret = requests.post(url = 'https://rda.ucar.edu/cgi-bin/login', data=payload)
		cookies = ret.cookies
	else:
		cookies = None

	req = requests.get(path_file, cookies=cookies, allow_redirects=True, stream=True)

	if req.ok is not True:
		print('Cannot download ' + str(file) + ' this way! Check spelling or go to: https://www.ncdc.noaa.gov/has/HAS.DsSelect')
		# sys.exit()

	filesize = int(req.headers['Content-length'])
	with open(out_dir + file, 'wb') as outfile:
		chunk_size = 1048576
		for chunk in req.iter_content(chunk_size=chunk_size):
			outfile.write(chunk)
			if chunk_size < filesize:
				check_file_status(out_dir + file, filesize)
	check_file_status(out_dir + file, filesize)

#################################################################################################################

def check_file_status(filepath=None, filesize=None):
	"""
	Check what percentage of a file has been downloaded

	Arguments
	=========
	filepath : string
		String containing the path to the file and the file name
	filesize : float
		Total size of file to be downloaded
	"""

	sys.stdout.write('\r')
	sys.stdout.flush()
	size = float(os.stat(filepath).st_size)
	percent_complete = (size/filesize)*100.
	sys.stdout.write(('Downloading ' + os.path.basename(filepath) + ', %.1f %s ' % (percent_complete, '% Completed')).ljust(60) + '\r')
	sys.stdout.flush()

#################################################################################################################

def get_closest_hr(datestr=None, utc_hour=None, hr_diff=0): 
	"""
	Find time & date that is nearest to the current time of trajectory

	Arguments
	=========
	datestr : string
		The current date of the trajectory (yyyymmdd)
	utc_hour : float
		Current time of trajectory
	hr_diff : float
		Number of hours in the past we wish the forecasts to be from (e.g. when we wish to download them in the morning for later in the day) needs to be a multiple of 6.
	Return:
		Time and date of determined closest hour
	"""

	hrs_6 = [0., 6., 12., 18.]

	if float(utc_hour) in hrs_6 and hr_diff == 0:

		closest_model = int(utc_hour)
		hhh3 = 0
		hhh6 = 0

	else:

		for hr in hrs_6:
			if hr <= float(utc_hour):
				closest_model = int(hr)

		diffs = np.array([np.abs(utc_hour - (closest_model + 3.*i)) for i in range(0, 129)])
		hhh3 = 3*np.where(diffs == min(diffs))[0][0]

		diffs = np.array([np.abs(utc_hour - (closest_model + 6.*i)) for i in range(0, 65)])
		hhh6 = 6*np.where(diffs == min(diffs))[0][0]

	closest_model -= hr_diff
	hhh3 += hr_diff
	days = 0

	if closest_model < 0:

		days = 0
		while closest_model < 0:
			closest_model += 24
			days += 1	

		year, month, day = int(datestr[:4]), int(datestr[4:6]), int(datestr[6:])
		date = dt.datetime(year, month, day) - dt.timedelta(days=days)
		datestr = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)

	return (closest_model, hhh3, hhh6, datestr)

#################################################################################################################

def get_interval(datestr=None, utc_hour=None, hr_diff=0, live=True):
	"""
	Find closest time & date left & right of current time of trajectory

	Arguments
	=========
	datestr : string
		The current date of the trajectory (yyyymmdd)
	utc_hour : float
		Current time of trajectory
	hr_diff : float
		Number of hours in the past we wish the forecasts to be from (e.g. when we wish to download them in the morning for later in the day) needs to be a multiple of 6
	Return:
		Times and dates of determined interval
	"""

	times = np.array([0., 6., 12., 18., 24.])

	if utc_hour in times and hr_diff == 0:
		left_hr, left_hhh, left_hhh6, left_datestr = get_closest_hr(datestr=datestr, utc_hour=utc_hour)
		right_hr, right_hhh, right_datestr = left_hr, left_hhh + 3, left_datestr
	else:
		times = np.array([0., 3., 6., 9., 12., 15., 18., 21., 24.])

		for i in range(len(times)):
			if utc_hour >= times[i] and times[i+1] > utc_hour:
				left_hour, right_hour = times[i], times[i+1]
				break

		left_hr, left_hhh, left_hhh6, left_datestr = get_closest_hr(datestr=datestr, utc_hour=left_hour)
		right_hr, right_hhh, right_hhh6, right_datestr = get_closest_hr(datestr=datestr, utc_hour=right_hour)

		# Comment for best possible weather forecasts!
		if live and right_hr == left_hr + 6:
			right_hr -= 6
			right_hhh += 6

		if hr_diff > 0:

			left_hr -= hr_diff
			right_hr -= hr_diff

			days1 = 0
			while left_hr < 0:
				left_hr += 24
				days1 += 1	

			days2 = 0
			while right_hr < 0:
				right_hr += 24
				days2 += 1	

			if live:
				right_hr = left_hr

			left_hhh += hr_diff
			right_hhh += hr_diff

			year, month, day = int(left_datestr[:4]), int(left_datestr[4:6]), int(left_datestr[6:])
			date1 = dt.datetime(year, month, day) - dt.timedelta(days=days1)

			year, month, day = int(right_datestr[:4]), int(right_datestr[4:6]), int(right_datestr[6:])
			date2 = dt.datetime(year, month, day) - dt.timedelta(days=days2)

			left_datestr = str(date1.year) + str(date1.month).zfill(2) + str(date1.day).zfill(2)
			right_datestr = str(date2.year) + str(date2.month).zfill(2) + str(date2.day).zfill(2)

	return left_hr, left_hhh, right_hr, right_hhh, left_datestr, right_datestr

#################################################################################################################

if __name__ == '__main__':

	# file = sys.argv[1]
	# get_gfs_files(weather_files=[file])
	print(get_interval(datestr='20180406', utc_hour=9.0, hr_diff=0, live=True))

#################################################################################################################
#!/usr/bin/python

from shutil import copyfile
import datetime as dt
import numpy as np
import requests
import time
import wget
import sys
import os

def nearest_hr(utc_hour=None):

	utc_hour = float(utc_hour)

	hrs = [0., 3., 6., 9., 12., 15., 18., 21., 24.]

	if utc_hour in hrs:
		return utc_hour

	diffs_abs = np.array([np.abs(hr - utc_hour) for hr in hrs])
	diffs = np.array([hr - utc_hour for hr in hrs])

	index = int(np.where(diffs_abs == min(diffs_abs))[0][0]) # if hr exactly between two use past hr / hr closest to past model
	hr = hrs[index]

	return hr

# method to find time & date left & right of time of ascent/descent
def get_interval(datestr=None, utc_hour=None):

	utc_hour = float(utc_hour)
	hhhh = [0., 6., 12., 18., 24.]

	if utc_hour in hhhh:
		min_hr, min_hhh = utc_hour, 0.

	else:

		hhh = [3.*i for i in range(0, 129)]
		for hr in hhhh:
			if hr < utc_hour:
				min_hr = hr

		diffs = [np.abs(min_hr + h - utc_hour) for h in hhh]
		min_hhh = hhh[int(np.where(diffs == min(diffs))[0])]

		if min_hr + min_hhh > utc_hour:
			min_hhh -= min_hhh - 3.

	return int(min_hr), int(min_hhh)

# method to find time & date that is nearest to the time of ascent/descent
def get_closest_hr(utc_hour=None):

	utc_hour = float(utc_hour)
	hr = nearest_hr(utc_hour)

	hhhh = [0., 6., 12., 18., 24.]

	if utc_hour in hhhh:
		return (int(utc_hour), 0)

	for hh in hhhh:
		if hh < utc_hour:
			min_hr = int(hh)

	diffs = np.array([np.abs(min_hr - utc_hour + 3.*i) for i in range(0, 129)])
	hhh = 3*np.where(diffs == min(diffs))[0][0]

	# if hrs[int(min_ind)] == 24.:
	# 	hrs[int(min_ind)] = 0.
	# 	datestr = date_check(datestr=datestr)

	return (min_hr, hhh)

# method to iterate date by 1 day and return new datestr
def date_check(datestr=None):

	year, month, day = int(datestr[:4]), int(datestr[4:6]), int(datestr[6:])
	date = dt.datetime(year, month, day) + dt.timedelta(days=1)
	new_datestr = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)

	return new_datestr

# method to download weather file if entire file string is given
def get_gfs_file(weather_file=None, verbose=False):

	out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

	if weather_file == None:

		print('Need a weather file to download!')
		sys.exit()

	else:

		print('Getting ' + weather_file + '...')

		if os.path.isfile(out_dir + weather_file):

			print(weather_file + ' already downloaded!')
			return [weather_file[:-5]]
			
		else:

			if 'gfs.' in weather_file:

				hr = weather_file[5:7]
				now = dt.datetime.now()
				datestr = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

				url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
				url = url_base + 'gfs.' + datestr + hr + '/' +  weather_file

				print('Downloading ' + url + '...')
				filename = wget.download(url, out=out_dir)

				second_file = 'gfs_4_' + str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2) + '_' + weather_file[5:7] + '00_' + weather_file[-3:] + '.grb2'
				copyfile('/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + weather_file, '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + second_file)
				os.remove(out_dir + weather_file)

				weather_file = second_file

			else:

				url_base = 'https://nomads.ncdc.noaa.gov/data/gfs4/'

				datestr = weather_file[6:14]
				month = datestr[0:6]

				url = url_base + month
				req = requests.get(url)

				if req.ok is not True:
					print("Could not connect! Year or month not available.")
					sys.exit()

				day = datestr

				url = url_base + month + '/' + day
				req = requests.get(url)

				if req.ok is not True:
					print("Could not connect! Day not available.")
					sys.exit()

				url = url_base + month + '/' + day + '/' + weather_file

				print('Downloading ' + weather_file + '...')
				filename = wget.download(url, out=out_dir)

	return [weather_file[:-5]]

# method to find & download weather_file nearest in time to time of ascent/descent
def get_closest_gfs_file(datestr=None, utc_hour=None, verbose=False):

	print('Finding closest gfs file...')

	url_base = 'https://nomads.ncdc.noaa.gov/data/gfs4/'
	out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

	hhhh, hhh = get_closest_hr(utc_hour=utc_hour)

	month = datestr[0:6]
	day = datestr

	hhhh = str(hhhh*100).zfill(4)
	hhh = str(hhh).zfill(3)

	weather_file = 'gfs_4_' + datestr + '_' + hhhh + '_' + hhh + '.grb2'

	if os.path.isfile('/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + weather_file):
		print(weather_file + ' already downloaded!')
		return [weather_file[:-5]]

	url = url_base + month 
	req = requests.get(url)

	if req.ok is not True:
		print("Could not connect! Year or month not available.")
		sys.exit()

	url = url_base + month + '/' + day
	req = requests.get(url)

	if req.ok is not True:
		print("Could not connect! Day not available.")
		sys.exit()

	text = req.content.split('</a>')
	available_files = []

	for t in text:
		if 'grb2' in t and not 'align' in t:
			available_files.append(t.split('>')[-1].split('>')[-1])

	if weather_file not in available_files:
		print("File not available.")
		sys.exit()
	else:

		url = url_base + month + '/' + day + '/' + weather_file

		print('Downloading ' + weather_file + '...')
		filename = wget.download(url, out=out_dir)

		return [weather_file[:-5]]

# method to find & download latest weather file
def get_latest_gfs_file(verbose=False):

	now = dt.datetime.now()
	now_hour = int(round(now.hour - 1 + now.minute/60.)) ### using UTC time

	url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
	out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

	req = requests.get(url_base)
	if req.ok is not True:
		print("Could not connect!")
		sys.exit()

	text = req.content.split('</a>')
	available_folders = []

	for t in text:

		t = t.split('>')[-1]
		if 'gfs.20' in t:
			available_folders.append(t)

	available_folders.sort()

	url = url_base + available_folders[-1]
	hhhh = url[-3:-1]

	req = requests.get(url)
	if req.ok is not True:
		print("Could not connect!")
		sys.exit()

	text = req.content.split('</a>')
	available_files = []

	for t in text:
		t = t.split('>')[-1]
		if 'pgrb2.0p50' in t and 'idx' not in t and 'anl' not in t:
			available_files.append(t)

	hhhh = int(url[-3:-1])
	diffs = [np.abs(hhhh + int(file[-3:]) - now_hour)  for file in available_files]

	weather_file = available_files[int(np.where(diffs == min(diffs))[0])]
		
	if not os.path.isfile('/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + weather_file):

		url = url + weather_file
		print('Downloading ' + url + '...')
		filename = wget.download(url, out=out_dir)

		second_file = 'gfs_4_' + str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2) + '_' + weather_file[5:7] + '00_' + weather_file[-3:] + '.grb2'
		copyfile('/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + weather_file, '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + second_file)
		os.remove(out_dir + weather_file)

	else:
		print(weather_file + ' already downloaded!')

	return [weather_file]

# method to find & download weather files needed for interpolation (3)
def get_interpolation_gfs_files(weather_file=None, datestr=None, utc_hour=None):

	time0 = time.time()
	now = dt.datetime.now()

	print('Getting interpolation files...')

	out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

	if weather_file == None:
		hhhh, hhh = get_interval(datestr=datestr, utc_hour=utc_hour)
	else:
		if 'gfs.'in weather_file:
			datestr, hhhh, hhh = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2), int(weather_file[5:7]), int(weather_file[-3:])
		else:
			datestr, hhhh, hhh = weather_file[6:14], int(int(weather_file[15:19])/100.), int(weather_file[20:23])

	files = []

	i = 0
	while i < 3:

		if datestr == str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2):
			weather_file = 'gfs.t' + str(hhhh).zfill(2) + 'z.pgrb2.0p50.f' + str(hhh).zfill(3)
		else: 
			weather_file = 'gfs_4_' + datestr + '_' + str(hhhh*100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2'

		file = get_gfs_file(weather_file=weather_file)[0]

		files.append(file)

		hhh += 3
		i+=1

	return files

if __name__ == '__main__':

	time0 = time.time()

	res = get_latest_gfs_file()
	print(res)

	print('Program finished in %.1f s' % (time.time() - time0))
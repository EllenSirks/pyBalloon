#!/usr/bin/python
from shutil import copyfile
import datetime as dt
import numpy as np
import requests
import time
import wget
import sys
import os

def get_gfs_file(weather_file=None, verbose=False):

	if weather_file == None:

		print('Need a weather file to download!')
		sys.exit()

	else:

		if os.path.isfile('/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + weather_file):

			print('File already downloaded!')
			return weather_file[15:19], weather_file[20:23]

		else:

			url_base = 'https://nomads.ncdc.noaa.gov/data/gfs4/'
			out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

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

			hhhh = weather_file[15:19]
			hhh = weather_file[20:23]

			return hhhh, hhh

def get_closest_gfs_file(datestr=None, utc_hour=None, verbose=False):

	url_base = 'https://nomads.ncdc.noaa.gov/data/gfs4/'
	out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

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

	text = req.content.split('</a>')
	available_files = []

	for t in text:
		if 'grb2' in t and not 'align' in t:
			available_files.append(t.split('>')[-1].split('>')[-1])

	check = 0
	for file in available_files:

		if check == 1:
			break

		yyyy_init, mm_init, dd_init, hhhh_init, hhh_init = int(file[6:10]), int(file[10:12]), int(file[12:14]), int(int(file[15:19])/100) , int(file[20:23])

		yyyy, mm, dd = int(datestr[0:4]), int(datestr[4:6]), int(datestr[6:8])
		request_time = dt.datetime(yyyy, mm, dd, int(utc_hour))

		init_dt = dt.datetime(yyyy_init, mm_init, dd_init, hhhh_init)
		delta_t = request_time - init_dt

		if delta_t.seconds < 0 or delta_t.days < 0 or np.abs(hhhh_init - int(utc_hour)) > 5 or np.abs(hhhh_init + hhh_init - int(utc_hour)) > 2:
			continue

		check += 1

		weather_file = file

	url = url_base + month + '/' + day + '/' + weather_file

	if not os.path.isfile('/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + weather_file):

		print('Downloading ' + weather_file + '...')
		filename = wget.download(url, out=out_dir)

	else:
		print('File already downloaded!')

	hhhh = weather_file[15:19]
	hhh = weather_file[20:23]

	return hhhh, hhh

def get_latest_gfs(verbose=False):

	now = dt.datetime.now()
	now_hour = int(round(now.hour - 5 + now.minute/60.)) ### check what time is being used!!

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
		print('Downloading ' + weather_file + '...')
		filename = wget.download(url, out=out_dir)

		if now.month < 10:
			now_month = '0' + str(now.month)
		else:
			now_month = str(now.month)

		if now.day < 10:
			now_day = '0' + str(now.day)
		else:
			now_day = str(now.day)

		second_file = 'gfs_4_' + str(now.year) + now_month + now_day + '_' + weather_file[5:7] + '00_' + weather_file[-3:] + '.grb2'
		copyfile('/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + weather_file, '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + second_file)

	else:
		print(weather_file + ' already downloaded!')

def get_future_gfs_files(ini_weather_file=None, no_files=1):

	print('Getting ' + str(no_files) + ' prediction(s)...')

	out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/'

	if 'gfs_4' in ini_weather_file:
		url_base = 'https://nomads.ncdc.noaa.gov/data/gfs4/' + ini_weather_file[6:12] + '/' + ini_weather_file[6:14] + '/'
	else:
		url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/' + ini_weather_file[6:12] + '/' + ini_weather_file[6:14] + '/'

	hhh = int(ini_weather_file[20:23])

	i = 0
	while i < no_files:

		hhh_next = hhh + 3

		if hhh_next < 10:
			hhh_str = '00' + str(hhh_next)
		elif hhh_next > 9 and hhh_next< 100:
			hhh_str = '0' + str(hhh_next)
		else:
			hhh_str = str(hhh_next)

		weather_file = ini_weather_file.replace(ini_weather_file[20:23], hhh_str)
		hhh = hhh_next

		if not os.path.isfile('/home/ellen/Desktop/SuperBIT/Weather_data/grb_files/' + weather_file):

			url = url_base + weather_file
		
			print('Downloading ' + weather_file + '...')
			filename = wget.download(url, out=out_dir)

		else:

			print(weather_file + ' already downloaded!')

		i+=1

if __name__ == '__main__':

	time0 = time.time()

	get_latest_gfs()
	# get_future_gfs_files(ini_weather_file='gfs_4_20180304_0600_000.grb2', no_files=5)

	print('Program finished in %.1f s' % (time.time() - time0))
#!/usr/bin/python

from shutil import copyfile
import datetime as dt
import numpy as np
import requests
import time
import wget
import sys
import os

import get_gfs

# method to download weather file if entire file string is given
def get_gfs_file(weather_files=None, file_type='GFS', verbose=False):

	out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/' + file_type + '/'

	now = dt.datetime.now()

	if weather_files == None:
		print('Need weather files to download!')
		sys.exit()

	if file_type == 'GFS':

		datestr = weather_files[0][6:14]
		month = datestr[0:6]
		day = datestr
		hhhh, hhh = weather_files[0][-13:-9], weather_files[0][-8:-5]

		if datestr in [str(now.year) + str(now.month).zfill(2) + str(now.day - i).zfill(2) for i in range(0, 10)]:
			url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
			url = url_base + 'gfs.' + datestr + hr + '/'
			weather_files_names = ['gfs.t' + str(int(hhhh)/100).zfill(2) +  'z.pgrb2.0p50.f' + hhh]
		else:
			url_base = 'https://nomads.ncdc.noaa.gov/data/gfs4/'
			url = url_base + month + '/' + day + '/'
			weather_files_names = weather_files

	elif file_type == 'GEFS':

		datestr = weather_files[0][11:19]
		hhhh, hhh = weather_files[0][-13:-9], weather_files[0][-8:-5]

		out_dir += datestr + '/'

		if datestr in [str(now.year) + str(now.month).zfill(2) + str(now.day - i).zfill(2) for i in range(0, 7)]:
			url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gens/prod/'
			url = url_base + 'gefs.' + datestr + '/' + str(int(hhhh)/100).zfill(2) + '/pgrb2ap5/'
			weather_files_names = ['geavg.t' + str(int(hhhh)/100).zfill(2) +  'z.pgrb2a.0p50.f' + hhh, 'gespr.t' + str(int(hhhh)/100).zfill(2) +  'z.pgrb2a.0p50.f' + hhh]
		else:
			print('Cannot download these GEFS files via ftp!')
			sys.exit()

	for i in range(len(weather_files)):

		print('Getting ' + weather_files[i] + '...')

		if os.path.isfile(out_dir + weather_files[i]):

			print(weather_files[i] + ' already downloaded!')
			return [weather_files[i][:-5]]
				
		else:

			print('Downloading ' + url + weather_files[i] + '...')
			filename = wget.download(url + weather_files_names[i], out=out_dir)

			second_file = weather_files[i]
			copyfile(out_dir + weather_files_names[i], out_dir + second_file)
			os.remove(out_dir + weather_files_names[i])

	return [weather_files[0][:-5]]

# method to find & download weather_file nearest in time to time of ascent/descent
def get_closest_gfs_file(datestr=None, utc_hour=None, file_type='GFS', verbose=False):

	print('Finding closest ' + file_type.lower() + ' file...')

	out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/' + file_type + '/'
	now = dt.datetime.now()

	h1, h2, h3 = get_gfs.get_closest_hr(utc_hour=utc_hour)

	month = datestr[0:6]
	day = datestr

	hhhh = str(h1*100).zfill(4)
	hhh1 = str(h2).zfill(3)
	hhh2 = str(h3).zfill(3)

	time_str1 = '_' + datestr + '_' + hhhh + '_' + hhh1 + '.grb2'
	time_str2 = '_' + datestr + '_' + hhhh + '_' + hhh2 + '.grb2'

	if file_type == 'GEFS' or file_type == 'gefs': 
		out_dir += datestr + '/'
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)
		ind = 7
	elif file_type == 'GFS' or file_type == 'gfs':
		ind = 10

	if datestr in [str(now.year) + str(now.month).zfill(2) + str(now.day - i).zfill(2) for i in range(0, ind+1)]:

		url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/'

		if file_type == 'GFS':
			url_base += 'gfs/prod/'
			weather_files1 = 'gfs_4' + time_str1

		elif file_type == 'GEFS': 
			url_base += 'gens/prod/'
			weather_files1 = ['geavg_0p50' + time_str2,  'gespr_0p50' + time_str2]

		if os.path.isfile(out_dir + weather_files1[0]):
			print(weather_files1[0] + ' already downloaded!')
			return [weather_files1[0][:-5]]

		url1 = url_base + file_type.lower() + '.' + datestr 

		if file_type == 'GFS':
			url1 += str(h1).zfill(2) + '/' 
			weather_files = ['gfs.t' + str(h1).zfill(2) + 'z.pgrb2.0p50.f' + hhh1]
		elif file_type == 'GEFS':
			url1 += '/' + str(h1).zfill(2) + '/pgrb2ap5/'
			weather_files = ['geavg.t' + str(h1).zfill(2) + 'z.pgrb2a.0p50.f' + hhh2, 'gespr.t' + str(int(hhhh)/100).zfill(2) +  'z.pgrb2a.0p50.f' + hhh2]

		req = requests.get(url1)
		if req.ok is not True:
			print("Could not connect! Month/day/hour not available.")
			sys.exit()

		text = req.content.split('</a>')
		available_files = []

		for t in text:
			t = t.split('>')[-1]
			if '.0p50' in t and 'pgrb2' in t and 'idx' not in t and 'anl' not in t:
				available_files.append(t)

		for i in range(len(weather_files)):
			if weather_files[i] not in available_files:
				print("File not available.")
				sys.exit()
			else:
				url2 = url1 + weather_files[i]
				print('Downloading ' + url1 + weather_files[i] + '...')
				filename = wget.download(url2, out=out_dir)

				second_file = weather_files1[i]
				copyfile(out_dir + weather_files[i], out_dir + second_file)
				os.remove(out_dir + weather_files[i])

		return [weather_files1[0]][:-5]

	else:

		if file_type == 'GFS':

			url_base = 'https://nomads.ncdc.noaa.gov/data/gfs4/'
			weather_file = 'gfs_4_' + datestr + '_' + hhhh + '_' + hhh1 + '.grb2'

			if os.path.isfile(out_dir + weather_file):
				print(weather_file + ' already downloaded!')
				return [weather_file[:-5]]

			url = url_base + month + '/' + day
			req = requests.get(url)

			if req.ok is not True:
				print("Could not connect! Date not available.")
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

		elif file_type == 'GEFS':
			print('Cannot download these GEFS files like this!')
			sys.exit()

# method to find & download latest weather file
def get_latest_gfs_file(file_type='GFS', verbose=False):

	now = dt.datetime.now()
	now_hour = int(round(now.hour - 1 + now.minute/60.)) ### using UTC time
	now_date_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

	out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/' + file_type + '/'
	url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/'

	if file_type == 'GFS':

		url_base += 'gfs/prod/'
		req = requests.get(url_base)
		if req.ok is not True:
			print("Could not connect!")
			sys.exit()

		text = req.content.split('</a>')
		available_folders = []

		for t in text:

			t = t.split('>')[-1]
			if 'gfs.20' in t and now_date_str in t:
				available_folders.append(t)

		folder_lengths = [len(folder) for folder in available_folders]

		if len(set(folder_lengths)) == len(folder_lengths):

			available_folders.sort()

			url_base += available_folders[-1]

			if len(available_folders[-1]) == 13:

				req = requests.get(url_base)
				if req.ok is not True:
					print("Could not connect!")
					sys.exit()

				text = req.content.split('</a>')
				available_times = []

				for t in text:
					t = t.split('>')[-1]
					if '/' in t:
						available_times.append(t)

				available_times.sort()
				url_base += available_times[-1]

			hhhh = int(url_base[-3:-1])

		else:

			times = {}

			for folder in available_folders:

				times[folder] = []

				if len(folder) == 13:

					req = requests.get(url_base + folder)

					text = req.content.split('</a>')
					available_times = []

					for t in text:
						t = t.split('>')[-1]
						if '/' in t:
							times[folder].append(int(t))

				else:

					times[folder].append(int(folder[-3:-1]))

			time_list = list(times.values())
			keys = list(times.keys())

			folder = keys[np.where(time_list == max(time_list))[0][0]]

			if len(folder) == 13:
				url_base += folder + str(max(time_list)) + '/'
				hhhh = max(time_list)
			else:
				url_base += folder
				hhhh = int(url_base[-3:-1])

		req = requests.get(url_base)
		if req.ok is not True:
			print("Could not connect!")
			sys.exit()

		text = req.content.split('</a>')
		available_files = []

		for t in text:
			t = t.split('>')[-1]
			if 'pgrb2.0p50' in t and 'idx' not in t and 'anl' not in t:
				available_files.append(t)

		diffs = [np.abs(hhhh + int(file[-3:]) - now_hour) for file in available_files]

		weather_files = [available_files[int(np.where(diffs == min(diffs))[0])]]

	else:

		out_dir += now_date_str + '/'

		if not os.path.exists(out_dir):
			os.makedirs(out_dir)

		url_base += 'gens/prod/gefs.' + now_date_str + '/'

		req = requests.get(url_base)
		if req.ok is not True:
			print("Could not connect!")
			sys.exit()

		text = req.content.split('</a>')
		available_times = []

		for t in text:
			t = t.split('>')[-1]
			if '/' in t:
				available_times.append(int(t[:-1]))

		available_times.sort()

		hhhh = available_times[-1]

		url_base += str(hhhh).zfill(2) + '/pgrb2ap5/'

		req = requests.get(url_base)
		if req.ok is not True:
			print("Could not connect!")
			sys.exit()

		text = req.content.split('</a>')
		available_files = []

		for t in text:
			t = t.split('>')[-1]
			if 'pgrb2a.0p50' in t and 'idx' not in t and 'anl' not in t and ('geavg' in t or 'gespr' in t):
				available_files.append(t)

		diffs = [np.abs(hhhh + int(file[-3:]) - now_hour) for file in available_files]
		weather_files = [available_files[int(np.where(diffs == min(diffs))[0][0])], available_files[int(np.where(diffs == min(diffs))[0][1])]]

	for weather_file in weather_files:

		if file_type == 'GFS':
			second_file = 'gfs_4_' + str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2) + '_' + weather_file[5:7] + '00_' + weather_file[-3:] + '.grb2'
		elif file_type == 'GEFS':
			second_file = weather_file[:5] + '_0p50_' + str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2) + '_' + weather_file[7:9] + '00_' + weather_file[-3:] + '.grb2'

		if not os.path.isfile(out_dir + second_file):

			url = url_base + weather_file
			print('Downloading ' + url + '...')
			filename = wget.download(url, out=out_dir)

			copyfile(out_dir + weather_file, out_dir + second_file)
			os.remove(out_dir + weather_file)

		else:
			print(second_file + ' already downloaded!')

	return weather_files

if __name__ == '__main__':

	get_latest_gfs_file('GFS')
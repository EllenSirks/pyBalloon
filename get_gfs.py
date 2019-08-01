from shutil import copyfile
import datetime as dt
import numpy as np
import requests
import sys, os
import math
import time

import param_file as p

#################################################################################################################

# method to find & download latest weather file
def get_latest_gfs_file(resolution=0.5):

	res1 = 'pgrb2.0p' + str(int(resolution*100))
	res2 = str(int(-4*resolution + 6))

	now = dt.datetime.now()
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

# method to find & download weather files needed for interpolation (3)
def get_interpolation_gfs_files(datestr=None, utc_hour=None, resolution=0.5):

	print('Getting interpolation files...')

	res = str(int(-4*resolution + 6))
	left_hr, left_hhh, right_hr, right_hhh = get_interval(utc_hour=utc_hour)
	weather_files = ['gfs_' + res + '_' + datestr + '_' + str(left_hr*100).zfill(4) + '_' + str(left_hhh).zfill(3) + '.grb2', 'gfs_' + res + '_' + datestr + '_' + str(right_hr*100).zfill(4) + '_' + str(right_hhh).zfill(3) + '.grb2']
	file = get_gfs_files(weather_files=weather_files)

	return [weather_files[i][:-5] for i in range(len(weather_files))]

#################################################################################################################

# method to find & download weather_file nearest in time to time of ascent/descent
def get_closest_gfs_file(datestr=None, utc_hour=None, resolution=0.5, hr_diff=0):

	sys.stdout.write('\r')
	sys.stdout.flush()
	sys.stdout.write('Finding closest gfs file...'.ljust(60) + '\r')
	sys.stdout.flush()
	time.sleep(1)

	res = str(int((-4*resolution + 6)))

	out_dir = p.path + p.weather_data_folder + p.GFS_folder
	now = dt.datetime.now()

	hrs = get_closest_hr(datestr=datestr, utc_hour=utc_hour, hr_diff=hr_diff) ### change here!!
	h1, h2, h3, datestr = hrs[0], hrs[1], hrs[2], hrs[3]

	month = datestr[0:6]
	day = datestr

	hhhh = str(h1*100).zfill(4)
	hhh = str(h2).zfill(3)

	files = ['gfs_' + res + '_' + datestr + '_' + hhhh + '_' + hhh + '.grb2']
	file = get_gfs_files(weather_files=files)

	return [files[0][:-5]]

#################################################################################################################

def get_gfs_files(weather_files=None): # files should be formatted like; gfs_x_datestr_hhhh_hhh.grb2

	out_dir = p.path + p.weather_data_folder + p.GFS_folder
	now = dt.datetime.now()

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

		if os.path.isfile(out_dir + second_file) and float(os.stat(out_dir + second_file).st_size) > 40000000.:
			# print(file + ' already downloaded!')
			continue

		download_file(path_file=filename, out_dir=out_dir)

		if os.path.basename(filename) != second_file:
			copyfile(out_dir + os.path.basename(filename), out_dir + second_file)
			os.remove(out_dir + os.path.basename(filename))

	return [weather_files[0][:-5]]

#################################################################################################################

def get_gefs_files(datestr=None, utc_hour=None): # gfs files we wish to have the gefs for, in format: gfs_x_datestr_hhhh_hhh.grb2

	out_dir = p.path + p.weather_data_folder + p.GFS_folder
	now = dt.datetime.now()

	url_base = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gens/prod/'

	hrs = get_closest_hr(utc_hour=utc_hour)
	closest_model, hhh3, hhh6 = hrs[0], hrs[1], hrs[2]
 
	if datestr in [str(now.year) + str(now.month).zfill(2) + str(now.day - i).zfill(2) for i in range(0, 10)]:

		url = url_base + 'gefs.' + datestr + '/' + str(closest_model).zfill(2)  + '/pgrb2ap5/' 

		dl_files = ['geavg.t' + str(closest_model).zfill(2) + 'z.pgrb2a.0p50.f' + str(hhh3).zfill(3), 'gespr.t' + str(closest_model).zfill(2) + 'z.pgrb2a.0p50.f' + str(hhh3).zfill(3)]
		second_files = ['geavg_4_' + datestr + '_' + str(int(closest_model*100)).zfill(4) + '_' + str(hhh3).zfill(3) + '.grb2', 'gespr_4_' + datestr + '_' + str(int(closest_model*100)).zfill(4) + '_' + str(hhh3).zfill(3) + '.grb2'] 

		for i in range(len(dl_files)):

			if os.path.isfile(out_dir + second_files[i]) and float(os.stat(out_dir + second_files[i]).st_size) > 40000000.:
				print(second_files[i] + ' already downloaded!')
				continue

			download_file(path_file=url + dl_files[i], out_dir=out_dir)
			copyfile(out_dir + dl_files[i], out_dir + second_files[i])
			os.remove(out_dir + dl_files[i])

	else:

		second_files = ['gens-a_3_' + datestr + '_' + str(int(closest_model*100)).zfill(4) + '_' + str(hhh6).zfill(3) + '_' + str(i).zfill(2) + '.grb2' for i in range(0, 21)]

		for file in second_files:
			
			if os.path.isfile(out_dir + file) and float(os.stat(out_dir + file).st_size) > 40000000.:
				print(file + ' already downloaded!')

			else:
				print('Cannot download ' + str(file) + ' this way! Go to: https://www.ncdc.noaa.gov/has/HAS.DsSelect')
				continue

#################################################################################################################

def download_file(path_file=None, out_dir=None):

	file = os.path.basename(path_file)

	if 'gfs' in path_file and '0p25' in file:
			values = {'email' : p.email, 'passwd' : p.password, 'action' : 'login'}
			ret = requests.post(url = 'https://rda.ucar.edu/cgi-bin/login',  data = values)
			cookies = ret.cookies
	else:
		cookies = None

	req = requests.get(path_file, cookies = cookies, allow_redirects=True, stream=True)

	if req.ok is not True:
		print('Cannot download ' + str(file) + ' this way! Check spelling or go to: https://www.ncdc.noaa.gov/has/HAS.DsSelect')
		sys.exit()

	filesize = int(req.headers['Content-length'])
	with open(out_dir + file, 'wb') as outfile:
		chunk_size = 1048576
		for chunk in req.iter_content(chunk_size=chunk_size):
			outfile.write(chunk)
			if chunk_size < filesize:
				check_file_status(out_dir + file, filesize)
	check_file_status(out_dir + file, filesize)

#################################################################################################################

# method to check what percentage of a file has been downloaded
def check_file_status(filepath=None, filesize=None):

	sys.stdout.write('\r')
	sys.stdout.flush()
	size = float(os.stat(filepath).st_size)
	percent_complete = (size/filesize)*100.
	sys.stdout.write(('Downloading ' + os.path.basename(filepath) + ', %.1f %s ' % (percent_complete, '% Completed')).ljust(40) + '\r')
	sys.stdout.flush()

#################################################################################################################

# method to find time & date that is nearest to the time of ascent/descent
# hr_diff needs to be a multiple of 6
def get_closest_hr(datestr=None, utc_hour=None, hr_diff=0): 

	hrs_6 = [0., 6., 12., 18., 24.]

	if utc_hour in hrs_6 and hr_diff == 0:

		closest_model = int(utc_hour)
		hhh3 = 0
		hhh6 = 0

	else:

		for hr in hrs_6:
			if hr < utc_hour:
				closest_model = int(hr)

		diffs = np.array([np.abs(utc_hour - (closest_model + 3.*i)) for i in range(0, 129)])
		hhh3 = 3*np.where(diffs == min(diffs))[0][0]

		diffs = np.array([np.abs(utc_hour - (closest_model + 6.*i)) for i in range(0, 65)])
		hhh6 = 6*np.where(diffs == min(diffs))[0][0]

	hour = closest_model + hhh3

	closest_model -= hr_diff
	hhh3 += hr_diff
	days = 0

	if closest_model < 0:

		days = int(math.ceil((-1*closest_model) / 24.))
		closest_model %= 24

		year, month, day = int(datestr[:4]), int(datestr[4:6]), int(datestr[6:])
		date = dt.datetime(year, month, day) - dt.timedelta(days=days)
		datestr = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)

	return (closest_model, hhh3, hhh6, datestr)

#################################################################################################################

# method to find time & date left & right of time of ascent/descent
def get_interval(utc_hour=None):

	hrs = get_closest_hr(utc_hour=utc_hour)
	left_hr, left_hhh, hhh6 = hrs[0], hrs[1], hrs[2]

	if float(left_hr) == utc_hour and left_hhh == 0:
		left_hr, left_hhh = int(utc_hour), 0
	else:
		if left_hr + left_hhh > utc_hour:
			left_hhh -= 3

	right_hr, right_hhh = left_hr, left_hhh + 3

	return left_hr, left_hhh, right_hr, right_hhh

#################################################################################################################

# method to iterate date by 1 day and return new datestr
def date_check(datestr=None):

	year, month, day = int(datestr[:4]), int(datestr[4:6]), int(datestr[6:])
	date = dt.datetime(year, month, day) + dt.timedelta(days=1)
	new_datestr = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)

	return new_datestr

#################################################################################################################

if __name__ == '__main__':

	file = sys.argv[1]

	get_gfs_files(weather_files=[file])

#################################################################################################################
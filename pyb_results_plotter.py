""" method to plot final results """

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from haversine import haversine
from astropy.io import ascii
import matplotlib as mpl
import datetime as dt
import numpy as np
import time
import csv
import os

import param_file as p
import pyb_io

def plot_results(run=None, params=None):

	time0 = time.time()
	now = dt.datetime.now()
	now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)

	print('Plotting final results...')

	if params == None:
		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		interpolate = p.interpolate
		drift_time = p.drift_time

	else:

		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		interpolate = bool(params[-2])
		drift_time = float(params[-1])

	if interpolate:
		interp = '_interpolated'
	else:
		interp = ''

	dir = '/home/ellen/Desktop/SuperBIT/'
	dir_pred =  dir + 'Output/'
	dir_gps =  dir + 'Flight_data/'

	if run == None:
		files = [filename for filename in os.listdir(dir_pred) if now_str in filename]
		folder = now_str + '_' + str(len(files)-1) + '/'

	else:
		folder = run + '/'

	fig_dir = dir_pred + folder + 'figs/Results/'
	dir_pred += folder + 'trajectories/'

	# check if the paths exist/make them
	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	d = {}
	temperatures = {}
	total_time = {}
	total_distance = {}
	ini_alts = {}
	ini_lats = {}
	ini_lons = {}

	for filename in os.listdir(dir_pred):

		if filename.endswith(str(int(drift_time)).zfill(4) + 'min.dat'):

			if not interpolate and 'interpolated' in filename:
				continue
			elif interpolate and 'interpolated' not in filename:
				continue

			data_pred = ascii.read(dir_pred + filename)

			end_lat_pred = data_pred['lats'][-1]
			end_lon_pred = data_pred['lons'][-1]

			if '20180406' in filename:
				if '18.' in filename:
					datestr_pred = filename[11:19] + '_2'
				else:
					datestr_pred = filename[11:19] + '_1'
			else:
				datestr_pred = filename[11:19]

			temperatures[datestr_pred] = data_pred['temps'][0]
			total_distance[datestr_pred] = np.sum(np.array(data_pred['dists']))
			total_time[datestr_pred] = data_pred['times'][-1]
			ini_alts[datestr_pred] = data_pred['alts'][0]/1000.
			ini_lats[datestr_pred] = data_pred['lats'][0]
			ini_lons[datestr_pred] = data_pred['lons'][0]

			datestr_gps = pyb_io.match_pred2gps(datestr_pred)
			gps_file = datestr_gps + '.csv'

			with open(dir_gps + gps_file) as csvfile:

				data_gps = np.array(list(csv.reader(csvfile)))
				j = -1

				if 'RB' in data_gps[j][0]:
					d[datestr_pred] = haversine((end_lat_pred, end_lon_pred), (float(data_gps[j][3]), float(data_gps[j][4])))
				else:
					d[datestr_pred] = haversine((end_lat_pred, end_lon_pred), (float(data_gps[j][2]), float(data_gps[j][3])))

	labels = list(d.keys())
	flight_no = np.arange(1, len(labels) + 1, 1)
	vals = list(d.values())
	temps = list(temperatures.values())
	dists = list(total_distance.values())
	times = list(total_time.values())
	alts = list(ini_alts.values())
	lats = list(ini_lats.values())
	lons = list(ini_lons.values())

	vals_sorted = [x for _,x in sorted(zip(labels,vals))]
	mean = np.mean(vals_sorted)
	labels.sort()

	fig = plt.figure()

	purple_dot = mlines.Line2D([], [], color='m', marker='o', linestyle='None', markersize=6, label='Morocco')
	green_dot = mlines.Line2D([], [], color='g', marker='o', linestyle='None', markersize=6, label='Greenland')
	red_dot = mlines.Line2D([], [], color='r', marker='o', linestyle='None', markersize=6, label='Switzerland')
	black_line = mlines.Line2D([], [], color='k', marker='None', linestyle='--', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	red_line = mlines.Line2D([], [], color='r', marker='None', linestyle='--', linewidth=1, label='5 km')

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	for i in range(len(flight_no)):
		if labels[i] in ['20180811', '20180803', '20180806', '20180813']:
			plt.plot(flight_no[i], vals_sorted[i], 'go')
		elif labels[i] == '20181028':
			plt.plot(flight_no[i], vals_sorted[i], 'mo')
		else:
			plt.plot(flight_no[i], vals_sorted[i], 'ro')

	plt.xticks(flight_no, labels, rotation=30)
	plt.xticks(fontsize=8)
	plt.xlabel('Flight date', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', handles=[purple_dot, green_dot, red_dot, black_line, red_line])
	plt.grid(True)
	plt.tight_layout()

	if descent_only:
		fig.savefig(fig_dir + 'dvsflight' + interp + '_startpoint' + next_point + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)
	else:
		fig.savefig(fig_dir + 'dvsflight' + interp + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)

	plt.clf()

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	plt.axvline(np.mean(temps), linestyle='--', color='blue', linewidth=1, label='Average T: ' + str(round(np.mean(temps), 1)) + ' K')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	labels = list(d.keys())
	for i in range(len(labels)):
		if labels[i] in ['20180811', '20180803', '20180806', '20180813']:
			plt.plot(temps[i], vals[i], 'go')
		elif labels[i] == '20181028':
			plt.plot(temps[i], vals[i], 'mo')
		else:
			plt.plot(temps[i], vals[i], 'ro')

	blue_line1 = mlines.Line2D([], [], color='b', marker='None', linestyle='--', linewidth=1, label='Average T: ' + str(round(np.mean(temps), 1)) + ' K')

	plt.xlabel('Initial Temperature [K]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', handles=[purple_dot, green_dot, red_dot, black_line, red_line, blue_line1])
	plt.grid(True)

	if descent_only:
		fig.savefig(fig_dir + 'dvsT' + interp + '_startpoint' + next_point + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)
	else:
		fig.savefig(fig_dir + 'dvsT' + interp + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)

	plt.clf()

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')
	plt.axvline(np.mean(alts), linestyle='--', color='blue', linewidth=1, label='Average altitude: ' + str(round(np.mean(alts), 1)) + ' km')

	labels = list(d.keys())
	for i in range(len(labels)):
		if labels[i] in ['20180811', '20180803', '20180806', '20180813']:
			plt.plot(alts[i], vals[i], 'go')
		elif labels[i] == '20181028':
			plt.plot(alts[i], vals[i], 'mo')
		else:
			plt.plot(alts[i], vals[i], 'ro')

	blue_line2 = mlines.Line2D([], [], color='b', marker='None', linestyle='--', linewidth=1, label='Average altitude: ' + str(round(np.mean(alts), 1)) + ' km')

	# labels = list(d.keys())
	# for i in range(len(alts)):
	# 	xy = (alts[i], vals[i][0] - 1)
	# 	plt.annotate(labels[i], xy = xy, size=5)

	plt.xlabel('Initial Altitude [km]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', handles=[purple_dot, green_dot, red_dot, black_line, red_line, blue_line2])
	plt.grid(True)
	# plt.tight_layout()

	if descent_only:
		fig.savefig(fig_dir + 'dvsh' + interp + '_startpoint' + next_point + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)
	else:
		fig.savefig(fig_dir + 'dvsh' + interp + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)

	plt.clf()

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	plt.axvline(np.mean(lats), linestyle='--', color='blue', linewidth=1, label='Average latitude: ' + str(round(np.mean(lats), 1)) + '$^\circ$')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	labels = list(d.keys())
	for i in range(len(labels)):
		if labels[i] in ['20180811', '20180803', '20180806', '20180813']:
			plt.plot(lats[i], vals[i], 'go')
		elif labels[i] == '20181028':
			plt.plot(lats[i], vals[i], 'mo')
		else:
			plt.plot(lats[i], vals[i], 'ro')

	blue_line3 = mlines.Line2D([], [], color='b', marker='None', linestyle='--', linewidth=1, label='Average latitude: ' + str(round(np.mean(lats), 1)) + '$^\circ$')

	plt.xlabel('Initial Latitude [$^\circ$]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', handles=[purple_dot, green_dot, red_dot, black_line, blue_line3])
	plt.grid(True)

	if descent_only:
		fig.savefig(fig_dir + 'dvslat' + interp + '_startpoint' + next_point + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)
	else:
		fig.savefig(fig_dir + 'dvslat' + interp + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)

	plt.clf()

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	plt.axvline(np.mean(lons), linestyle='--', color='blue', linewidth=1, label='Average longitude: ' + str(round(np.mean(lons), 1)) + '$^\circ$')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	labels = list(d.keys())
	for i in range(len(labels)):
		if labels[i] in ['20180811', '20180803', '20180806', '20180813']:
			plt.plot(lons[i], vals[i], 'go')
		elif labels[i] == '20181028':
			plt.plot(lons[i], vals[i], 'mo')
		else:
			plt.plot(lons[i], vals[i], 'ro')

	blue_line4 = mlines.Line2D([], [], color='b', marker='None', linestyle='--', linewidth=1, label='Average longitude: ' + str(round(np.mean(lons), 1)) + '$^\circ$')

	plt.xlabel('Initial Longitude [$^\circ$]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', handles=[purple_dot, green_dot, red_dot, black_line, red_line, blue_line4])
	plt.grid(True)

	if descent_only:
		fig.savefig(fig_dir + 'dvslon' + interp + '_startpoint' + next_point + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)
	else:
		fig.savefig(fig_dir + 'dvslon' + interp + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)

	plt.clf()

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	plt.axvline(np.mean(times), linestyle='--', color='blue', linewidth=1, label='Average time: ' + str(round(np.mean(times), 1)) + ' min')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	labels = list(d.keys())
	for i in range(len(labels)):
		if labels[i] in ['20180811', '20180803', '20180806', '20180813']:
			plt.plot(times[i], vals[i], 'go')
		elif labels[i] == '20181028':
			plt.plot(times[i], vals[i], 'mo')
		else:
			plt.plot(times[i], vals[i], 'ro')

	blue_line5 = mlines.Line2D([], [], color='b', marker='None', linestyle='--', linewidth=1, label='Average time: ' + str(round(np.mean(times), 1)) + ' min')

	plt.xlabel('Total Flight Time [min]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', handles=[purple_dot, green_dot, red_dot, black_line, red_line, blue_line5])
	plt.grid(True)

	if descent_only:
		fig.savefig(fig_dir + 'dvstime' + interp + '_startpoint' + next_point + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)
	else:
		fig.savefig(fig_dir + 'dvstime' + interp + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)

	plt.clf()

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	plt.axvline(np.mean(dists), linestyle='--', color='blue', linewidth=1, label='Average distance: ' + str(round(np.mean(dists), 1)) + ' km')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	labels = list(d.keys())
	for i in range(len(labels)):
		if labels[i] in ['20180811', '20180803', '20180806', '20180813']:
			plt.plot(dists[i], vals[i], 'go')
		elif labels[i] == '20181028':
			plt.plot(dists[i], vals[i], 'mo')
		else:
			plt.plot(dists[i], vals[i], 'ro')

	blue_line6 = mlines.Line2D([], [], color='b', marker='None', linestyle='--', linewidth=1, label='Average dist. travelled: ' + str(round(np.mean(dists), 1)) + ' km')

	plt.xlabel('Total Distance Travelled [km]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', handles=[purple_dot, green_dot, red_dot, black_line, red_line, blue_line6], prop={'size': 10})
	plt.grid(True)

	if descent_only:
		fig.savefig(fig_dir + 'dvsdist' + interp + '_startpoint' + next_point + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)
	else:
		fig.savefig(fig_dir + 'dvsdist' + interp + '_+' + str(int(drift_time)).zfill(4) + 'min.png', dpi=1000)

	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':

	plot_results()
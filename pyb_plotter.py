""" method to plot vertical speeds and final results """

from scipy import interpolate
from haversine import haversine
from astropy.io import ascii
import datetime as dt
import numpy as np
import time
import csv
import sys
import os

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib as mpl

import param_file as p
import pyb_io

def get_rates(params=None, run=None):

	FMT = "%H:%M:%S"

	descent_rates_pred = {}
	descent_rates_gps = {}
	gps_indices = {}
	alts_pred = {}
	alts_gps = {}

	dir_pred = '/home/ellen/Desktop/SuperBIT/Output/'

	if params == None:
		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		else:
			next_point = '0'
		interpolate = p.interpolate
		drift_time = p.drift_time

		params = [descent_only, next_point, interpolate, drift_time]

	else:

		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		else:
			next_point = '0'
		interpolate = bool(params[2])
		drift_time = float(params[3])

	if run == None:

		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)
		files = [filename for filename in os.listdir(dir_pred) if now_str in filename]
		folder = now_str + '_' + str(len(files)-1) + '/'

	else:
		folder = run + '/'

	dir_pred += folder + 'trajectories/'
	dir_gps = '/home/ellen/Desktop/SuperBIT/Flight_data/'

	add = int(next_point)

	for filename in os.listdir(dir_pred):

		data = ascii.read(dir_pred + filename)

		alts = np.array(data['alts'])
		times = np.array(data['times'])
		descent_speeds = np.array(data['speeds'])

		if '20180406' in filename:
			if '8.9' in filename:
				datestr_pred = filename[11:19] + '_1'
			else:
				datestr_pred = filename[11:19] + '_2'
		else:
			datestr_pred = filename[11:19]

		descent_rates_pred[datestr_pred] = descent_speeds
		alts_pred[datestr_pred] = np.array(alts)

		datestr_gps = pyb_io.match_pred2gps(datestr_pred)
		gps_file = datestr_gps + '.csv'

		with open(dir_gps + gps_file) as csvfile:

			data_gps = np.array(list(csv.reader(csvfile)))

			alts, times = [], []

			for i in range(len(data_gps)):
				if 'RB' in data_gps[i][0]:
					alts.append(int(data_gps[i][5]))
					times.append(data_gps[i][2])
				else:
					alts.append(int(data_gps[i][4]))
					times.append(data_gps[i][1])

			ind = int(np.where(alts == np.max(alts))[0]) + add
			alt0 = np.max(alts)

			descent_rates_gps[datestr_pred] = np.array([float((alts[i+1] - alts[i]))/(dt.datetime.strptime(times[i+1], FMT) - dt.datetime.strptime(times[i], FMT)).seconds for i in range(ind, len(alts)-2)])
			alts_gps[datestr_pred] = np.array([alts[i] for i in range(ind, len(alts) -2)])

	return (descent_rates_gps, descent_rates_pred, alts_gps, alts_pred), dir_pred

def plot_rates(data=None, params=None, dir_pred=None, run=None):

	time0 = time.time()

	print('Plotting descent speeds...')

	if data == None or dir_pred == None:
		(descent_rates_gps, descent_rates_pred, alts_gps, alts_pred), dir_pred = get_rates(params=params, run=run)

	fig_dir = dir_pred + '../figs/DescentRates/'

	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	pred_keys = list(descent_rates_pred.keys())
	gps_keys = list(descent_rates_gps.keys())

	descent_rates_gps_vals = list(descent_rates_gps.values())
	descent_rates_pred_vals = list(descent_rates_pred.values())

	alts_gps_vals = list(alts_gps.values())
	alts_pred_vals = list(alts_pred.values())

	fig = plt.figure()
	plt.rc('axes', axisbelow=True)

	for i in range(len(descent_rates_gps_vals)):

		### x-coordinates at which to evaluate the interpolated values, x&y coordinates used in the interpolation
		descent_rates_interp_pred = np.interp(alts_gps_vals[i][::-1], alts_pred_vals[i][::-1], descent_rates_pred_vals[i][::-1])[::-1]

		plt.plot(descent_rates_gps_vals[i], descent_rates_interp_pred, 'ro--', markersize=2, label='Results')
		plt.plot(descent_rates_gps_vals[i], descent_rates_gps_vals[i], 'b--', linewidth=1, label=r'$v_{pred}=v_{true}$')

		plt.ylabel(r'$v_{pred}$ [m/s]', fontsize=15)
		plt.xlabel(r'$v_{true}$ [m/s]', fontsize=15)
		plt.legend(loc='best')
		plt.grid(True)
		fig.savefig(fig_dir + 'vdescent_' + pred_keys[i] + '.png')

		plt.clf()

		plt.plot(alts_pred_vals[i], descent_rates_pred_vals[i], 'bo--', markersize=2, label='Predicted')
		plt.plot(alts_gps_vals[i], descent_rates_gps_vals[i], 'ro--', markersize=2, label='True')

		plt.ylabel('Descent rate [m/s]', fontsize=15)
		plt.xlabel('Altitude [m]', fontsize=15)
		plt.grid(True)
		plt.legend(loc='best')
		fig.savefig(fig_dir + 'vdescent_' + pred_keys[i] + '.png', dpi=1000)

		plt.clf()

	for i in range(len(descent_rates_gps_vals)):

		descent_rates_interp_pred = np.interp(alts_gps_vals[i][::-1], alts_pred_vals[i][::-1], descent_rates_pred_vals[i][::-1])[::-1]
		plt.plot(descent_rates_gps_vals[i], descent_rates_interp_pred, 'ro-', markersize=1, linewidth=0.5)

	minimum = min([min(descent_rates_gps_vals[i]) for i in range(len(descent_rates_gps_vals))])
	testx = np.arange(minimum, 0, 0.5)

	plt.plot(testx, testx, 'b--', linewidth=1, label=r'$v_{pred}=v_{true}$')
	plt.ylabel(r'$v_{pred}$ [m/s]', fontsize=15)
	plt.xlabel(r'$v_{true}$ [m/s]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best', prop={'size':10})

	fig.savefig(fig_dir + 'vdescent_all.png', dpi=1000)

	print('Total time elapsed: %.1f s' % (time.time() - time0))

def plot_results(run=None, params=None):

	time0 = time.time()

	print('Plotting final results...')

	if params == None:
		descent_only = p.descent_only
		if descent_only:
			next_point = p.next_point
		else:
			next_point = '0'
		interpolate = p.interpolate
		drift_time = p.drift_time

		params = [descent_only, next_point, interpolate, drift_time]

	else:
		descent_only = bool(params[0])
		if descent_only:
			next_point = str(params[1])
		else:
			next_point = '0'
		interpolate = bool(params[2])
		drift_time = float(params[3])

	dir_base = '/home/ellen/Desktop/SuperBIT/'
	dir_pred =  dir_base + 'Output/'
	dir_gps =  dir_base + 'Flight_data/'

	if run == None:
		now = dt.datetime.now()
		now_str = str(now.year) + str(now.month).zfill(2) + str(now.day).zfill(2)
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

			if 'RB' in data_gps[-1][0]:
				d[datestr_pred] = haversine((end_lat_pred, end_lon_pred), (float(data_gps[-1][3]), float(data_gps[-1][4])))
			else:
				d[datestr_pred] = haversine((end_lat_pred, end_lon_pred), (float(data_gps[-1][2]), float(data_gps[-1][3])))

	#### get values ####

	labels = list(d.keys())
	flight_no = np.arange(1, len(labels) + 1, 1)

	ds = list(d.values())
	temps = list(temperatures.values())
	dists = list(total_distance.values())
	times = list(total_time.values())
	alts = list(ini_alts.values())
	lats = list(ini_lats.values())
	lons = list(ini_lons.values())

	ds_sorted = [x for _,x in sorted(zip(labels,ds))]
	mean = np.mean(ds_sorted)
	labels.sort()

	### make plots ###

	fig = plt.figure()

	#### error in distance vs. date ####

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Avg. error: ' + str(round(mean, 1)) + ' km')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	plt.plot(flight_no, ds_sorted, 'ro')

	plt.xticks(flight_no, labels, rotation=30)
	plt.xticks(fontsize=8)
	plt.xlabel('Flight date', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', prop={'size': 10})
	plt.grid(True)
	plt.tight_layout()

	fig.savefig(fig_dir + 'dvsflight.png', dpi=1000)

	#### error in distance vs. initial temperature ####

	plt.clf()

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Avg. error: ' + str(round(mean, 1)) + ' km')
	plt.axvline(np.mean(temps), linestyle='--', color='blue', linewidth=0.5, label='Avg. T: ' + str(round(np.mean(temps), 1)) + ' K')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	plt.plot(temps, ds, 'ro')

	plt.xlabel('Initial Temperature [K]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', prop={'size': 10})
	plt.grid(True)

	fig.savefig(fig_dir + 'dvsT.png', dpi=1000)

	#### error in distance vs. total flight time ####

	plt.clf()

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	plt.axvline(np.mean(times), linestyle='--', color='blue', linewidth=0.5, label='Avg. time: ' + str(round(np.mean(times), 1)) + ' min')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	plt.plot(times, ds, 'ro')

	plt.xlabel('Total Flight Time [min]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', prop={'size': 10})
	plt.grid(True)

	fig.savefig(fig_dir + 'dvstime.png', dpi=1000)

	#### error in distance vs. total flight distance ####

	plt.clf()

	plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
	plt.axvline(np.mean(dists), linestyle='--', color='blue', linewidth=0.5, label='Avg. distance: ' + str(round(np.mean(dists), 1)) + ' km')
	plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	plt.plot(dists, ds, 'ro')

	plt.xlabel('Total Distance Travelled [km]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.legend(loc='best', prop={'size': 10})
	plt.grid(True)

	fig.savefig(fig_dir + 'dvsdist.png', dpi=1000)

	#### error in distance vs. initial altitude ####

	# plt.clf()

	# plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Avg. error: ' + str(round(mean, 1)) + ' km')
	# plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')
	# plt.axvline(np.mean(alts), linestyle='--', color='blue', linewidth=0.5, label='Avg. altitude: ' + str(round(np.mean(alts), 1)) + ' km')

	# plt.plot(alts, ds, 'ro')

	# plt.xlabel('Initial Altitude [km]', fontsize=15)
	# plt.ylabel('Error in distance [km]', fontsize=15)
	# plt.legend(loc='best', prop={'size': 10})
	# plt.grid(True)

	# fig.savefig(fig_dir + 'dvsh.png', dpi=1000)

	#### error in distance vs. latitude ####

	# plt.clf()

	# plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Avg. error: ' + str(round(mean, 1)) + ' km')
	# plt.axvline(np.mean(lats), linestyle='--', color='blue', linewidth=0.5, label='Avg. latitude: ' + str(round(np.mean(lats), 1)) + '$^\circ$')
	# plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	# plt.plot(lats, ds, 'ro')

	# plt.xlabel('Initial Latitude [$^\circ$]', fontsize=15)
	# plt.ylabel('Error in distance [km]', fontsize=15)
	# plt.legend(loc='best', prop={'size': 10})
	# plt.grid(True)

	# fig.savefig(fig_dir + 'dvslat.png', dpi=1000)

	#### error in distance vs. longitude ####

	# plt.clf()

	# plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Avg. error: ' + str(round(mean, 1)) + ' km')
	# plt.axvline(np.mean(lons), linestyle='--', color='blue', linewidth=0.5, label='Avg. longitude: ' + str(round(np.mean(lons), 1)) + '$^\circ$')
	# plt.axhline(5, linestyle='--', color='red', linewidth=1, label='5 km')

	# plt.plot(lons, ds, 'ro')

	# plt.xlabel('Initial Longitude [$^\circ$]', fontsize=15)
	# plt.ylabel('Error in distance [km]', fontsize=15)
	# plt.legend(loc='best', prop={'size': 10})
	# plt.grid(True)

	# fig.savefig(fig_dir + 'dvslon.png', dpi=1000)

	print('Total time elapsed: %.1f s' % (time.time() - time0))

if __name__ == '__main__':
	plot_rates()
	plot_results()
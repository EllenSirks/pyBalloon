from scipy.optimize import curve_fit
from collections import defaultdict
from haversine import haversine
from astropy.io import ascii
import numpy as np
import sys, os
import time
import csv

import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.rcParams['axes.axisbelow'] = True
plt.rcParams["font.family"] = "serif"

import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'dejavuserif'

import pyb_info_searcher
import pyb_plotter
import pyb_aux

import param_file as p

def func(x, m, c):
	return m*x + c

def forecast_flights_plotter(run=None):

	##################################################################################################

	in_dir = p.path + 'Output/'
	dir_gps = p.path + 'Flight_data/'
	fig_dir = p.path + 'Output/' + run + '/Figs/Results/'

	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	err_dists = {}
	tot_dists = {}

	dir_pred = in_dir + run + '/Trajectories/'

	fname = 'params.txt'
	lines = [line.rstrip('\n').split(' ') for line in open(in_dir + run + '/' + fname)]

	for i in range(len(lines)):
		if lines[i][0] == 'difference':
			hr_diffs = lines[i][6:-1]
			hr_diffs[-1] = hr_diffs[-1][:-1]
			hr_diffs = np.array([float(j) for j in hr_diffs])
			break

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

		if filename[-5] == ')':
			no = hr_diffs[0]
		else:
			no = hr_diffs[int(filename[-5])-1]

		err_dists.setdefault(datestr_pred, {})
		tot_dists.setdefault(datestr_pred, {})

		tot_dist = np.sum(np.array(data_pred['dists']))

		datestr_gps = pyb_plotter.match_pred2gps(datestr_pred)
		gps_file = datestr_gps + '.csv'

		with open(dir_gps + gps_file) as csvfile:

			data_gps = np.array(list(csv.reader(csvfile)))

			if 'RB' in data_gps[-1][0]:
				err_dist = haversine((end_lat_pred, end_lon_pred), (float(data_gps[-1][3]), float(data_gps[-1][4])))
			else:
				err_dist = haversine((end_lat_pred, end_lon_pred), (float(data_gps[-1][2]), float(data_gps[-1][3])))

		err_dists[datestr_pred][no] = err_dist
		tot_dists[datestr_pred][no] = tot_dist

	datestrs = list(err_dists.keys())
	nos = list(err_dists[datestrs[0]].keys())

	err_dists_2 = {}
	tot_dists_2 = {}

	for no in nos:
		err_dists_2[no] = []
		tot_dists_2[no] = []
		for datestr in datestrs:
			err_dists_2[no].append(err_dists[datestr][no])
			tot_dists_2[no].append(tot_dists[datestr][no])

	colors = ['deeppink', 'firebrick', 'darkorange', 'cornflowerblue', 'limegreen'][::-1]
	# markers = ['o', '^', 's', 'P', '*', 'D', 'v', '<', '>', 'X', 'd', r'$o$', 'p', '.', 'h', r'$v$', r'$u$', r'$n$']
	markers = ['$' + str(i) + '$' for i in range(len(datestrs))]#

	fig = plt.figure()

	plt.xlabel('Total Flight Time [min]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)

	for c, hr_diff in zip(colors, hr_diffs):

		popt, pcov = curve_fit(func, tot_dists_2[hr_diff], err_dists_2[hr_diff])
		model_y = func(np.array(tot_dists_2[hr_diff]), popt[0], popt[1])

		plt.axhline(np.mean(err_dists_2[hr_diff]), linewidth=1, color=c, linestyle='--')
		plt.plot(tot_dists_2[hr_diff], model_y, '-', color=c, linewidth=1)

	for c, hr_diff in zip(colors, hr_diffs):
		for i in range(len(datestrs)):
			if i == 0:
				plt.plot(tot_dists_2[hr_diff][i], err_dists_2[hr_diff][i], 'o', color=c, markersize=6, marker=markers[i], label=str(hr_diff) + ' hrs, avg. = ' + str(round(np.mean(err_dists_2[hr_diff]), 2)) + ' km')
			else:
				plt.plot(tot_dists_2[hr_diff][i], err_dists_2[hr_diff][i], 'o', color=c, markersize=6, marker=markers[i])

	plt.legend(loc='best')

	plt.grid(True)
	plt.tight_layout()


	fig.savefig(fig_dir + 'forecast_results.png')

	plt.clf()

	plt.xlabel('Total Flight Time [min]', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)

	for c, hr_diff in zip(colors, hr_diffs):

		popt, pcov = curve_fit(func, tot_dists_2[hr_diff], err_dists_2[hr_diff])
		model_y = func(np.array(tot_dists_2[hr_diff]), popt[0], popt[1])

		plt.axhline(np.mean(err_dists_2[hr_diff]), linewidth=1, color=c, linestyle='--', label=str(hr_diff) + ' hrs, avg. = ' + str(round(np.mean(err_dists_2[hr_diff]), 2)) + ' km')
		plt.plot(tot_dists_2[hr_diff], model_y, '-', color=c, linewidth=1)
		plt.plot(tot_dists_2[hr_diff], err_dists_2[hr_diff], 'o', color=c, linewidth=1, markersize=3)

	plt.legend(loc='best')

	plt.grid(True)
	plt.tight_layout()

	fig_dir = p.path + 'Output/' + run + '/Figs/Results/'

	fig.savefig(fig_dir + 'forecast_results_points.png')

#################################################################################################

def forecast_test_plotter(run=None, loc0=None):

	lats = {}
	lons = {}

	in_dir = p.path + 'Output/'
	dir_pred = in_dir + run + '/Trajectories/'
	fig_dir = p.path + 'Output/' + run + '/Figs/Results/'

	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	lat0, lon0 = loc0

	if lon0 > 180:
		lon0 -= 360

	fname = 'params.txt'
	lines = [line.rstrip('\n').split(' ') for line in open(in_dir + run + '/' + fname)]

	for i in range(len(lines)):
		if lines[i][0] == 'difference':
			hr_diffs = lines[i][6:-1]
			hr_diffs[-1] = hr_diffs[-1][:-1]
			hr_diffs = np.array([float(j) for j in hr_diffs])
			break

	for filename in os.listdir(dir_pred):

		if filename[-5] == ')':
			nr = 1
		else:
			nr = int(filename[-5])

		data_pred = ascii.read(dir_pred + filename)
		datestr_pred = filename[11:19]

		lats.setdefault(datestr_pred, {})
		lons.setdefault(datestr_pred, {})

		lats[datestr_pred][nr] = data_pred['lats'][-1]
		lons[datestr_pred][nr] = data_pred['lons'][-1] - 360

	datestrs = list(lats.keys())
	datestrs.sort()

	colors = cm.rainbow(np.linspace(0, 1, len(datestrs)))
	markers = ['o', '^', 's', '*', 'P']

	fig = plt.figure()
	plt.xlabel('Latitude', fontsize=15)
	plt.ylabel('Longitude', fontsize=15)

	for c, datestr in zip(colors, datestrs):
		plt.plot(np.array(lats[datestr].values()), np.array(lons[datestr].values()), '--', color=c, linewidth=1)
		for i in range(1, len(markers) + 1):
			if i == 0:
				plt.plot(lats[datestr][i], lons[datestr][i], 'o', color=c, markersize=4, marker=markers[i-1], label=str(hr_diffs[i-1]) + ' hrs. diff.')
			else:
				plt.plot(lats[datestr][i], lons[datestr][i], 'o', color=c, markersize=4, marker=markers[i-1])

	plt.plot(lat0, lon0, 'ko', markersize=5, label='Starting Location')

	plt.legend(loc='best')
	plt.grid(True)
	plt.tight_layout()

	fig.savefig(fig_dir + 'forecast_results_2018Sept.png')

	plt.clf()

	plt.xlabel('Latitude', fontsize=15)
	plt.ylabel('Longitude', fontsize=15)

	for datestr, c in zip(datestrs, colors):
		if datestr == datestrs[0]:
			plt.plot(np.array(lats[datestr].values())[0], np.array(lons[datestr].values())[0], 'o', color=c, markersize=4, label=str(hr_diffs[0]) + ' hrs. diff.')
			plt.plot(np.array(lats[datestr].values())[-1], np.array(lons[datestr].values())[-1], 'o', color=c, markersize=4, marker='s', label=str(hr_diffs[-1]) + ' hrs. diff.')
		else:
			plt.plot(np.array(lats[datestr].values())[0], np.array(lons[datestr].values())[0], 'o', color=c, markersize=4,)
			plt.plot(np.array(lats[datestr].values())[-1], np.array(lons[datestr].values())[-1], 'o', color=c, markersize=4, marker='s')

		plt.plot(np.array(lats[datestr].values()), np.array(lons[datestr].values()), '--', color=c, linewidth=1)

	plt.plot(lat0, lon0, 'ko', markersize=5, label='Starting Location')

	plt.legend(loc='lower right')
	plt.grid(True)
	plt.tight_layout()

	fig.savefig(fig_dir + 'forecast_results_2018Sept_2.png')

	# ##################################################################################################

	delta_xs = {}

	for datestr in datestrs:

		loc0 = (lats[datestr][1], lons[datestr][1])
		delta_x = {}

		for i in range(1, len(lats[datestr])+1):
			loc = (lats[datestr][i], lons[datestr][i])
			d = haversine(loc0, loc)
			delta_x[i] = d

		delta_xs[datestr] = delta_x

	plt.clf()

	plt.xlabel('$\Delta t_{forecast}$ [hours]', fontsize=15)
	plt.ylabel('$\Delta x$ [km]', fontsize=15)

	step = hr_diffs[1]-hr_diffs[0]
	plt.xticks(np.arange(min(hr_diffs), max(hr_diffs) + step, step=step), hr_diffs)

	for datestr, c in zip(datestrs, colors):
		plt.plot(np.array(hr_diffs), np.array(delta_xs[datestr].values()), 'o', color=c, markersize=3)

	plt.grid(True)
	plt.tight_layout()

	fig.savefig(fig_dir + 'forecast_results_2018Sept_xvst.png')

	plt.clf()

	plt.xlabel('$\Delta t_{forecast}$ [hours]', fontsize=15)
	plt.ylabel('$\Delta x$ [km]', fontsize=15)

	step = hr_diffs[1]-hr_diffs[0]
	plt.xticks(np.arange(min(hr_diffs), max(hr_diffs) + step, step=step), hr_diffs)

	for datestr, c in zip(datestrs, colors):
		plt.plot(np.array(hr_diffs), np.array(delta_xs[datestr].values()), 'o--', color=c, linewidth=1, markersize=2)

	line1 = matplotlib.lines.Line2D([0], [0], color=colors[0], linewidth=1, linestyle='--')
	line2 = matplotlib.lines.Line2D([0], [0], color=colors[-1], linewidth=1, linestyle='--')

	min_date = datestrs[0][6:] + '-' + datestrs[0][4:6] + '-' + datestrs[0][:4]
	max_date = datestrs[-1][6:] + '-' + datestrs[-1][4:6] + '-' + datestrs[-1][:4]

	plt.legend([line1, line2], [min_date, max_date])
	plt.grid(True)
	plt.tight_layout()

	fig.savefig(fig_dir + 'forecast_results_2018Sept_xvst_lines.png')

	maxs = []
	mins = []

	for i in range(len(hr_diffs)):

		mins.append(min([delta_xs[datestr][i+1] for datestr in datestrs]))
		maxs.append(max([delta_xs[datestr][i+1] for datestr in datestrs]))

	plt.clf()

	plt.xlabel('$\Delta t_{forecast}$ [hours]', fontsize=15)
	plt.ylabel('$\Delta x$ [km]', fontsize=15)

	step = hr_diffs[1]-hr_diffs[0]
	plt.xticks(np.arange(min(hr_diffs), max(hr_diffs) + step, step=step), hr_diffs)

	plt.plot(np.array(hr_diffs), np.array(mins), '-', color='r', linewidth=1)
	plt.plot(np.array(hr_diffs), np.array(maxs), '-', color='r', linewidth=1)

	plt.fill_between(np.array(hr_diffs), mins, maxs, color='red', alpha='0.2')

	plt.grid(True)
	plt.tight_layout()

	fig.savefig(fig_dir + 'forecast_results_2018Sept_xvst_minmax.png')

##################################################################################################

if __name__ == '__main__':

	time0 = time.time()

	run = sys.argv[1]

	# forecast_flights_plotter(run=run)

	loc0 = float(sys.argv[2]), float(sys.argv[3])
	forecast_test_plotter(run=run, loc0=loc0)

	print('Total time elapsed: %.1f s' % (time.time() - time0))


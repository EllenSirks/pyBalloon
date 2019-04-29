from haversine import haversine
import matplotlib.pyplot as plt
from astropy.io import ascii
import match_files as mf
import endpoints as ep
import numpy as np
import csv
import sys
import os

def calc_results(trunc='True', next_point='0', re_endpoints='no'):

	### Which point after the max. altitude do you want to check ###

	if next_point == '2':
		fext = 'next2_points/'
	elif next_point == '1':
		fext = 'next_points/'
	else:
		fext = 'all_points/'

	dir = '/home/ellen/Desktop/SuperBIT/'

	date_keys = ['20180406', '20180304', '20181512', '20190101', '20180813', '20180803', '20180506', '20180630', '20180616', '20180511', '20181028', '20180519', '20180806', '20180811', '20190202'. '20190209', '20190210']
	# date_keys = ['20180406', '20181512', '20190101', '20180813', '20180506', '20180630', '20180511', '20180519', '20180806', '20180811']

	dir_pred =  dir + 'Weather_data/Endpoints/' + fext
	dir_gps =  dir + 'Flight_data/'

	d = {}
	deltat = {}

	data_times = ascii.read(dir_gps + 'burst_times.txt')

	dates = data_times['date']
	burst_times = data_times['burst_time']

	for filename in os.listdir(dir_pred):

		if filename.endswith('dat') and filename[:8] in date_keys:

			### calculate difference in landing point ###

			if re_endpoints == 'yes' and '201808' not in filename:
				end_lat_pred, end_lon_pred = ep.get_endpoint(filename=filename[:-13] + '_trajectory.dat', fext=fext)

			elif re_endpoints == 'no' or '201808' in filename:
				data_pred = ascii.read(dir_pred + filename)

				end_lat_pred = data_pred['lat'][-1]
				end_lon_pred = data_pred['lon'][-1]

			weather_time = int(float(filename[9:13]) / 100) + int(filename[14:17])

			if '20180406' in filename:
				if '0600' in filename:
					datestr_pred = filename[:8] + '_1'
				else:
					datestr_pred = filename[:8] + '_2'
			else:
				datestr_pred = filename[:8]

			datestr_gps = mf.match_pred2gps(datestr_pred)
			gps_file = datestr_gps + '.csv'

			with open(dir_gps + gps_file, newline='') as csvfile:

				data_gps = np.array(list(csv.reader(csvfile)))

				if trunc == 'False':

					ext ='_ep_' + re_endpoints 
					j = -1

				### gps tracks sometimes include moving the tracker to the road, need to cut this part off ###
				elif trunc == 'True':

					alts = []
					ext = '_trunc_ep_' + re_endpoints 

					for i in range(len(data_gps)):
						if 'RB' in data_gps[i][0]:
							alts.append(int(data_gps[i][5]))
						else:
							alts.append(int(data_gps[i][4]))

					alt0 = np.max(alts)
					ind, = np.where(alts == np.max(alts))[0]

					for j in range(ind + 1, len(alts)):
						if alts[j] < alt0:
							alt0 = alts[j]
						elif alts[j] >= alt0:
							j -= 1
							break

		if 'RB' in data_gps[j][0]:
			d.setdefault(datestr_pred, []).append(haversine((end_lat_pred, end_lon_pred), (float(data_gps[j][3]), float(data_gps[j][4]))))
		else:
			d.setdefault(datestr_pred, []).append(haversine((end_lat_pred, end_lon_pred), (float(data_gps[j][2]), float(data_gps[j][3]))))

		burst_time = float(burst_times[np.where(dates == datestr_pred)[0]])
		deltat.setdefault(datestr_pred, []).append(round(burst_time - weather_time, 2))

	return deltat, d, fext

def plot_results(deltat=None, d=None, trunc='True', fext='all_points', re_endpoints = 'no'):

	if deltat == None or d == None:
		deltat, d, fext = calc_results(trunc=trunc, re_endpoints=re_endpoints)
	else:
		deltat, d, fext = deltat, d, fext

	ext = '_'

	if trunc == 'True':
		ext += 'trunc_'

	ext += 'ep_' + re_endpoints

	fig_dir = '/home/ellen/Desktop/SuperBIT/figs/Results/' + fext

	if not os.path.exists(fig_dir):
	    os.makedirs(fig_dir)

	keys = list(d.keys())
	vals = list(d.values())
	tims = list(deltat.values())

	values = [item for sublist in list(d.values()) for item in sublist]
	times = [item for sublist in list(deltat.values()) for item in sublist]

	min_vals = [vals[i][np.where(np.array(vals[i]) == min(vals[i]))[0][0]] for i in range(len(vals))]
	min_tims = [tims[i][np.where(np.array(vals[i]) == min(vals[i]))[0][0]] for i in range(len(vals))]

	fig = plt.figure()

	plt.axhline(np.mean(values), linestyle='--', linewidth=1, label='Average: ' + str(round(np.mean(values), 3)))

	plt.plot(times, values, 'ro')
	plt.xlabel('Delta t [hr]', fontsize=15)
	plt.ylabel('Distance [km]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best')

	fig.savefig(fig_dir + 'd_vs_t' + ext + '.png')

	plt.clf()

	plt.axhline(np.mean(min_vals), linestyle='--', linewidth=1, label='Average: ' + str(round(np.mean(min_vals), 3)))

	plt.plot(min_tims, min_vals, 'ro')
	plt.xlabel('Delta t [hr]', fontsize=15)
	plt.ylabel('Distance [km]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best')

	fig.savefig(fig_dir + 'mind_vs_t' + ext + '.png')

	plt.clf()

	plt.axhline(np.mean(values), linestyle='--', linewidth=1, label='Average: ' + str(round(np.mean(values), 3)))

	plt.plot(np.abs(times), values, 'ro')
	plt.xlabel('Delta t [hr]', fontsize=15)
	plt.ylabel('Distance [km]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best')

	fig.savefig(fig_dir + 'd_vs_tabs' + ext + '.png')

	plt.clf()

	plt.axhline(np.mean(min_vals), linestyle='--', linewidth=1, label='Average: ' + str(round(np.mean(min_vals), 3)))

	plt.plot(np.abs(min_tims), min_vals, 'ro')
	plt.xlabel('Delta t [hr]', fontsize=15)
	plt.ylabel('Distance [km]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best')

	fig.savefig(fig_dir + 'mind_vs_tabs' + ext + '.png')

	plt.clf()

	plt.axhline(np.mean([v for val in vals for v in val]), linestyle='--', linewidth=1, label='Average: ' + str(round(np.mean([v for val in vals for v in val]), 3)))

	for i in range(len(vals)):

		if tims[i][0] >= tims[i][1]:
			xy = (tims[i][0], vals[i][0])
			xytext = (tims[i][0] + 0.1, vals[i][0])
		else:
			xy = (tims[i][1], vals[i][1])
			xytext = (tims[i][1] - 0.1, vals[i][1])

		plt.plot(tims[i] , vals[i], 'ro--', markersize=4, linewidth=0.5)
		plt.annotate(keys[i], xy = xy, xytext=xytext, size=5)

	plt.xlabel('Delta t [hr]', fontsize=15)
	plt.ylabel('Distance [km]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best')

	fig.savefig(fig_dir + 'd_vs_t_con' + ext + '.png')

	plt.clf()

	plt.axhline(np.mean([v for val in vals for v in val]), linestyle='--', linewidth=1, label='Average: ' + str(round(np.mean([v for val in vals for v in val]), 3)))

	for i in range(len(vals)):

		if np.abs(tims[i][0]) >= np.abs(tims[i][1]):
			xy = (np.abs(tims[i][0]), vals[i][0])
			xytext = (tims[i][0] + 0.1, vals[i][0])
		else:
			xy = (np.abs(tims[i][1]), vals[i][1])
			xytext = (np.abs(tims[i][1]) - 0.1, vals[i][1])

		plt.plot(np.abs(tims[i]), vals[i], 'ro--', markersize=4, linewidth=0.5)
		plt.annotate(keys[i], xy = xy, xytext=xytext, size=5)

	plt.xlabel('Delta t [hr]', fontsize=15)
	plt.ylabel('Distance [km]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best')

	fig.savefig(fig_dir + 'd_vs_tabs_con' + ext + '.png')

	plt.clf()

if __name__ == '__main__':

	next_point = input('Start at point after max. altitude: ')
	re_endpoints = input('Check if prediction endpoints go below ground: ')

	deltat, d, fext = calc_results(next_point=next_point ,re_endpoints=re_endpoints)

	plot_results(deltat=deltat, d=d, fext=fext, re_endpoints=re_endpoints)

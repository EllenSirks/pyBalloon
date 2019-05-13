from haversine import haversine
import matplotlib.pyplot as plt
from astropy.io import ascii
import match_files as mf
import endpoints as ep
import numpy as np
import csv
import sys
import os

next_point = '0'

ext = 'start_point' + next_point + '/'

dir_pred = '/home/ellen/Desktop/SuperBIT/Weather_data/Endpoints/' + ext
dir_gps = '/home/ellen/Desktop/SuperBIT/Flight_data/'
fig_dir = '/home/ellen/Desktop/SuperBIT/figs/Results/' + ext + 'extra/'

def calc_results(trunc, re_endpoints='no'):

	d = {}
	deltat = {}

	data_times = ascii.read(dir_gps + 'burst_times.txt')

	dates = data_times['date']
	burst_times = data_times['burst_time']

	for filename in os.listdir(dir_pred):

		if filename.endswith('dat'):

			### calculate difference in landing point ###

			if re_endpoints == 'yes' and '201808' not in filename:
				end_lat_pred, end_lon_pred = ep.get_endpoint(filename[:-13] + '_trajectory.dat')

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

					for j in range(len(alts)):
						if j > ind and alts[j] < alt0:
							alt0 = alts[j]
						elif j > ind and alts[j] > alt0:
							break

		if 'RB' in data_gps[j][0]:
			d.setdefault(datestr_pred, []).append(haversine((end_lat_pred, end_lon_pred), (float(data_gps[j][3]), float(data_gps[j][4]))))
		else:
			d.setdefault(datestr_pred, []).append(haversine((end_lat_pred, end_lon_pred), (float(data_gps[j][2]), float(data_gps[j][3]))))

		burst_time = float(burst_times[np.where(dates == datestr_pred)[0]])
		deltat.setdefault(datestr_pred, []).append(round(burst_time - weather_time, 2))

	labels = list(d.keys())
	flight_no = np.arange(1, len(labels) + 1, 1)
	d = list(d.values())
	deltats = list(deltat.values())

	deltatss = [x for _,x in sorted(zip(labels, deltats))]
	ds = [x for _,x in sorted(zip(labels, d))]

	labels.sort()

	fig = plt.figure()
	plt.rc('axes', axisbelow=True)

	cm = plt.cm.get_cmap('RdYlBu')

	x = [flight_no[key] for key in range(len(labels)) for i in range(len(ds[key]))]
	sc = plt.scatter(x, [item for sublist in ds for item in sublist], c=[item for sublist in deltatss for item in sublist], s=100, cmap=cm)
	clb = plt.colorbar(sc)
	clb.set_label('delta t [hr]', fontsize=15)

	plt.xlabel('Flight date', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)

	plt.xticks(np.arange(1, len(d) + 1, step = 1 ), labels=labels, rotation=60)
	plt.grid(True)
	# plt.legend(loc='best')
	plt.tight_layout()

	fig.savefig(fig_dir + 'results' + ext + '.png')

	plt.close()

	return deltatss, ds, labels

if __name__ == '__main__':

	re_endpoints = input('Check if prediction endpoints go below ground: ')

	deltat1, d1, labels1 = calc_results('True', re_endpoints=re_endpoints)
	deltat2, d2, labels2 = calc_results('False', re_endpoints=re_endpoints)

	flight_no = np.arange(1, len(labels1) + 1, 1)

	cm = plt.cm.get_cmap('RdYlBu')

	fig = plt.figure()
	plt.rc('axes', axisbelow=True)

	x = [flight_no[i] for i in range(len(labels1)) for j in range(len(d1[i])) if np.abs(deltat1[i][j]) == min(np.abs(deltat1[i]))]
	y = [d1[i][j] for i in range(len(d1)) for j in range(len(d1[i])) if np.abs(deltat1[i][j]) == min(np.abs(deltat1[i]))]
	c = [deltat1[i][j] for i in range(len(d1)) for j in range(len(d1[i])) if np.abs(deltat1[i][j]) == min(np.abs(deltat1[i]))]

	sc = plt.scatter(x, y, c=c, s=100, cmap=cm)
	clb = plt.colorbar(sc)
	plt.axhline(np.mean(y), linestyle='--', linewidth=1, label='Average: ' + str(round(np.mean(y), 2)))
	clb.set_label('delta t [hr]', fontsize=15)
	plt.xlabel('Flight date', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.xticks(np.arange(1, len(d1) + 1, step = 1), labels=labels1, rotation=60)
	plt.grid(True)
	plt.legend(loc='best')
	plt.tight_layout()

	fig.savefig(fig_dir + 'results_minT_ep_' + re_endpoints + '.png')
	plt.close()

	fig = plt.figure()
	plt.rc('axes', axisbelow=True)

	x = [flight_no[i] for i in range(len(labels1)) for j in range(len(d1[i])) if d1[i][j] == min(d1[i])]
	y = [d1[i][j] for i in range(len(d1)) for j in range(len(d1[i])) if d1[i][j] == min(d1[i])]
	c = [deltat1[i][j] for i in range(len(d1)) for j in range(len(d1[i])) if d1[i][j] == min(d1[i])]

	sc = plt.scatter(x, y, c=c, s=100, cmap=cm)
	clb = plt.colorbar(sc)
	plt.axhline(np.mean(y), linestyle='--', linewidth=1, label='Average: ' + str(round(np.mean(y), 2)))
	clb.set_label('delta t [hr]', fontsize=15)
	plt.xlabel('Flight date', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.xticks(np.arange(1, len(d1) + 1, step = 1), labels=labels1, rotation=60)
	plt.grid(True)
	plt.legend(loc='best')
	plt.tight_layout()

	fig.savefig(fig_dir + 'results_minD_ep_' + re_endpoints + '.png')
	plt.close()

	fig = plt.figure()
	plt.rc('axes', axisbelow=True)

	x = [flight_no[i] for i in range(len(labels1)) for j in range(len(d1[i])) if d1[i][j] == min(d1[i])]
	y = [d1[i][j] for i in range(len(d1)) for j in range(len(d1[i])) if d1[i][j] == min(d1[i])]
	c = [deltat1[i][j] for i in range(len(d1)) for j in range(len(d1[i])) if d1[i][j] == min(d1[i])]

	plt.axhline(np.mean(y), color='orange', linestyle='--', linewidth=1, label='Mean min. d: ' + str(round(np.mean(y), 2)))

	sc = plt.scatter(x, y, c=c, s=100, cmap=cm, label='Min. d')

	x = [flight_no[i] for i in range(len(labels1)) for j in range(len(d1[i])) if np.abs(deltat1[i][j]) == min(np.abs(deltat1[i]))]
	y = [d1[i][j] for i in range(len(d1)) for j in range(len(d1[i])) if np.abs(deltat1[i][j]) == min(np.abs(deltat1[i]))]
	c = [deltat1[i][j] for i in range(len(d1)) for j in range(len(d1[i])) if np.abs(deltat1[i][j]) == min(np.abs(deltat1[i]))]

	sc = plt.scatter(x, y, c=c, s=100, cmap=cm, marker='^', label='Min. dt')

	clb = plt.colorbar(sc)
	plt.axhline(np.mean(y), linestyle='--', linewidth=1, label='Mean min. delta t: ' + str(round(np.mean(y), 2)))
	clb.set_label('delta t [hr]', fontsize=15)
	plt.xlabel('Flight date', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.xticks(np.arange(1, len(d1) + 1, step = 1), labels=labels1, rotation=60)
	plt.grid(True)
	plt.legend(loc='best')
	plt.tight_layout()

	fig.savefig(fig_dir + 'results_minD+T_ep_' + re_endpoints + '.png')
	plt.close()


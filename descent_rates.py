from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from datetime import datetime
from scipy import interpolate
from astropy.io import ascii
import match_files as mf
import numpy as np
import csv
import sys
import os

def get_rates(next_point='0'):

	FMT = "%H:%M:%S"

	descent_rates_pred = {}
	descent_rates_gps = {}
	gps_indices = {}
	alts_pred = {}
	alts_gps = {}

	dir_pred = '/home/ellen/Desktop/SuperBIT/Weather_data/Trajectories/'

	if next_point == '2':
		ext = 'next2_points/'
	elif next_point == '1':
		ext = 'next_points/'
	else:
		ext = 'all_points/'

	dir_pred = dir_pred + ext
	dir_gps = '/home/ellen/Desktop/SuperBIT/Flight_data/'

	add = int(next_point)

	for filename in os.listdir(dir_pred):

		if filename.endswith('.dat'):

			data = ascii.read(dir_pred + filename)

			alts = np.array(data['alts'])
			times = np.array(data['times'])
			descent_speeds = np.array(data['descent_speeds'])

			if '20180406' in filename:
				if '0600' in filename:
					datestr_pred = filename[:8] + '_1'
				else:
					datestr_pred = filename[:8] + '_2'
			else:
				datestr_pred = filename[:8]

			descent_rates_pred[filename[:17]] = descent_speeds
			alts_pred[filename[:17]] = np.array(alts)

			datestr_gps = mf.match_pred2gps(datestr_pred)	
			gps_file = datestr_gps + '.csv'

			with open(dir_gps + gps_file, newline='') as csvfile:

				data_gps = np.array(list(csv.reader(csvfile)))

				alts = []
				times = []

				for i in range(len(data_gps)):
					if 'RB' in data_gps[i][0]:
						alts.append(int(data_gps[i][5]))
						times.append(data_gps[i][2])
					else:
						alts.append(int(data_gps[i][4]))
						times.append(data_gps[i][1])

				ind = int(np.where(alts == np.max(alts))[0]) + add

				alt0 = np.max(alts)

				for j in range(ind + 1, len(alts) - 1):
					if alts[j] < alt0:
						alt0 = alts[j]
					elif alts[j] >= alt0:
						j -= 1
						break

				gps_indices[filename[:17]] = j - ind
				descent_rates_gps[filename[:17]] = np.array([(alts[i+1] - alts[i])/(datetime.strptime(times[i+1], FMT) - datetime.strptime(times[i], FMT)).seconds for i in range(ind, len(alts) -1)])
				alts_gps[filename[:17]] = np.array([alts[i] for i in range(ind, len(alts) -1)])

	return (descent_rates_gps, descent_rates_pred, alts_gps, alts_pred), gps_indices, ext

def plot_results(data=None, next_point='0', ext='all_points'):

	gps_indices = None

	if data == None:
		(descent_rates_gps, descent_rates_pred, alts_gps, alts_pred), gps_indices, ext = get_rates(next_point=next_point)

	fig_dir = '/home/ellen/Desktop/SuperBIT/figs/DescentRates/' + ext

	if not os.path.exists(fig_dir):
	    os.makedirs(fig_dir)

	pred_keys = list(descent_rates_pred.keys())
	gps_keys = list(descent_rates_gps.keys())

	descent_rates_gps_vals = list(descent_rates_gps.values())
	descent_rates_pred_vals = list(descent_rates_pred.values())

	gps_indices_l = list(gps_indices.values())

	alts_gps_vals = list(alts_gps.values())
	alts_pred_vals = list(alts_pred.values())

	fig = plt.figure()
	plt.rc('axes', axisbelow=True)

	for i in range(len(descent_rates_gps_vals)):

		### x-coordinates at which to evaluate the interpolated values, x&y coordinates used in the interpolation

		descent_rates_interp_pred = np.interp(alts_gps_vals[i][::-1], alts_pred_vals[i][::-1], descent_rates_pred_vals[i][::-1]) 
		testx = np.arange(min(descent_rates_gps_vals[i]), max(descent_rates_gps_vals[i]), 0.5)

		# plt.plot(descent_rates_gps_vals[i], descent_rates_interp_pred[::-1], 'ro--', markersize=2, label='Results')
		plt.plot(descent_rates_gps_vals[i][:gps_indices_l[i]+1], descent_rates_interp_pred[::-1][:gps_indices_l[i]+1], 'ro--', markersize=2, label='Results')
		plt.plot(descent_rates_gps_vals[i][gps_indices_l[i]:], descent_rates_interp_pred[::-1][gps_indices_l[i]:], 'go--', markersize=2, label='Results - after end point')

		plt.plot(testx, testx, 'b--', linewidth=1, label=r'$v_{pred}=v_{true}$')

		plt.plot(descent_rates_gps_vals[i][gps_indices_l[i]], descent_rates_interp_pred[::-1][gps_indices_l[i]], 'bo', markersize=4, label='Endpoint GPS')
		# plt.axvline(descent_rates_gps_vals[i][gps_indices_l[i]], linestyle='--', linewidth=1, color='green')

		plt.ylabel(r'$v_{pred}$ [m/s]', fontsize=15)
		plt.xlabel(r'$v_{true}$ [m/s]', fontsize=15)
		plt.legend(loc='best')
		plt.grid(True)
		fig.savefig(fig_dir + 'vdescent_' + pred_keys[i] + '_interp.png')

		plt.clf()

		plt.plot(alts_gps_vals[i][:gps_indices_l[i]+1], descent_rates_gps_vals[i][:gps_indices_l[i]+1], 'ro--', markersize=2, label='True')
		plt.plot(alts_gps_vals[i][gps_indices_l[i]:], descent_rates_gps_vals[i][gps_indices_l[i]:], 'go--', markersize=2, label='True - after endpoint')

		plt.plot(alts_pred_vals[i], descent_rates_pred_vals[i], 'bo--', markersize=2, label='Predicted')

		plt.plot(alts_gps_vals[i][gps_indices_l[i]], descent_rates_gps_vals[i][gps_indices_l[i]], 'bo', markersize=4,  label='Endpoint GPS, ' + str(alts_gps_vals[i][gps_indices_l[i]]) + ' m')
		# plt.axvline(alts_gps_vals[i][gps_indices_l[i]], linestyle='--', linewidth=1, color='green')

		plt.ylabel('Descent rate [m/s]', fontsize=15)
		plt.xlabel('Altitude [m]', fontsize=15)
		plt.grid(True)
		plt.legend(loc='best')
		fig.savefig(fig_dir + 'vdescent_' + pred_keys[i] + '.png')

		plt.clf()

	minimum = np.inf
	pred_all_drates_tot = []

	for i in range(len(descent_rates_gps_vals)):

		if min(descent_rates_gps_vals[i]) < minimum:
			minimum = min(descent_rates_gps_vals[i])

		descent_rates_interp_pred = np.interp(alts_gps_vals[i][::-1], alts_pred_vals[i][::-1], descent_rates_pred_vals[i][::-1]) 
		pred_all_drates_tot.append(descent_rates_interp_pred[::-1])


		if i == 0:
			plt.plot(descent_rates_gps_vals[i], descent_rates_interp_pred[::-1], 'ro-', markersize=1, linewidth=0.5, label='Predictions')
			plt.plot(descent_rates_gps_vals[i][gps_indices_l[i]], descent_rates_interp_pred[::-1][gps_indices_l[i]], 'bo', markersize=4, label='End Point GPS')

		else:
			plt.plot(descent_rates_gps_vals[i], descent_rates_interp_pred[::-1], 'ro-', markersize=1, linewidth=0.5)
			plt.plot(descent_rates_gps_vals[i][gps_indices_l[i]], descent_rates_interp_pred[::-1][gps_indices_l[i]], 'bo', markersize=4)


	testx = np.arange(minimum, 0, 0.5)

	plt.plot(testx, testx, 'b--', linewidth=1, label=r'$v_{pred}=v_{true}$')
	plt.ylabel(r'$v_{pred}$ [m/s]', fontsize=15)
	plt.xlabel(r'$v_{true}$ [m/s]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best')

	fig.savefig(fig_dir + 'vdescent_all.png')

	plt.clf()

	gps_all_drates_tot = [rates for rates in descent_rates_gps_vals]

	gps_all_drates = [gps_all_drates_tot[i][j] for i in range(len(gps_all_drates_tot)) for j in range(len(gps_all_drates_tot[i])) if j < gps_indices_l[i]]
	pred_all_drates = [pred_all_drates_tot[i][j] for i in range(len(pred_all_drates_tot)) for j in range(len(pred_all_drates_tot[i])) if j < gps_indices_l[i]]

	gps_all_drates_rej = [gps_all_drates_tot[i][j] for i in range(len(gps_all_drates_tot)) for j in range(len(gps_all_drates_tot[i])) if j > gps_indices_l[i]]
	pred_all_drates_rej = [pred_all_drates_tot[i][j] for i in range(len(pred_all_drates_tot)) for j in range(len(pred_all_drates_tot[i])) if j > gps_indices_l[i]]

	epx = [gps_all_drates_tot[i][gps_indices_l[i]] for i in range(len(gps_all_drates_tot))]
	epy = [pred_all_drates_tot[i][gps_indices_l[i]] for i in range(len(pred_all_drates_tot))]

	def func(x, m, c):
		return m*x + c

	lim = np.inf

	popt, pcov = curve_fit(func, np.array(gps_all_drates)[np.array(gps_all_drates) <= lim,...], np.array(pred_all_drates)[np.array(gps_all_drates) <= lim,...])

	plt.plot(gps_all_drates, func(np.array(gps_all_drates), popt[0], popt[1]), 'm-', linewidth=1, label='Best Fit')
	plt.axvline(lim, linestyle='--', linewidth=1)
	
	plt.plot(testx, testx, 'b--', linewidth=1, label=r'$v_{pred}=v_{true}$')
	plt.scatter(gps_all_drates_rej, pred_all_drates_rej, color='green', s=1, label='After endpoint')
	plt.scatter(gps_all_drates, pred_all_drates, color='red', s=2, label='Before endpoint')
	plt.plot(epx, epy, 'bo', markersize=1, label='Endpoints')

	plt.xlabel(r'$v_{true}$ [m/s]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best')

	fig.savefig(fig_dir + 'vdescent_all2.png')

	plt.clf()

	return popt, pcov

if __name__ == '__main__':
	next_point = input('Start at point after max. altitude: ')
	popt, pcov = plot_results(next_point=next_point)


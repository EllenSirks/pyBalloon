from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from datetime import datetime
from scipy import interpolate
from astropy.io import ascii
import numpy as np
import csv
import sys
import os

import param_file as p
import pyb_aux


alt_err = 2.

def get_rates(params=None):

	FMT = "%H:%M:%S"

	descent_rates_pred = {}
	descent_rates_gps = {}
	gps_indices = {}
	alts_pred = {}
	alts_gps = {}

	dir_pred = '/home/ellen/Desktop/SuperBIT/Weather_data/Trajectories/'

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

	if descent_only:
		ext = 'descent_only/start_point' + next_point + '/'
	else:
		ext = 'ascent+descent/'

	dir_pred += ext
	dir_gps = '/home/ellen/Desktop/SuperBIT/Flight_data/'

	add = int(next_point)

	for filename in os.listdir(dir_pred):
		if filename.endswith(str(int(drift_time)).zfill(4) + 'min.dat') and '20180811' not in filename and '20180803' not in filename and 'interpolate' not in filename:
		# if ilename.endswith(str(int(drift_time)).zfill(4) + 'min.dat') and 'interpolate' not in filename:

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

			descent_rates_pred[datestr_pred] = descent_speeds[1:]
			alts_pred[datestr_pred] = np.array(alts)[1:]

			datestr_gps = pyb_aux.match_pred2gps(datestr_pred)
			gps_file = datestr_gps + '.csv'

			with open(dir_gps + gps_file) as csvfile:

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

				ind = int(np.where(alts == np.max(alts))[0])
				ind += add

				alt0 = np.max(alts)

				descent_rates_gps[datestr_pred] = np.array([float((alts[i+1] - alts[i]))/(datetime.strptime(times[i+1], FMT) - datetime.strptime(times[i], FMT)).seconds for i in range(ind, len(alts)-2)])
				alts_gps[datestr_pred] = np.array([alts[i] for i in range(ind, len(alts) -2)])

	return (descent_rates_gps, descent_rates_pred, alts_gps, alts_pred)

def plot_rates(data=None, params=None):

	if data == None:
		(descent_rates_gps, descent_rates_pred, alts_gps, alts_pred) = get_rates()

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

	if descent_only:
		ext = 'descent_only/start_point' + next_point + '/'
	else:
		ext = 'ascent+descent/'

	fig_dir = '/home/ellen/Desktop/SuperBIT/figs/DescentRates/' + ext + str(drift_time) + 'drift/'

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
		descent_rates_interp_pred = np.interp(alts_gps_vals[i][::-1], alts_pred_vals[i][::-1], descent_rates_pred_vals[i][::-1]) 
		testx = np.arange(-0.5 + min(descent_rates_gps_vals[i]), max(descent_rates_gps_vals[i]) + 0.5, 0.5)

		plt.plot(descent_rates_gps_vals[i], descent_rates_interp_pred[::-1], 'ro--', markersize=2, label='Results')

		plt.plot(testx, testx, 'b--', linewidth=1, label=r'$v_{pred}=v_{true}$')

		plt.ylabel(r'$v_{pred}$ [m/s]', fontsize=15)
		plt.xlabel(r'$v_{true}$ [m/s]', fontsize=15)
		plt.legend(loc='best')
		plt.grid(True)
		fig.savefig(fig_dir + 'vdescent_' + pred_keys[i] + '_interp.png')

		plt.clf()

		plt.plot(alts_pred_vals[i], descent_rates_pred_vals[i], 'bo--', markersize=2, label='Predicted')
		plt.plot(alts_gps_vals[i], descent_rates_gps_vals[i], 'ro--', markersize=2, label='True')
		# plt.errorbar(x=alts_gps_vals[i], y=descent_rates_gps_vals[i], xerr=alt_err, fmt='ro', markersize=2, label='True')

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
			# plt.errorbar(x=descent_rates_gps_vals[i], y=descent_rates_interp_pred[::-1], xerr=alt_err*np.sqrt(2), fmt='ro-', markersize=1, linewidth=0.5, label='Predictions')

		else:
			plt.plot(descent_rates_gps_vals[i], descent_rates_interp_pred[::-1], 'ro-', markersize=1, linewidth=0.5)
			# plt.errorbar(x=descent_rates_gps_vals[i], y=descent_rates_interp_pred[::-1], xerr=alt_err*np.sqrt(2), fmt='ro-', markersize=1, linewidth=0.5)


	testx = np.arange(minimum, 0, 0.5)

	plt.plot(testx, testx, 'b--', linewidth=1, label=r'$v_{pred}=v_{true}$')
	plt.ylabel(r'$v_{pred}$ [m/s]', fontsize=15)
	plt.xlabel(r'$v_{true}$ [m/s]', fontsize=15)
	plt.grid(True)
	plt.legend(loc='best')

	fig.savefig(fig_dir + 'vdescent_all.png')

	plt.clf()

	# return

def plot_rho(descent_only=True, next_point='0'):

	descent_only = p.descent_only
	
	if descent_only:
		next_point = p.next_point
		ext = 'descent_only/start_point' + next_point + '/'
	else:
		ext = 'ascent+descent/'


	dir = '/home/ellen/Desktop/SuperBIT/properties/' + ext
	fig_dir = '/home/ellen/Desktop/SuperBIT/figs/properties/' + ext

	fig = plt.figure()

	for filename in os.listdir(dir):

		if filename.startswith('prop_preinterp'):

			fpre = dir + filename
			fafter = dir + 'prop_afterinterp' + filename[14:]

			name = filename[14:-4]

			data_pre = ascii.read(fpre)
			data_after = ascii.read(fafter)

			alt_pre = data_pre['alt']
			rho_pre = data_pre['rho']
			u_pre = data_pre['u']
			v_pre = data_pre['v']
			T_pre = data_pre['T']

			alt_after = data_after['alt']
			rho_after = data_after['rho']
			u_after = data_after['u']
			v_after = data_after['v']
			T_after = data_after['T']

			alt0 = data_after['alt0'][0]
			grids = data_after['grid']

			plt.plot(alt_after, rho_after, 'bo', markersize=1, label='After Interpolation')
			plt.plot(alt_pre, rho_pre, 'ro', markersize=3, label='Before Interpolation')
			plt.xlabel('Altitude [m]', fontsize=15)
			plt.ylabel(r'Density [kg $m^{-3}$]', fontsize=15)

			plt.axvline(alt0, linestyle='--', linewidth=1, label='alt0 = ' + str(alt0) + ' m')

			plt.ylim([0, 1.4])
			plt.xlim([0, 50000])

			plt.grid(True)
			plt.legend(loc='best')
			plt.tight_layout()
			fig.savefig(fig_dir + 'rho/' + 'rho_check' + str(name) + '.png')

			plt.clf()

			plt.plot(alt_after, u_after, 'bo', markersize=1, label='After Interpolation')
			plt.plot(alt_pre, u_pre, 'ro', markersize=3, label='Before Interpolation')
			plt.xlabel('Altitude [m]', fontsize=15)
			plt.ylabel(r'U wind', fontsize=15)

			plt.axvline(alt0, linestyle='--', linewidth=1, label='alt0 = ' + str(alt0) + ' m')

			plt.xlim([0, 50000])

			plt.grid(True)
			plt.legend(loc='best')
			plt.tight_layout()
			fig.savefig(fig_dir + 'u_wind/' + 'u_check' + str(name) + '.png')

			plt.clf()

			plt.plot(alt_after, v_after, 'bo', markersize=1, label='After Interpolation')
			plt.plot(alt_pre, v_pre, 'ro', markersize=3, label='Before Interpolation')
			plt.xlabel('Altitude [m]', fontsize=15)
			plt.ylabel(r'V wind', fontsize=15)

			plt.axvline(alt0, linestyle='--', linewidth=1, label='alt0 = ' + str(alt0) + ' m')

			plt.xlim([0, 50000])

			plt.grid(True)
			plt.legend(loc='best')
			plt.tight_layout()
			fig.savefig(fig_dir + 'v_wind/' + 'v_check' + str(name) + '.png')

			plt.clf()

			plt.plot(alt_after, T_after, 'bo', markersize=1, label='After Interpolation')
			plt.plot(alt_pre, T_pre, 'ro', markersize=3, label='Before Interpolation')
			plt.xlabel('Altitude [m]', fontsize=15)
			plt.ylabel(r'Temperature [K]', fontsize=15)

			plt.axvline(alt0, linestyle='--', linewidth=1, label='alt0 = ' + str(alt0) + ' m')

			plt.xlim([0, 50000])

			plt.grid(True)
			plt.legend(loc='best')
			plt.tight_layout()
			fig.savefig(fig_dir + 'temperature/' +'T_check' + str(name) + '.png')

			plt.clf()

	return

if __name__ == '__main__':
	plot_rates()
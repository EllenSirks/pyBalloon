import matplotlib.pyplot as plt
from astropy.io import ascii
import numpy as np
import sys
import os

def plot_rho(next_point='0'):

	if next_point == '2':
		ext = 'next2_points/'
	elif next_point == '1':
		ext = 'next_points/'		
	else:
		ext = 'all_points/'

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
	next_point = input('Start at point after max. altitude: ')
	plot_rho(next_point=next_point)


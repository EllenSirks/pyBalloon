import matplotlib.pyplot as plt
from haversine import haversine
from astropy.io import ascii
import match_files as mf
import numpy as np
import csv
import os

next_point = '1'
re_endpoints = 'no'
interpolate = False

if interpolate:
	interp = '_interpolated_'
else:
	interp = ''

date_keys = ['20180304', '20180406', '20181512', '20190101', '20180813','20180803', '20180506', '20180616',  '20180630', '20180511', '20180519']

ext = 'start_point' + next_point + '/'

dir = '/home/ellen/Desktop/SuperBIT/'
dir_pred =  dir + 'Weather_data/Endpoints/' + ext
dir_gps =  dir + 'Flight_data/'
fig_dir = dir + 'figs/Results/' + ext 

d = {}
for filename in os.listdir(dir_pred):

	if filename.endswith('0000min.dat') and filename[9:17] in date_keys:

		if re_endpoints == 'yes' and '201808' not in filename:
			end_lat_pred, end_lon_pred = ep.get_endpoint(filename[:-13] + '_trajectory.dat')

		elif re_endpoints == 'no' or '201808' in filename:
			data_pred = ascii.read(dir_pred + filename)

			end_lat_pred = data_pred['lat'][-1]
			end_lon_pred = data_pred['lon'][-1]

		if '20180406' in filename:
			if '8.9' in filename:
				datestr_pred = filename[9:17] + '_1'
			else:
				datestr_pred = filename[9:17] + '_2'
		else:
			datestr_pred = filename[9:17]

		datestr_gps = mf.match_pred2gps(datestr_pred)
		gps_file = datestr_gps + '.csv'

		with open(dir_gps + gps_file, newline='') as csvfile:

			data_gps = np.array(list(csv.reader(csvfile)))

			ext ='_ep_' + re_endpoints 
			j = -1

			if 'RB' in data_gps[j][0]:
				d.setdefault(datestr_pred, []).append(haversine((end_lat_pred, end_lon_pred), (float(data_gps[j][3]), float(data_gps[j][4]))))
			else:
				d.setdefault(datestr_pred, []).append(haversine((end_lat_pred, end_lon_pred), (float(data_gps[j][2]), float(data_gps[j][3]))))

labels = list(d.keys())
flight_no = np.arange(1, len(labels) + 1, 1)
vals = list(d.values())

vals = [x for _,x in sorted(zip(labels,vals))]
mean = np.mean(vals)
labels.sort()

fig = plt.figure()

plt.plot(flight_no, vals, 'ro')
plt.axhline(mean, linestyle='--', color='black', linewidth=1, label='Average error: ' + str(round(mean, 1)) + ' km')
plt.xticks(flight_no, labels, rotation='vertical')
plt.xlabel('Flight date', fontsize=15)
plt.ylabel('Error in distance [km]', fontsize=15)
plt.legend(loc='best')
plt.grid(True)
plt.tight_layout()
plt.show()

fig.savefig(fig_dir + 'results' + interp + '_startpoint' + next_point +'.png')
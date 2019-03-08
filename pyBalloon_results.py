import matplotlib.pyplot as plt
from haversine import haversine
from astropy.io import ascii
import numpy as np
import csv
import sys
import os

dir_pred = '/home/ellen/Desktop/SuperBIT/Weather_data/Endpoints/'
dir_res = '/home/ellen/Desktop/SuperBIT/Flight_data/'
fig_dir = '/home/ellen/Desktop/SuperBIT/figs/'

def calc_results(trunc):

	endpts_res = {}
	endpts_pred = {}

	for filename in os.listdir(dir_pred):

		data = ascii.read(dir_pred + filename)
		endpts_pred[filename[:17]] = (data['lat'][1], data['lon'][1])

	for filename in os.listdir(dir_res):

		if filename.endswith(".csv"):

			for string in range(len(filename)):
				if filename[string] == '.':
					s = string

			with open(dir_res + filename, newline='') as csvfile:

				data = np.array(list(csv.reader(csvfile)))

				alts = []

				if trunc == 'False':

					ext = ''

					if 'RB' in data[-1][0]:

						endpts_res[filename[:s]] = (float(data[-1][3]), float(data[-1][4]))
					else:

						endpts_res[filename[:s]] = (float(data[-1][2]), float(data[-1][3]))

				elif trunc == 'True':

					ext = '_trunc'

					for i in range(len(data)):
						if 'RB' in data[i][0]:
							alts.append(int(data[i][5]))
						else:
							alts.append(int(data[i][4]))

					alt0 = np.max(alts)
					ind, = np.where(alts == np.max(alts))[0]

					for j in range(len(alts)):
						if j > ind and alts[j] < alt0:
							alt0 = alts[j]
						elif j > ind and alts[j] > alt0:
							break

					if 'RB' in data[-1][0]:
						endpts_res[filename[:s]] = (float(data[j][3]), float(data[j][4]))
					else:
						endpts_res[filename[:s]] = (float(data[j][2]), float(data[j][3]))

	d = []
	matched_keys = {}

	for key1 in endpts_pred.keys():

		for key2 in endpts_res.keys():

			if key1[2:4] == key2[:2]:

				if not '_' in key2:

					if key2[4] == '-':

						if len(key2) == 6: ##-#-#

							if key1[4:6] == '0' + key2[3] and key1[6:8] == '0' + key2[5]:

								d.append(haversine(endpts_pred[key1], endpts_res[key2]))
								key2 = key2.replace('-', '0')
								key2 = key2.replace(key2[:2], '20' + key2[:2])
								matched_keys[key2] = key1


						else: ##-#-##

							if key1[4:6] == '0' + key2[3] and key1[6:8] == key2[5:]:

								d.append(haversine(endpts_pred[key1], endpts_res[key2]))
								key2 = key2.replace('-' + key2[3], '0' + key2[3])
								key2 = key2.replace('-', '')
								key2 = key2.replace(key2[:2], '20' + key2[:2])
								matched_keys[key2] = key1


					elif key2[5] == '-':

						if len(key2) == 7: ##-##-#

							if key1[4:6] == key2[3:5] and key1[6:8] == '0' + key2[6]:

								d.append(haversine(endpts_pred[key1], endpts_res[key2]))
								key2 = key2.replace('-' + key2[6], '0' + key2[6])
								key2 = key2.replace('-', '')
								key2 = key2.replace(key2[:2], '20' + key2[:2])
								matched_keys[key2] = key1

						else: ##-##-##

							if key1[4:6] == key2[3:5] and key1[6:8] == key2[6:]: 

								d.append(haversine(endpts_pred[key1], endpts_res[key2]))
								key2 = key2.replace('-', '')
								key2 = key2.replace(key2[:2], '20' + key2[:2])
								matched_keys[key2] = key1

				elif '_' in key2:

					if key1 == '20180406_0600_000' and key2 == '18-4-6_1':

						d.append(haversine(endpts_pred[key1], endpts_res[key2]))
						key2 = key2.replace('-', '0')
						key2 = key2.replace(key2[:2], '20' + key2[:2])
						matched_keys[key2] = key1

					elif key1 == '20180406_1800_000' and key2 == '18-4-6_2':

						d.append(haversine(endpts_pred[key1], endpts_res[key2]))
						key2 = key2.replace('-', '0')
						key2 = key2.replace(key2[:2], '20' + key2[:2])
						matched_keys[key2] = key1

	flight_no = np.arange(1, len(d)+1, 1)

	labels = list(matched_keys.keys())

	ds = [x for _,x in sorted(zip(labels, d))]
	labels.sort()

	fig = plt.figure()
	plt.rc('axes', axisbelow=True)
	plt.axhline(np.mean(ds), linestyle='--', linewidth=1, label='Average')
	plt.plot(flight_no, ds, 'ro')
	plt.xlabel('Flight date', fontsize=15)
	plt.ylabel('Error in distance [km]', fontsize=15)
	plt.xticks(np.arange(1, len(ds) + 1, step = 1 ), labels=labels, rotation=60)
	plt.grid(True)
	plt.legend(loc='best')
	plt.tight_layout()
	fig.savefig(fig_dir + 'results' + ext + '.png')

	return d, matched_keys

d1, labels1 = calc_results('True')
d2, labels2 = calc_results('False')

labels1 = list(labels1.keys())
labels2 = list(labels2.keys())

s = sorted(zip(labels1, d1))
labels1, d1 = map(list, zip(*s))

s = sorted(zip(labels2, d2))
labels2, d2 = map(list, zip(*s))

flight_no = np.arange(1, len(d1) + 1, 1)

fig = plt.figure()
plt.rc('axes', axisbelow=True)
plt.axhline(np.mean(d1), linestyle='--', color='red', linewidth=1, label='Average - truncated')
plt.axhline(np.mean(d2), linestyle='--', color='blue', linewidth=1, label='Average - not truncated')
plt.plot(flight_no, d1, 'ro', label='Truncated')
plt.plot(flight_no, d2, 'bo', label='Not Truncated')
plt.xlabel('Flight date', fontsize=15)
plt.ylabel('Error in distance [km]', fontsize=15)
plt.xticks(np.arange(1, len(d1) + 1, step = 1), labels=labels1, rotation=90)
plt.grid(True)
plt.legend(loc='best')
plt.tight_layout()
fig.savefig(fig_dir + 'results_trunc_vs_untrunc.png')

fig = plt.figure()
plt.rc('axes', axisbelow=True)

d1 = np.array(d1)
d2 = np.array(d2)
diff = d2 - d1

plt.plot(flight_no, diff, 'ro', label='Diff = not truncated - truncated')
plt.xlabel('Flight date', fontsize=15)
plt.ylabel('Difference in error [km]', fontsize=15)
plt.xticks(np.arange(1, len(d1) + 1, step = 1), labels=labels1, rotation=90)
plt.grid(True)
plt.legend(loc='best')
plt.tight_layout()
fig.savefig(fig_dir + 'diff_trunc_vs_untrunc.png')
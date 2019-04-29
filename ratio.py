from match_files import match_pred2gps, match_gps2pred
import numpy as np
import csv
import sys
import re

def get_ratio(date, point):

	dir_res = '/home/ellen/Desktop/SuperBIT/Flight_data/'
	filename = date + '.csv'

	with open(dir_res + filename) as csvfile:

		data = np.array(list(csv.reader(csvfile)))

		alts = []
		lons = []
		lats = []

		for i in range(len(data)):
			if 'RB' in data[i][0]:
				alts.append(int(data[i][5]))
				lons.append(float(data[i][4]))
				lats.append(float(data[i][3]))
			else:
				alts.append(int(data[i][4]))
				lons.append(float(data[i][3]))
				lats.append(float(data[i][2]))

		alt0 = np.max(alts)
		ind, = np.where(alts == np.max(alts))[0]

		start = lats[ind + point], lons[ind + point], alts[ind + point]
		end = lats[ind + point+1], lons[ind + point+1], alts[ind + point+1]

		Dlat = end[0] - start[0] # positive: moving north
		Dlon = end[1] - start[1] # positive: moving east
		f = np.abs(Dlat/Dlon)

	return (Dlat, Dlon), f

def additives(Dlat=None, Dlon=None, f=None, date=None, point=0):

	factor = 1000.

	if date != None and Dlat == None and Dlon == None and f == None:
		datestr_gps = match_pred2gps(date)
		(Dlat, Dlon), f = get_ratio(date=datestr_gps, point=point)

	ex1 = (Dlon/np.abs(Dlon))*factor
	ex2 = (Dlat/np.abs(Dlat))*f*np.abs(ex1)

	return ex1, ex2

if __name__ == '__main__':
	(Dlat, Dlon), f = get_ratio(sys.argv[1])
	ex1, ex2 = additives(Dlat, Dlon, f)
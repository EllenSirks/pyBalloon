from astropy.io import ascii
import numpy as np
import pyb_io
import time
import sys
import os

time0 = time.time()

# next_point = input('Starting point: ')
# datestr = input('Date string: ')
# interpolate = input('Interpolated: ')
# if interpolate == 'Yes':
# 	interpolated = True
# else:
# 	interpolated = False

next_point = '2'
datestr = sys.argv[1]
interpolated = False
flight_nr = 1
ext_str = '.'

if datestr == '20180406':
	flight_nr = int(input('flight 1 or 2? '))
	if flight_nr == 1:
		ext_str = '8.9'
	else:
		ext_str = '18.6'

dir_base = '/home/ellen/Desktop/SuperBIT/Weather_data/'
ext = 'start_point' + next_point + '/'
dir1 = dir_base + 'kml_files/' + ext

lines = {}
for fname in os.listdir(dir1):
	if datestr in fname and fname.endswith('min.kml') and not 'merged' in fname:
		if interpolated and 'interpolated' in fname or not interpolated and not 'interpolated' in fname:
			if ext_str in fname:
				fname0 = fname
				lines[int(fname[-10:-7])] =  [line.rstrip('\n').split(' ') for line in open(dir1 + fname)]
		
keys = list(lines.keys())

kml_str1 = '<?xml version="1.0" encoding="UTF-8"?>\n'
kml_str1 += '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">\n'
kml_str1 += '<Document id="feat_2">\n'

for key in keys:
	info = lines[key]
	for i in range(len(info)):
		if info[i] == ['<name>pyballoon', 'trajectory</name>']:
			ind1 = i
		elif info[i] ==  ['</Document>']:
			ind2 = i - 1
	for i in range(len(info)):
		if i >= ind1 and i <= ind2:
			for j in range(len(info[i])):
				if j == len(info[i]) - 1:
					kml_str1 += str(info[i][j]) + '\n'
				else:
					kml_str1 += str(info[i][j]) + ' '

kml_str1 += '</Document>\n'
kml_str1 += '</kml>\n'

kml_fname1 = 'merged_' + fname0[:-13] + '_' + str(min(keys)) + '-' + str(max(keys)) + 'min.kml'

fid = open(dir1 + kml_fname1, 'w')
fid.write(kml_str1)
fid.close()

dir2 = dir_base + 'Endpoints/' + ext

endpoints = {}
for fname in os.listdir(dir2):
	if datestr in fname:
		if interpolated and 'interpolated' in fname or not interpolated and not 'interpolated' in fname:
			data = ascii.read(dir2 + fname)
			endpoint = [data['lat'][1], data['lon'][1], data['alt'][1]]
			endpoints[int(fname[-11:-7])] = endpoint

kml_fname2 = 'endpoints_merged_' + fname0[:-13] + '_' + str(min(keys)) + '-' + str(max(keys)) + 'min.kml'

kml_str2 = '<?xml version="1.0" encoding="UTF-8"?>\n'
kml_str2 += '<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">\n'
kml_str2 += '<Document id="feat_2">\n'
kml_str2 += '<name>pyballoon trajectory</name>\n'
kml_str2 += '<Style id="stylesel_362">\n'
kml_str2 += '<LineStyle id="substyle_363">\n'
kml_str2 += '<color>BF0000DF</color>\n'
kml_str2 += '<colorMode>normal</colorMode>\n'
kml_str2 += '<width>5</width>\n'
kml_str2 += '</LineStyle>\n'
kml_str2 += '<PolyStyle id="substyle_364">\n'
kml_str2 += '<color>BF0000DF</color>\n'
kml_str2 += '<colorMode>normal</colorMode>\n'
kml_str2 += '<fill>1</fill>\n'
kml_str2 += '<outline>1</outline>\n'
kml_str2 += '</PolyStyle>\n'
kml_str2 += '</Style>\n'
kml_str2 += '<Placemark id="feat_91">\n'
kml_str2 += '<styleUrl>#stylesel_362</styleUrl>\n'
kml_str2 += '<LineString id="geom_86">\n'
kml_str2 += '<coordinates>\n'

keys = list(endpoints.keys())
keys.sort()

for key in keys:
	kml_str2 += str(endpoints[key][1]) + ',' + str(endpoints[key][0])  + ',' + str(endpoints[key][2]) + '\n'

kml_str2 += '</coordinates>\n'
kml_str2 += '<extrude>1</extrude>\n'
kml_str2 += '<tessellate>1</tessellate>\n'
kml_str2 += '<altitudeMode>relativeToGround</altitudeMode>\n'
kml_str2 += '</LineString>\n'
kml_str2 += '</Placemark>\n'

for key in keys:
	kml_str2 += '<Placemark>\n'
	kml_str2 += '<name>End, drift = '+ str(key) + ' min.</name>\n'
	kml_str2 += '<Point>\n'
	kml_str2 += '<coordinates>' + str(endpoints[key][1]) + ',' + str(endpoints[key][0]) + '</coordinates>' + '\n'
	kml_str2 += '</Point>\n'
	kml_str2 += '</Placemark>\n'

	kml_str_add = pyb_io.create_circle(lon = endpoints[key][1], lat = endpoints[key][0])
	kml_str2 += kml_str_add

kml_str2 += '</Document>\n'
kml_str2 += '</kml>'

fid = open(dir1 + kml_fname2, 'w')
fid.write(kml_str2)
fid.close()

print('Program finished in %.3f s' % (time.time() - time0))

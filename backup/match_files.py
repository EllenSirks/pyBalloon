import re

def match_pred2gps(date):

	datestr = date[2:4] + '-' + date[4:6] + '-' + date[6:]

	if datestr[3] == '0':
		datestr = datestr[:3] + datestr[4:]
		if datestr[5] == '0':
			datestr = datestr[:5] + datestr[6:]
	else:
		if datestr[6] == '0':
			datestr = datestr[:6] + datestr[7:]

	return datestr

def match_gps2pred(date):

	if '_' in date:
		i0 = [m.start() for m in re.finditer('_', date)][0]
	else:
		i0 = len(date)

	inds = [m.start() for m in re.finditer('-', date)]
	i1, i2 = inds[0], inds[1]

	if i2 - i1 == 2:
		month = '0' + date[i1+1]
	else:
		month = date[i1+1:i1+3]

	if i0 - i2 == 2:
		day = '0' + date[i2+1:]
	else:
		day =  date[i2+1:]

	return '20' + date[:2] + month + day
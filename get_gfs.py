#!/usr/bin/python

import sys
import requests
import datetime as dt
import numpy as np
import os

def get_gfs_data(datestr, utc_hour, verbose = False):

    """

    Download GFS data from NOMADS-NOAA for requested time and area.

    Required arguments:

        - datestr -- requested date as a string, eg. '20131803'
        - utc_hour -- requested UTC hour with two digits as a string
        - area -- 4-tuple of coordinates limiting the area to be retrieved.

    Coordinates are given as floats: (blat, tlat, llon, rlon), where:

        - blat -- bottom latitude (southern) limit (float)
        - tlat -- top latitude (northern) limit (float)
        - llon -- left longitude (western) limit (float)
        - rlon -- right longitude (eastern) limit (float)

    Latitudes are north-positive, longitudes east-positive. Limits are
    inclusive. Data from closest available time will be downloaded for
    the latest available GFS main run, EPS main and 20 EPS members.

    """

    # Read requested date, time and  area

    yyyy, mm, dd = int(datestr[0:4]), int(datestr[4:6]), int(datestr[6:8])
    request_time = dt.datetime(yyyy, mm, dd, int(utc_hour))

    out_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/'

    # GFS Ensemble Forecasts (1 degree grid)

    url = 'https://nomads.ncdc.noaa.gov/data/gfs4'

    req = requests.get(url)
    if req.ok is not True:
        print("Could not connect! Error code: " % req.error)
        sys.exit()

    text = req.content.split('/</a>')

    available_months = []
    for t in text:
        if ">20" in t:
            available_months.append(t.split('>')[-1])

    if str(datestr[0:6]) not in available_months:
        print('Date not available. No data from this month.')
        sys.exit()

    month = datestr[0:6]

    # GFS Ensemble Forecasts (1 degree grid)

    url_base = 'https://nomads.ncdc.noaa.gov/data/gfs4/'

    url = url_base + month
    req = requests.get(url)

    if req.ok is not True:
        print("Could not connect! Error code: " % req.error)
        sys.exit()

    text = req.content.split('/</a>')

    available_days = []
    for t in text:
        if '>20' in t:
            available_days.append(t.split('>')[-1].split('>')[-1])

    if str(datestr) not in available_days:
        print('Date not available. No data from this day.')
        sys.exit()

    day = datestr

    url = url_base + month + '/' + day
    req = requests.get(url)

    if req.ok is not True:
        print("Could not connect! Error code: " % req.error)
        sys.exit()

    text = req.content.split('</a>')

    available_files = []

    for t in text:
        if 'grb2' in t and not 'align' in t:
        # if ('grb2' in t or 'inv' in t) and not 'align' in t:
            available_files.append(t.split('>')[-1].split('>')[-1])

    good_hhhh2 = None
    good_hhh2 = None

    for file in available_files:

        if good_hhhh2 is not None and good_hhh2 is not None:
            break

        yyyy2, mm2, dd2, hhhh2, hhh2 = int(file[6:10]), int(file[10:12]), int(file[12:14]), int(int(file[15:19])/100) , int(file[20:23])

        init_dt = dt.datetime(yyyy2, mm2, dd2, hhhh2)
        delta_t = request_time - init_dt

        if delta_t.seconds < 0 or delta_t.days < 0 or np.abs(hhhh2 - int(utc_hour)) > 4 or np.abs(hhhh2 + hhh2 - int(utc_hour)) > 2:
            continue

        url = url_base + month + '/' + day + '/' + file

        if not os.path.isfile('/home/ellen/Desktop/SuperBIT/Weather_data/' + file):

            print('Beginning file download...')

            req = requests.get(url)

            if req.ok is not True:
                print("Could not connect! Error code: " % req.error)
                sys.exit()

            fid = open(out_dir+file, 'wb')
            fid.write(req.content)
            fid.close()

        else:
            print('File already downloaded!')

        good_hhhh2 = file[15:19]
        good_hhh2 = file[20:23]

    return good_hhhh2, good_hhh2

if __name__ == '__main__':

    """

    Command-line interface for downloading GFS data from NOMADS-NOAA for requested time and area.

    Usage: python get_gfs.py <yyyymmdd> <utc_hour> <blat,tlat,llon,rlon>

    where:

        - yyyymmdd -- requested date
        - utc_hour -- requested UTC hour with two digits
        - blat -- bottom latitude (southern) limit
        - tlat -- top latitude (northern) limit
        - llon -- left longitude (western) limit
        - rlon -- right longitude (eastern) limit

    Latitudes are north-positive, longitudes east-positive. Limits are 
    inclusive. Data from closest available time will be downloaded for
    the latest available GFS main run, EPS main and 20 EPS members.

    Example: python get_gfs.py 20181024 17 44.0,48.0,4.0,8.0

    When spanning the meridian, use: python get_gfs.py 20181024 17 50.0,60.0,355.0,365.0

    """

    datestr = sys.argv[1]
    utc_hour = sys.argv[2]

    get_gfs_data(datestr, utc_hour, verbose=True)

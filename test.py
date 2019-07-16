from scipy import interpolate
import datetime as dt
import pygrib as pg
import numpy as np

import sys
import os
 
import matplotlib.pyplot as plt

import param_file as p
import get_gfs
import pyb_aux
import pyb_io

g_0 = 9.80665 # m/s surface acc.
R0 = 8.3144621 # Ideal gas constant, J/(mol*K)
R_e = 6371009 # mean Earth radius in meters
M_air = 0.0289644 # molar mass of air [kg/mol], altitude dependence
T0 = 288.15 # K

def read_gefs_file(fname=None, area=None, alt0=0, t_0=None, extra_data=None, descent_only=False, step=100):

    indir = '/home/ellen/Desktop/SuperBIT/Weather_data/GEFS/'

    if area is not None:
        tlat, llon, blat, rlon = area
    else:
        print('Do you really wish to search the entire planet?')
        tlat, llon, blat, rlon = 90., 0., -90., 360.

    grib = pg.open(indir + fname)
    grib.seek(0)
    u_msgs = grib.select(name='U component of wind')
    v_msgs = grib.select(name='V component of wind')
    g_msgs = grib.select(name='Geopotential Height')

    lats, lons, = u_msgs[0].latlons() # lats: -90 -> 90, lons: 0 -> 360
    lats2, lons2 = u_msgs[0].latlons()

    for i in range(len(lons2)):
        for j in range(len(lons2[i])):
            if lons2[i][j] > 180:
                lons2[i][j] -= 360

    locs = pyb_aux.all_and([lats <= tlat, lats >= blat, lons2 <= rlon, lons2 >= llon])

    row_idx, col_idx = np.where(locs)
    lats = lats[row_idx, col_idx]
    lons = lons[row_idx, col_idx]

    if len(lats) == 0: print( 'Warning! lats is empty!')
    if len(lons) == 0: print( 'Warning! lons is empty!')

    u_wind, v_wind, altitude = {}, {}, {}
    for msg in u_msgs:
        if msg.typeOfLevel == 'isobaricInhPa':
            u_wind[msg.level] = msg.values[row_idx, col_idx]
    
    for msg in v_msgs:
        if msg.typeOfLevel == 'isobaricInhPa':
            v_wind[msg.level] = msg.values[row_idx, col_idx]

    for msg in g_msgs:
        if msg.typeOfLevel == 'isobaricInhPa':
            altitude[msg.level] = msg.values[row_idx, col_idx]

    pressures = list(u_wind.keys())
    pressures.sort()
    pressures.reverse()

    alt_keys, u_winds, v_winds, altitudes, alt_interp = [], [], [], [], []

    for key in pressures:
        uwnd, vwnd, alt = [], [], []
        uwnd.append(u_wind[key])
        vwnd.append(v_wind[key])

        u_winds.append(np.hstack(uwnd))
        v_winds.append(np.hstack(vwnd))

        if key in altitude.keys():

            alt_keys.append(key)
            alt.append(altitude[key])
            altitudes.append(np.hstack(alt))

    p_interp_val = list(set(pressures).symmetric_difference(set(alt_keys)))[0]
    
    for i in range(len(lats)):
        f = interpolate.interp1d(alt_keys, np.array(altitudes)[:, i])
        val = float(f(p_interp_val))
        alt_interp.append(val)

    altitudes.insert(np.where(np.array(pressures) == p_interp_val)[0][0], np.array(alt_interp))

    data = {}
    data['lats'] = np.array(lats)
    data['lons'] = np.array(lons)
    data['u_winds'] = np.array(u_winds)
    data['v_winds'] = np.array(v_winds)
    data['altitudes'] = np.array(altitudes)
    all_pressures = []

    for dat in data['lats']:
        all_pressures.append(100*np.array(pressures)) # Convert hPa to Pa

    data['pressures'] = np.array(all_pressures).transpose()

    return data

def calc_gefs_std(weather_file=None, current_time=None, loc0=None, descent_only=False):

    indir = '/home/ellen/Desktop/SuperBIT/Weather_data/GEFS/'

    tile_size = 10.
    lat0, lon0, alt0 = loc0
    area = lat0 + tile_size/2, lon0 - tile_size/2, lat0 - tile_size/2, lon0 + tile_size/2
    tlat, llon, blat, rlon = area

    datestr = weather_file[6:14]
    hhhh, hhh = int(int(weather_file[15:19])/100), int(weather_file[20:23])
    hour = hhhh + hhh
    closest_hhhh, closest_hhh1, closest_hhh2 = get_gfs.get_closest_hr(utc_hour=hour)

    files = [filename for filename in os.listdir(indir)]

    count = 0

    while count < 2:

        files = [filename for filename in os.listdir(indir)]

        if 'gespr_4_' + datestr + '_' + str(hhhh*100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2' in files \
        and 'geavg_4_' + datestr + '_' + str(hhhh*100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2' in files:

            data1 = read_gefs_file(fname='gespr_4_' + datestr + '_' + str(hhhh*100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2', area=area, alt0=alt0, descent_only=descent_only)
            data2 = read_gefs_file(fname='geavg_4_' + datestr + '_' + str(hhhh*100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2', area=area, alt0=alt0, descent_only=descent_only)

            data = {}
            data['lats'] = data1['lats']
            data['lons'] = data1['lons']
            data['u_winds'] = data1['u_winds']
            data['v_winds'] = data1['v_winds']
            data['altitudes'] = data2['altitudes']
            data['pressures'] = data1['pressures']

            data = pyb_aux.data_interpolation(data=data, alt0=alt0, step=100, descent_only=descent_only, output_figs=False)

            return data

        elif 'gens-a_3_' + datestr + '_' + str(hhhh*100).zfill(4) + '_' + str(closest_hhh2).zfill(3) + '_00.grb2' in files:

            u_winds, v_winds, altitudes = {}, {}, {}

            for no in range(0, 21):

                ens_fname = 'gens-a_3_' + datestr + '_' + str(hhhh*100).zfill(4) + '_' + str(closest_hhh2).zfill(3) + '_' + str(no).zfill(2) + '.grb2'
                data = read_gefs_file(fname=ens_fname, area=area, alt0=alt0, descent_only=descent_only)
                u_winds[no], v_winds[no], altitudes[no] = data['u_winds'], data['v_winds'], data['altitudes']

            u_winds_std, v_winds_std, altitudes_mean = {}, {}, {}

            lats, lons, pressures = data['lats'], data['lons'], data['pressures']

            for p in range(len((pressures))):
                key = pressures[p][0]
                u_winds_std[key] = []
                v_winds_std[key] = []
                altitudes_mean[key] = []
                for j in range(len(u_winds[0][p])):
                    u_winds_std[key].append(np.std([u_winds[i][p][j] for i in range(0, 21)]))
                    v_winds_std[key].append(np.std([v_winds[i][p][j] for i in range(0, 21)]))
                    altitudes_mean[key].append(np.mean([altitudes[i][p][j] for i in range(0, 21)]))

            u_winds, v_winds, altitudes = [], [], []

            for key in pressures[:, 0]:
                uwnd, vwnd, alt = [], [], []
                uwnd.append(u_winds_std[key])
                vwnd.append(v_winds_std[key])
                alt.append(altitudes_mean[key])

                u_winds.append(np.hstack(uwnd))
                v_winds.append(np.hstack(vwnd))
                altitudes.append(np.hstack(alt))

            data = {}
            data['lats'] = lats
            data['lons'] = lons
            data['u_winds'] = np.array(u_winds)
            data['v_winds'] = np.array(v_winds)
            data['altitudes'] = np.array(altitudes)
            data['pressures'] = pressures

            data = pyb_aux.data_interpolation(data=data, alt0=alt0, step=100, descent_only=descent_only, output_figs=False)

            return data

        else:

            get_gfs.get_gefs_files(datestr=datestr, utc_hour=current_time)

        count += 1

def add_uv_errs(main_data=None, err_data=None):

    print('Adding u/v wind errors...\n')

    main_data['u_wind_errs'] = err_data['u_winds']
    main_data['v_wind_errs'] = err_data['u_winds']
    main_data['lats_err'] = err_data['lats']
    main_data['lons_err'] = err_data['lons']

    return main_data

if __name__ == '__main__':

    loc0 = 30.776012, -5.613083, 25892

    tile_size = 10.
    lat0, lon0, alt0 = loc0
    area = lat0 + tile_size/2, lon0 - tile_size/2, lat0 - tile_size/2, lon0 + tile_size/2
    tlat, llon, blat, rlon = area

    # err_data = read_gefs_file(fname='gespr_4_20190624_1200_003.grb2', area=area, alt0=alt0, descent_only=True, step=100)

    data = calc_gefs_std(weather_file='gfs_4_20190514_1200_003.grb2', loc0=loc0,  current_time=15, descent_only=True)
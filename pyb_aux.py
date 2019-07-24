"""Auxiliary functions used in pyBalloon"""

from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.io import ascii
from shutil import copyfile
import datetime as dt
import pygrib as pg
import numpy as np
import gdal
import json
import sys
import os
import re

import get_gfs
import pyb_io

import param_file as p

# necessary constants
g_0 = 9.80665 # Earth gravitational acceleration at surface
R_e = 6371009 # mean Earth radius in meters
R = 8.3144621 # Ideal gas constant
M_air = 0.0289644 # molar mass of air [kg/mol], altitude dependence
M_helium = 4.002602
Cd_sphere = 0.47 # Drag coefficient for a sphere

def all_and(data):
    """Logical and for a list of arrays.
    
    Required arguments:
        - data -- list of Numpy boolean arrays

    Return:
        - result -- Logical and of all given arrays
    """
    result = data.pop()
    for d in data:
        result = np.logical_and(result, d)

    return result

# method to calculate Earth WGS84 based radius on a given latitude.
# see http://en.wikipedia.org/wiki/Earth_radius#Radius_at_a_given_geodetic_latitude
def earth_radius(lat_rad):
    
    # WGS84 reference ellipsoid
    a = 6378.137     # Earth equatorial radius, km
    b = 6356.7523142 # Earth polar radius, km
    
    cos_lat = np.cos(lat_rad)
    sin_lat = np.sin(lat_rad)

    r = np.sqrt(((a*a*cos_lat)**2 + (b*b*sin_lat)**2) / ((a*cos_lat)**2 + (b*sin_lat)**2))

    return r

# method to calculate the air density
def air_density(data):

    p = data['pressures'] # Pa
    T = data['temperatures'] # K

    if p.shape != T.shape:
        x, y = p.shape
        rho = [np.array((p[:, i] * M_air)/(R * T)) for i in range(0, y)]
        rho = np.array(rho).transpose()
    else:
        rho = (p * M_air)/(R * T)

    return rho # kg m-3

# method to interpolate data to altitude steps
def data_interpolation(data, alt0, step, mode='spline', descent_only=False, output_figs=True):

    altitudes = data['altitudes']

    new_data = {}

    if descent_only:
        # new_data['altitudes'] = np.arange(alt0 % step,  altitudes.max(), step)
        new_data['altitudes'] = np.arange(alt0 % step,  alt0 + step, step)
    else:
        new_data['altitudes'] = np.arange(alt0, altitudes.max(), step)

    new_data['lats'] = data['lats']
    new_data['lons'] = data['lons']

    checks = ['altitudes', 'lats', 'lons']

    for key in data.keys():        
        if key not in checks:

            arr = []
            d = data[key]

            try:
                x, y = d.shape

            except ValueError:

                x = 0
                y, = d.shape

            if mode == 'spline':
                if not descent_only:
                    if x > 0:
                        for i in range(0, y):
                            ok_idxs = altitudes[:, i] >= alt0
                            tck = interpolate.splrep(altitudes[ok_idxs, i], d[ok_idxs, i])
                            arr.append(np.array(interpolate.splev(new_data['altitudes'], tck)))
                    else:
                        tck = interpolate.splrep(altitudes, d)
                        arr.append(np.array(interpolate.splev(new_data['altitudes'], tck)))

                elif descent_only:                       
                    if x > 0:
                        for i in range(0, y):
                            ok_idxs = altitudes[:, i] <= alt0
                            # ok_idxs = altitudes[:, i] <= alt0 + 500
                            # ok_idxs = altitudes[:, i] <= 60000
                            tck = interpolate.splrep(altitudes[ok_idxs, i], d[ok_idxs, i])
                            arr.append(np.array(interpolate.splev(new_data['altitudes'], tck)))
                    else:
                        tck = interpolate.splrep(altitudes, d)
                        arr.append(np.array(interpolate.splev(new_data['altitudes'], tck)))

            else: # use linear interpolation 
                # There's something wrong here:
                for i in range(0, y):
                    for i in range(0, len(d)):
                        tck = interpolate.interp1d(altitudes[:, i], d[:, i])
                        arr.append(tck(new_data['altitudes']))
            new_data[key] = np.array(arr)

    if output_figs:

        figs = {}

        ind1 = np.where(altitudes[:, 0] == alt0)[0][0]
        ind2 = np.where(new_data['altitudes'] == alt0)[0][0]

        for key in data.keys():
            if key not in checks:

                fig = plt.figure()

                plt.axvline(alt0, linewidth=1, linestyle='--', label='Initial alt.')
                plt.plot(altitudes[:, 0], data[key][:, 0], 'ro--', label='Before interp.', markersize=3.5)
                plt.plot(altitudes[:, 0][ind1], data[key][:, 0][ind1], 'go', markersize=5, label='Inserted at alt. 0')
                plt.plot(new_data['altitudes'][:ind2+1], new_data[key][0][:ind2+1], 'bo', markersize=0.5, label='After interp.')
                plt.ylabel(key.capitalize().replace('_', ' '), fontsize=15)
                plt.xlabel('Altitude [m]', fontsize=15)
                plt.legend(loc='best')
                plt.grid(True)
                plt.tight_layout()

                figs[key] = fig

        return new_data, figs

    else:

        return new_data

def lift(data, mass):
    """Calculate effective lift (force, in Newtons) caused by the balloon.

    Required arguments:
        - data -- Dictionary containing 'altitudes' (meters), balloon
        'volumes' (m^3) and 'air_densities' (kg/m^3)
        - mass -- Mass of the whole balloon

    Return resultant lift force.
    """

    h = data['altitudes']
    V_b = data['balloon_volumes']
    rho_air = data['air_densities']
    g = g_0 * (R_e / (R_e + h))**2 # Gravitational acceleration at height h

    F_lift = g * (V_b*rho_air - mass)

    return F_lift


def balloon_volume(data):
    """Calculate volume of a sphere.

    Required argument:
        - data -- Dictionary containing 'balloon_radii' (in meters)

    Return:
        - Volume of the balloon at each level in data
    """

    r = data['balloon_radii']
    V = 4/3. * np.pi * r**3

    return V


def balloon_volume_ideal_gas(data, gas_mass, gas_molar_mass=M_helium):
    """Calculate gas (~balloon) volume based on ideal gas law: pV = nRT.

    Required arguments:
        - data -- Dictionary returned by read_gfs_data()
        - gas_mass -- Mass of the gas in question

    Optional arguments:
        - gas_molar_mass -- Gas molar mask (g/mol). Default: 4.002602 (Helium)

    Return:
        - Gas volume at each level in data
    """

    m = gas_mass # 
    M = gas_molar_mass/1000. # default to Helium, convert to kg/mol

    T = data['temperatures']
    p = data['pressures'] # pressure in Pascals

    # gas volume without the balloon
    V = m*R*T/(M*p) # pV = nRT, n = m/M

    return V


def burst_altitude(data, burst_radius):

    radii = data['balloon_radii']

    alt_out = []
    i_out = []
    i = 0
    for radius in radii:
        diff = np.abs(radius-burst_radius)
        idx = np.where(diff == diff.min())    
        i_out.append(idx[0][0])
        alt_out.append(data['altitudes'][idx])
        i+=1

    return np.array(alt_out), np.array(i_out)


def neutral_buoyancy_level(data):
    """Find the level of neutral buoyancy (or more precise, the level where effective lift is closest to zero).

    Required arguments:
        - data -- Dictionary containing 'altitudes' and corresponding 'lifts'

    Return:
        - Altitude of neutral buoyancy and corresponding array index
    """

    alt_out = []
    i_out = []
    idx = 0
    for lft in data['lifts']:
        lft = np.abs(lft)
        idx2 = np.where(lft == lft.min())
        i_out.append(idx2[0][0])
        alt_out.append(data['altitudes'][idx, idx2])
        idx += 1

    return np.array(alt_out), np.array(i_out)


def ascent_speed(data, mass, Cd=Cd_sphere):
    """Calculate the rate of ascent (in m/s) for the inflated balloon at given levels.

    Required arguments:
        - data -- Dictionary with 'altitudes' and corresponding
        'air_densities', 'balloon_radii' and 'balloon_volumes'
        - mass -- Full balloon mass (in kg)

    Optional arguments:
        - Cd -- Coefficient of drag. Default: 0.47 (sphere)

    Return:
        - Ascent speed (m/s) at every level in input data.
    """

    m = mass
    rho = data['air_densities']
    A = np.pi*data['balloon_radii']**2
    V = data['balloon_volumes']
    h = data['altitudes']

    g = g_0 * (R_e / (R_e + h))**2

    Fb = V*rho*g # byoyance
    Fg = m*g     # gravity
    F = Fb-Fg

    # Set negative buyoyancies to zero (we won't get there)
    idxs = np.where(F <= 0)
    if len(idxs[0]) > 0:
        idx = idxs[0][0]
        F[idx:] = 1e-30

    v = np.sqrt(2*F/(rho*Cd*A))

    return v

def descent_speed(data, mass, Cd, areas, alt_step, change_alt=None):
    """Calculate the rate of descent for deflated (burst) balloon with
    1 or 2 different sized parachutes with given areas, change
    altitude and drag-coefficient.

    Required arguments:
        - data -- Dictionary with 'altitudes', and corresponding 'air_densities'
        - mass -- Mass of the payload + assumed remnants of the balloon
        - Cd -- Coefficients of drag for one or two parachutes in a tuple
        - areas -- Effective areas (in m^2) of one or two parachutes in a tuple

    Optional arguments:
        - change_alt -- Altitude where first parachute is changed to
        the second one. If None, only one parachute is used. Default:
        None

    Return:
        - Rate of descent (m/s) for each level in input data
    """

    m = mass
    h = data['altitudes']
    g = g_0 * (R_e / (R_e + h))**2 # Gravitational acceleration at height h

    speeds = []
    if change_alt is not None:
        idxs = h < change_alt
        for rho in data['air_densities']:
            v = np.sqrt(2*m*g/(rho*Cd*areas[0]))
            v[idxs] = np.sqrt(2*m*g[idxs]/(rho[idxs]*Cd*areas[1]))
            speeds.append(v)
    else:
        factor = 1
        for rho in data['air_densities']:
            v = np.sqrt(2*m*g/(rho*Cd*areas))
            speeds.append(v)
    return -1*np.array(speeds)
    

def mooney_rivlin(data, radius_empty, radius_filled, thickness_empty, gas_molar_mass=M_helium):

    r0 = radius_empty # in meters
    r1 = radius_filled # in meters
    t0 = thickness_empty # in meters
    M = gas_molar_mass / 1000. # convert to kg/mol
    p = data['pressures']
    p0 = p[0]
    T = data['temperatures']
    T0 = T[0]

    mu = 300000. # balloon shear modulus in Pascals
    alfa = 10./11.

    # Amount of gas in moles
    n = (4/3. * np.pi * r1**3)/(R*T0) * (p0 - 2*mu*(t0/r0)*((r0/r1) - (r0/r1)**7) * (1 + (1/alfa - 1) * (r1/r0)**2))
    gas_mass = n*M

    # Solve balloon-radius-roots for each height level
    
    # Constants for the polynomials
    a8 = (1/alfa-1)*2*mu*t0/(r0**2)
    a6 = 2*mu*t0
    a2 = -2*mu*t0*(1/alfa-1)*(r0**4)
    a0 = -2*mu*t0*(r0**6)
    
    all_radii = []

    try:
        x, y = p.shape
        
        for i in range(0, x):
            radii = []
            for j in range(0, y):
                a4 = -3*n[j]*R/(4*np.pi)
                # 8th degree polynomial
                poly = [a8,        # r^8
                        p[i,j],    # r^7
                        a6,        # r^6
                        0,         # r^5
                        a4*T[i,j], # r^4
                        0,         # r^3
                        a2,        # r^2
                        0,         # r^1
                        a0]        # r^0

                roots = np.roots(poly)
        
                for r in roots:
                    if r.real > 0 and r.imag == 0:
                        radii.append(r.real)
            all_radii.append(np.array(radii))


    except ValueError:        
        for i in range(0, len(p)):
            a4 = -3*n*R/(4*np.pi)
            # 8th degree polynomial
            poly = [a8,        # r^8
                    p[i],    # r^7
                    a6,        # r^6
                    0,         # r^5
                    a4*T[i], # r^4
                    0,         # r^3
                    a2,        # r^2
                    0,         # r^1
                    a0]        # r^0

            roots = np.roots(poly)
        
            for r in roots:
                if r.real > 0 and r.imag == 0:
                    all_radii.append(r.real)

    all_radii = np.array(all_radii)

    return all_radii, gas_mass

def get_elevation(lon, lat):

    SAMPLES = 1201  # Change this to 3601 for SRTM1
    HGTDIR = '/home/ellen/Desktop/SuperBIT/SRTM_data/'

    if lat > 180:
        lat -= 360
    if lon > 180:
        lon -= 360

    if lat >= 0:
        ns = 'N'
    elif lat < 0:
        ns = 'S'

    if lon >= 0:
        ew = 'E'
    elif lon < 0:
        ew = 'W'

    hgt_file = "%(ns)s%(lat)02d%(ew)s%(lon)03d.hgt" % {'lat': abs(lat), 'lon': abs(lon), 'ns': ns, 'ew': ew}
    hgt_file_path = os.path.join(HGTDIR, hgt_file)

    if not os.path.isfile(hgt_file_path):
        print(str(hgt_file) + ' does not exist!')
        return -32768

    else:
        with open(hgt_file_path, 'rb') as hgt_data:
            # HGT is 16bit signed integer(i2) - big endian(>)
            elevations = np.fromfile(hgt_data, np.dtype('>i2'), SAMPLES*SAMPLES).reshape((SAMPLES, SAMPLES))

            lat_row = int(round((lat - int(lat)) * (SAMPLES - 1), 0))
            lon_row = int(round((lon - int(lon)) * (SAMPLES - 1), 0))

            return elevations[SAMPLES - 1 - lat_row, lon_row].astype(int)

# method to get more accurate endpoint for predictions as they can go underground.
def get_endpoint(data=None, filename=None, params=None):

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
        fext = 'descent_only/start_point' + next_point + '/'
    else:
        fext = 'ascent+descent/'

    dir_pred = '/home/ellen/Desktop/SuperBIT/Weather_data/Trajectories/' + fext

    if data == None:

        data  = ascii.read(dir_pred + filename)

        lats = data['lats']
        lons = data['lons']
        alts = data['alts']
        dists = data['dists']

    else:

        lats, lons, alts, dists = data

    elevations = []

    print('Checking altitudes & elevations...')

    for i in range(1, len(lats)):
        elevation = get_elevation(lat=lats[-i], lon=lons[-i])
        elevations.append(elevation)
        if elevation < alts[-i]:
            break

    if i == 1:
        return (lats[-i], lons[-i]), elevation

    inds = (-i, -i + 1)

    dlon = lons[inds[0]] - lons[inds[1]]
    dlat = lats[inds[0]] - lats[inds[1]]

    x1 = dists[inds[1]]

    y1 = alts[inds[0]] - alts[inds[1]]
    y2 = elevations[-2] - alts[inds[1]]

    dx = x1*(y2/y1)
    x2 = x1 - dx
    f = dx/x1

    newlon = lons[inds[1]] + f*dlon
    newlat = lats[inds[1]] + f*dlat

    return (newlat, newlon), get_elevation(lat=newlat, lon=newlon)

def calc_gefs_errs(weather_file=None, current_time=None, loc0=None, descent_only=False):

    indir = '/home/ellen/Desktop/SuperBIT/Weather_data/GEFS/'

    tile_size = 10.0
    lat0, lon0, alt0 = loc0
    area = (lat0 + tile_size / 2, lon0 - tile_size / 2, lat0 - tile_size / 2, lon0 + tile_size / 2)
    tlat, llon, blat, rlon = area

    datestr = weather_file[6:14]
    hhhh, hhh = int(int(weather_file[15:19]) / 100), int(weather_file[20:23])
    hour = hhhh + hhh

    hrs = get_gfs.get_closest_hr(utc_hour=hour)
    closest_hhhh, closest_hhh1, closest_hhh2 = hrs[0], hrs[1], hrs[2]

    files = [filename for filename in os.listdir(indir)]

    count = 0
    while count < 2:

        files = [filename for filename in os.listdir(indir)]

        if 'gespr_4_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2' in files and 'geavg_4_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2' in files:

            data1 = pyb_io.read_gefs_file(fname='gespr_4_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2', area=area, alt0=alt0, descent_only=descent_only)
            data2 = pyb_io.read_gefs_file(fname='geavg_4_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(hhh).zfill(3) + '.grb2', area=area, alt0=alt0, descent_only=descent_only)
            
            data = {}
            data['lats'] = data1['lats']
            data['lons'] = data1['lons']
            data['u_winds'] = data1['u_winds']
            data['v_winds'] = data1['v_winds']
            data['altitudes'] = data2['altitudes']
            data['pressures'] = data1['pressures']

            data = data_interpolation(data=data, alt0=alt0, step=100, descent_only=descent_only, output_figs=False)

            return data

        if 'gens-a_3_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(closest_hhh2).zfill(3) + '_00.grb2' in files:

            u_winds, v_winds, altitudes = {}, {}, {}

            print('Reading 20 GEFS files...')

            for no in range(0, 20):

                ens_fname = 'gens-a_3_' + datestr + '_' + str(hhhh * 100).zfill(4) + '_' + str(closest_hhh2).zfill(3) + '_' + str(no + 1).zfill(2) + '.grb2'
                data = pyb_io.read_gefs_file(fname=ens_fname, area=area, alt0=alt0, descent_only=descent_only)
                u_winds[no], v_winds[no], altitudes[no] = data['u_winds'], data['v_winds'], data['altitudes']

            u_winds_std, v_winds_std, altitudes_mean = {}, {}, {}
            lats, lons, pressures = data['lats'], data['lons'], data['pressures']

            for p in range(len(pressures)):

                key = pressures[p][0]
                u_winds_std[key] = []
                v_winds_std[key] = []
                altitudes_mean[key] = []

                for j in range(len(u_winds[0][p])):

                    u_winds_std[key].append(np.std([u_winds[i][p][j] for i in range(0, 20)]))
                    v_winds_std[key].append(np.std([v_winds[i][p][j] for i in range(0, 20)]))
                    altitudes_mean[key].append(np.mean([altitudes[i][p][j] for i in range(0, 20)]))

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

            data = data_interpolation(data=data, alt0=alt0, step=100, descent_only=descent_only, output_figs=False)

            return data

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

    fname1 = 'gfs_4_20181028_0600_006.grb2'
    fname2 = 'gfs_4_20190530_0000_000.grb2'
    descent_only = True
    loc0 = 30.776012, -5.613083, 25892

    new_data = calc_uv_errs(weather_file=fname1, loc0=loc0, descent_only=descent_only, current_time=6)
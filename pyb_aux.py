from scipy import interpolate
import numpy as np
import re

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

    r = np.sqrt(((a*a*cos_lat)**2 + (b*b*sin_lat)**2) / 
                ((a*cos_lat)**2 + (b*sin_lat)**2))

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
def data_interpolation(data, alt0, step, mode='spline', descent_only=False):

    altitudes = data['altitudes']

    new_data = {}

    if descent_only:
        new_data['altitudes'] = np.arange(alt0 % step,  altitudes.max(), step)
        # new_data['altitudes'] = np.arange(alt0 % step,  alt0 + step, step)
    else:
        new_data['altitudes'] = np.arange(alt0, altitudes.max(), step)

    new_data['lats'] = data['lats']
    new_data['lons'] = data['lons']

    checks = ['altitudes', 'lats', 'lons']

    # print(data.keys())

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

    return new_data

def lift(data, mass):
    """Calculate effective lift (force, in Newtons) caused by the
    balloon.

    Required arguments:
        - data -- Dictionary containing 'altitudes' (meters), balloon
        'volumes' (m^3) and 'air_densities' (kg/m^3)
        - mass -- Mass of the whole balloon

    Return:
        - Resultant lift force.
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
    """Calculate gas (~balloon) volume based on ideal gas law: pV =
    nRT.

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
    """Find the level of neutral buoyancy (or more precise, the level
    where effective lift is closest to zero).

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
    """Calculate the rate of ascent (in m/s) for the inflated balloon
    at given levels.

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
    

def mooney_rivlin(data, radius_empty, radius_filled, 
                  thickness_empty, gas_molar_mass=M_helium):

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

# method to find the elevation at the location (the loc closest in the grid)
def get_elevation(lat, lon, srtm_file):

    srtm_dir = '/home/ellen/Desktop/SuperBIT/SRTM_data/'

    if lon < 0:
        lon = lon + 360 % 360

    if srtm_file == None:
        data = ascii.read(srtm_dir + 'srtm_data_limits.txt')
        files, tlat, blat, llon, rlon = data['file'], data['latN'], data['latS'], data['lonW'], data['lonE']

        file = None

        for i in range(len(files)):

            if lat < tlat[i] and lat > blat[i]:
                if lon < rlon[i] and lon > llon[i]:
                    file = files[i]

        if file != None:
            srtm_file = str(file)
        else:
            print('Correct SRTM file data is not here!')
            return

    file = srtm_dir + srtm_file

    ds = gdal.Open(file)

    band = ds.GetRasterBand(1)
    elevations = band.ReadAsArray()
    nrows, ncols = elevations.shape
    x0, dx, dxdy, y0, dydx, dy = ds.GetGeoTransform()

    lons = np.arange(x0, x0 + dx*(ncols-1) + dx, dx) 
    lon_arr = np.array([lons for i in range(nrows)])
    lats = np.arange(y0, y0 + dy *nrows, dy)                                                   
    lat_arr = np.array([np.array([l for i in range(ncols)]) for l in lats])

    diff = np.sqrt((lon_arr - lon)**2 + (lat_arr - lat)**2)
    min_diff = np.min(diff)
    i0, i1 = np.where(diff == min_diff)

    elevation = elevations[int(i0)][int(i1)]

    elevations, diff, lons, lats = None, None, None, None

    return elevation

# method to match prediction to gps file
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

# method to match gps to prediction file
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
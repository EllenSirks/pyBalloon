#!/usr/bin/python

import numpy as np
import sys
import time
import pyb_io
import pyb_traj
import get_gfs
import os

##import matplotlib as mpl
##from mpl_toolkits.mplot3d import Axes3D
##import matplotlib.pyplot as plt
##from mpl_toolkits.basemap import Basemap

if __name__ == "__main__":

    in_dir = '/home/ellen/Desktop/SuperBIT/Weather_data/'

    datestr = sys.argv[1]
    utc_hour = sys.argv[2]

    file_name = in_dir + 'gfs_4_' + datestr +'_*'

    # if os.path.isfile(file_name):
    #     print('yeah')
    #     get_gfs.get_gfs_data(datestr, utc_hour)

    # yyyy, mm, dd = int(datestr[0:4]), int(datestr[4:6]), int(datestr[6:8])

    time0 = time.time()

    # example balloon:
    balloon = {}
    balloon['altitude_step'] = 50.0 # meters
    balloon['equip_mass'] = 1.554 # kg
    balloon['balloon_mass'] = 1.50 # kg
    balloon['fill_radius'] = 2.122/2 # meters
    balloon['radius_empty'] = 1.61/2 # meters
    balloon['burst_radius'] = 9.44/2 # meters
    balloon['thickness_empty'] = 0.2 * 10**-3 # mm -> meters
    balloon['Cd_balloon'] = 0.35
    balloon['Cd_parachute'] = 0.8
    balloon['parachute_areas'] = np.pi * np.array([1.07])**2 # m^2
    balloon['parachute_change_altitude'] = None # meters

    lat0 = 55.0 # degrees
    lon0 = -1.5 # degrees
    alt0 = 100. # meters

    # note top, left, bottom, right ordering for area
    model_data = pyb_io.read_gfs_single(file_name, 
                                    (lat0+1.5, lon0-1.5, 
                                     lat0-1.5, lon0+1.5))
                                    #     ,
                                    # ens_main=None,
                                    # ens_member_pattern=None, alt0=alt0)

    print( 'GFS data read, %.1f s elapsed' % (time.time() - time0))

    loc0 = (lat0, lon0, alt0)

    trajectories = []

    for data in model_data:
        trajectories.append(pyb_traj.calc_movements(data, loc0, balloon))

    print( 'Trajectories calculated, %.1f s elapsed' % (time.time() - time0))

    # highest point in main-run trajectory
    idx, = np.where(trajectories[0]['alts'] == np.max(trajectories[0]['alts']))
    latx, _ = trajectories[0]['lats'][idx]
    lonx, _ = trajectories[0]['lons'][idx]
    altx, _ = trajectories[0]['alts'][idx]
    timex, _ = trajectories[0]['times'][idx]
    print( latx, lonx, altx, '%.0f minutes' % (timex))
    other_info = [(latx, lonx, altx, 'Burst point', '%.0f minutes, %.0f meters' % (timex, altx))]

    kml_fname = 'pyballoon_trajectories.kml'
    pyb_io.save_kml(kml_fname, trajectories, other_info=other_info)
    
    print( 'Program finished in %.1f s' % (time.time() - time0))


##    fig = plt.figure()
##    m = Basemap(projection='stere', lon_0=lon0, lat_0=lat0, lat_ts=lat0, 
##                llcrnrlat=lat0-2, urcrnrlat=lat0+2, 
##                llcrnrlon=lon0-2, urcrnrlon=lon0+2,
##                rsphere=6371200.,resolution='l',area_thresh=10000)
##    m.drawcoastlines()
##    m.drawcountries()
##    i = 0
##    trajectories.reverse()
##    for t in trajectories:
##        x, y = m(t['lons'], t['lats'])
##        if i == len(trajectories)-1:
##            plt.plot(x, y, 'r')
##        else:
##            plt.plot(x, y, 'k.')
##        i += 1
##    #plt.ylabel('Latitude (deg)')
##    #plt.xlabel('Longitude (deg)')
##    #plt.legend(['U','V','Total'])
##    plt.show()


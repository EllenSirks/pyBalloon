pyBalloon
=========

Python scripts that can be used to simulate trajectories of tropospheric weather balloons and parachutes.

Dependencies
------------
- Python 2.x (tested with 2.7.15+)
- Numpy
- Scipy
- astropy
- matplotlib
- requests
- pygrib
- shutil
- datetime
- haversine
- shapely
- functools
- pyproj

Folders
------------
- pyBalloon: contains the code. 
- Output: here the output from pyBalloon is saved (this folder is created automatically).
- Weather_data: contains the GEFS and GFS folders, which contain the ensemble weather forecasts and the weather forecasts models respectively.
- SRTM_data: contains the elevation data.

Notes
------------
- In param_file.py you can set the path to the folder containing pyBalloon, if you wish to change the folder names you will have to go through all the files.
- Elevation data needs to be downloaded from: http://viewfinderpanoramas.org/dem3.html#hgt
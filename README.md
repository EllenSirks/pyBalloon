pyBalloon
=========

Python scripts that can be used to simulate trajectories of tropospheric weather balloons and parachutes.

Dependencies
------------
- Tested with Python 2.7.15+ and Python 3.6.8
- Numpy
- Scipy
- astropy
- matplotlib
- requests
- pygrib
- shutil
- datetime
- shapely
- functools
- pyproj

Each function has documentation at the beginning of said function. These docstrings are also in html format in the docs/ subdirectory.

Folders
------------
- pyBalloon: contains the code. 
- Output: here the output from pyBalloon is saved.
- Weather_data: contains the GEFS and GFS folders, which contain the ensemble weather forecasts and the weather forecasts models respectively.
- SRTM_data: contains the elevation data.

Before running the program, please run the script setup.py.
This script while create all the necessary folders if you have changed the path and folder names in param_file.py

Notes
------------
- Elevation data needs to be downloaded from: http://viewfinderpanoramas.org/dem3.html#hgt
- Most weather GFS data should be downloaded automatically. However, sometimes a given date does not exist on the server and needs to be downloaded from somewhere else.
from shapely.geometry import Point, mapping
from shapely.ops import transform
from functools import partial
import pyproj
import json
import sys

def create_circle(lon=None, lat=None):

	point = Point(lon, lat)

	local_azimuthal_projection = '+proj=aeqd +R=6371000 +units=m +lat_0={point.y} +lon_0={point.x}'

	wgs84_to_aeqd = partial(pyproj.transform, pyproj.Proj('+proj=longlat +datum=WGS84 +no_defs'), pyproj.Proj(local_azimuthal_projection), )
	aeqd_to_wgs84 = partial(pyproj.transform, pyproj.Proj(local_azimuthal_projection), pyproj.Proj('+proj=longlat +datum=WGS84 +no_defs'), )

	point_transformed = transform(wgs84_to_aeqd, point)

	buffer = point_transformed.buffer(10000)

	buffer_wgs84 = transform(aeqd_to_wgs84, buffer)
	coords = json.dumps(mapping(buffer_wgs84))

	from io import StringIO
	io = StringIO(unicode(coords))
	coords = json.load(io)['coordinates'][0]

	kml_str = '<name>error circle</name>'
	kml_str += '<description>error circle</description>'
	kml_str += '<Style id="stylesel_362">'
	kml_str += '<LineStyle id="substyle_363">'
	kml_str += '<color>BF0000DF</color>'
	kml_str += '<colorMode>normal</colorMode>'
	kml_str += '<width>5</width>'
	kml_str += '</LineStyle>'
	kml_str += '<PolyStyle id="substyle_364">'
	kml_str += '<color>BF0000DF</color>'
	kml_str += '<colorMode>normal</colorMode>'
	kml_str += '<fill>1</fill>'
	kml_str += '<outline>1</outline>'
	kml_str += '</PolyStyle>'
	kml_str += '</Style>'
	kml_str += '<Placemark id="feat_91">'
	kml_str += '<styleUrl>#stylesel_362</styleUrl>'
	kml_str += '<LineString id="geom_86">'
	kml_str += '<coordinates>'

	for i in range(len(coords)):
		kml_str += str(coords[i][0]) + ',' + str(coords[i][1]) + ',' +  '0' +'\n'

	kml_str += '</coordinates>\n'
	kml_str += '<extrude>1</extrude>\n'
	kml_str += '<tessellate>1</tessellate>\n'
	kml_str += '<altitudeMode>relativeToGround</altitudeMode>\n'
	kml_str += '</LineString>\n'
	kml_str += '</Placemark>\n'

	return kml_str
	
if __name__ == '__main__':

	kml_str = create_circle(float(sys.argv[1]), float(sys.argv[2]))
	print(kml_str)
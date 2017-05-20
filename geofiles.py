import gdal
import gdalconst
import pdb

class gdalFile(object):
    """
    This is the class used to open a georeferenced raster file and make
    gridded data available. Raster values for a coordinate are extracted
    using the extract_coord_val method. This works for .bil and GeoTiff (.tif) 
    files. 
    """
    def __init__(self, geofile):
        """
        Construct the gdalFile object and initialize its
        properties.

        Args:
            geofile: path to file. May include "/vsizip/"
        """
        self.fname = geofile
        # Open the file and put in dataset
        dataset = gdal.Open(self.fname, gdalconst.GA_ReadOnly)
        self.nbands = dataset.RasterCount
        self.ncol = dataset.RasterXSize
        self.nrow = dataset.RasterYSize
        print('Opened {0}\n {1} raster bands\n{2} columns, {3} rows'.format(
            fname, str(self.nbands), str(self.ncol), str(self.nrow)))

        # Get georeferencing information and assign to object
        # properties
        geotransform = dataset.GetGeoTransform()
        self.originX = geotransform[0]
        self.originY = geotransform[3]
        self.pixelWidth = geotransform[1]
        self.pixelHeight = geotransform[5]
        def getbands(ds, nbands):
            return [ds.GetRasterBand(x+1).ReadAsArray() for x in range(nbands)]
        self.data = getbands(dataset, self.nbands)

    def extract_coord_val(self, lat, lon, bandnum=0):
        """
        Method for extracting raster values at given coordinates

        Args:
            lat: latitude, float
            lon: longitude, float
            bandnum: if raster has multiple bands identify band to extract from
        Returns:
            numeric: pixel value at given lat/lon
        """
        y = int((lat - self.originY)/self.pixelHeight)
        x = int((lon - self.originX)/self.pixelWidth)
        return self.data[bandnum][y, x]

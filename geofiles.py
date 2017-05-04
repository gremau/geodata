import gdal
import gdalconst

class BilFile(object):
# This is the class used to open a .bil file and make georeferenced data
# available. Raster values for a coordinate are extracted using the 
# extract_coord_val method

    def __init__(self, bil_file):
        # Construct the BilFile object and initialize its
        # properties
        self.bil_file = bil_file
        self.hdr_file = bil_file.split('.')[0]+'.hdr'
        # Register EHdr driver needed for these files - automatic in python?
        # gdal.GetDriverByName('EHdr').Register()
        # Open the file and put in dataset
        dataset = gdal.Open(self.bil_file, gdalconst.GA_ReadOnly)
        band = dataset.GetRasterBand(1)
        self.ncol = dataset.RasterXSize
        self.nrow = dataset.RasterYSize
        # Get georeferencing information and assign to object
        # properties
        geotransform = dataset.GetGeoTransform()
        self.originX = geotransform[0]
        self.originY = geotransform[3]
        self.pixelWidth = geotransform[1]
        self.pixelHeight = geotransform[5]
        self.data = band.ReadAsArray()

    def extract_coord_val(self, lat, lon):
        # Method for extracting raster values at given coordinates
        y = int((lat - self.originY)/self.pixelHeight)
        x = int((lon - self.originX)/self.pixelWidth)
        return self.data[y, x]

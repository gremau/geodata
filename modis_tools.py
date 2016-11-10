"""
This was adapted from the MOD17 python example found here:

http://hdfeos.org/zoo/index_openLPDAAC_Examples.php#MOD
"""

import os
import re
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import pyproj
import numpy as np
from pyhdf.SD import SD, SDC
import pdb as pdb

def list_datasets(FILE_NAME):
    # Open file.
    hdf = SD(FILE_NAME, SDC.READ)
    # List available SDS datasets.
    print(hdf.datasets())

def quickplot_2d_dataset(FILE_NAME, DATASET_NAME):
    hdf = SD(FILE_NAME, SDC.READ)
    data2D = hdf.select(DATASET_NAME)
    data = data2D[:,:].astype(np.double)
    # Make a figure
    fig, ax = plt.subplots(figsize=(6,6))
    ax.imshow(data[:,:], cmap=plt.cm.Greys, vmin=0, vmax=6000)
    plt.show()

class Mod17File(object):
# This is the class used to open a .hdr file and make georeferenced data
# available. Raster values for a coordinate are extracted using the 
# extract_coord_val method

    def __init__(self, mod_file, grid_name, data_field):
        import gdal
        import gdalconst
        # Construct dataset name
        dname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(mod_file,
                grid_name, data_field)

        # Construct the Mod17File object and initialize its
        # properties
        self.mod_file = mod_file
        self.grid_name = grid_name
        self.data_field = data_field
        self.dataset_name = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(
                self.mod_file, self.grid_name, self.data_field)
        # Open the file and put in dataset
        dataset = gdal.Open(self.dataset_name, gdalconst.GA_ReadOnly)
        #band = dataset.GetRasterBand(1)
        self.ncol = dataset.RasterXSize
        self.nrow = dataset.RasterYSize
        # Get georeferencing information and assign to object
        # properties
        geotransform = dataset.GetGeoTransform()
        self.originX = geotransform[0]
        self.originY = geotransform[3]
        self.pixelWidth = geotransform[1]
        self.pixelHeight = geotransform[5]

        # Read othher attributes.
        meta = dataset.GetMetadata()
        long_name = meta['long_name']        
        units = meta['units']
        _FillValue = np.float(meta['_FillValue'])
        scale_factor = np.float(meta['scale_factor'])
        add_offset = np.float(meta['add_offset'])
        valid_range = [np.float(x) for x in meta['valid_range'].split(', ')] 

        # Construct the grid (sinusoidal projection)
        x = np.linspace(self.originX, self.originX + self.pixelWidth * 
                self.ncol, self.ncol)
        y = np.linspace(self.originY, self.originY + self.pixelHeight * 
                self.nrow, self.nrow)
        self.xv, self.yv = np.meshgrid(x, y)


        #self.data = band.ReadAsArray()
        self.data = dataset.ReadAsArray().astype(np.float64)
    
    def reproj_wgs84(self):
        # In basemap, the sinusoidal projection is global, so we won't use it.
        # Instead we'll convert the grid back to lat/lons.
        sinu = pyproj.Proj("+proj=sinu +R=6371007.181000 +nadgrids=@null +wktext")
        wgs84 = pyproj.Proj("+init=EPSG:4326") 
        self.Lon, self.Lat= pyproj.transform(sinu, wgs84, self.xv, self.yv)


    def extract_coord_val(self, lon_sel, lat_sel):
        # Convert requested lat/lon to sinusoidal coords
        sinu = pyproj.Proj("+proj=sinu +R=6371007.181000 +nadgrids=@null +wktext")
        wgs84 = pyproj.Proj("+init=EPSG:4326") 
        lon_sel_sinu, lat_sel_sinu= pyproj.transform(wgs84, sinu, lon_sel, lat_sel)

        # Method for extracting raster values at given coordinates
        y_sel = int((lat_sel_sinu - self.originY)/self.pixelHeight)
        x_sel = int((lon_sel_sinu - self.originX)/self.pixelWidth)
        return self.data[y_sel, x_sel]

def get_Mod17File(mod_file, grid_name, data_field):

    return Mod17File(mod_file, grid_name, data_field)

# Function for extracting daily PRISM data
def getDailyPrism(year, metdata, data_path, coords_file):
    # Read in site coordinates, get date range and create a DataFrame
    #to fill
    pnts = pd.read_csv(coords_file)
    drange =  pd.date_range('1-1-{0}'.format(year),
            '12-31-{0}'.format(year), freq='D')
    df = pd.DataFrame(index=drange, columns=pnts.sitecode)
    for i in range(len(drange)):
        # Create a tuple to fill the file dates in,
        # pad month & day with zeros
        ymd_tuple = (str(drange.year[i]),
                str(drange.month[i]).zfill(2),
                str(drange.day[i]).zfill(2))
        # Create GDAL vsizip file path (read directly from zip archive)
        # See https://trac.osgeo.org/gdal/wiki/UserDocs/ReadInZip
        # If for some reason this vsizip interface won't work, extract the
        # archive and remove the 'vsizip' part of pathname
        bil_file = (r'/vsizip/' + data_path +
        r'PRISM_{0}_stable_4kmD2_{1}0101_{1}1231_bil.zip/'.format(
            metdata, year) +
        r'PRISM_{0}_stable_4kmD2_{1}{2}{3}_bil.bil'.format(metdata, *ymd_tuple))
        bil_ds = BilFile(bil_file)
        for j in range(len(pnts.index)):
            pt_val = bil_ds.extract_coord_val(pnts.lat[j], pnts.lon[j])
            df.iloc[i, j] = pt_val
    return df

    
def get_dataset_wgs84(FILE_NAME, GRID_NAME, DATAFIELD_NAME, plot=False):
    '''
    FILE_NAME = Name of the HDF-EOS file
    GRID_NAME, and DATAFIELD_NAME parameters can all be found either by
        examining the .hdf file in HDFView or running list_mod_datasets
    '''
    import gdal
    
    gname = 'HDF4_EOS:EOS_GRID:"{0}":{1}:{2}'.format(FILE_NAME,
            GRID_NAME, DATAFIELD_NAME)
    gdset = gdal.Open(gname)
    data = gdset.ReadAsArray().astype(np.float64)

    # Construct the grid.
    x0, xinc, _, y0, _, yinc = gdset.GetGeoTransform()
    nx, ny = (gdset.RasterXSize, gdset.RasterYSize)
    x = np.linspace(x0, x0 + xinc*nx, nx)
    y = np.linspace(y0, y0 + yinc*ny, ny)
    xv, yv = np.meshgrid(x, y)
    #pdb.set_trace()

    # In basemap, the sinusoidal projection is global, so we won't use it.
    # Instead we'll convert the grid back to lat/lons.
    sinu = pyproj.Proj("+proj=sinu +R=6371007.181000 +nadgrids=@null +wktext")
    wgs84 = pyproj.Proj("+init=EPSG:4326") 
    lon, lat= pyproj.transform(sinu, wgs84, xv, yv)

    # Read the attributes.
    meta = gdset.GetMetadata()
    long_name = meta['long_name']        
    units = meta['units']
    _FillValue = np.float(meta['_FillValue'])
    scale_factor = np.float(meta['scale_factor'])
    add_offset = np.float(meta['add_offset'])
    valid_range = [np.float(x) for x in meta['valid_range'].split(', ')] 
    # Close the hdf file
    del gdset
        
    invalid = np.logical_or(data > valid_range[1],
                            data < valid_range[0])
    invalid = np.logical_or(invalid, data == _FillValue)
    data[invalid] = np.nan
    data = scale_factor * (data - add_offset)
    data = np.ma.masked_array(data, np.isnan(data))

    if plot:
        # Make a simple plot
        m = Basemap(projection='cyl', resolution='h',
                llcrnrlat= 25.5, urcrnrlat = 42.5,
                llcrnrlon=-122.5, urcrnrlon = -95.5)
        #m = Basemap(projection='cyl', resolution='h',
        #        llcrnrlat= 32.45, urcrnrlat = 32.5,
        #        llcrnrlon=-105.5, urcrnrlon = -105.45)
        m.drawcoastlines(linewidth=0.5)
        m.drawparallels(np.arange(-10, 5, 5), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-70, -55, 5), labels=[0, 0, 0, 1])
        m.pcolormesh(lon, lat, data, latlon=True)
        cb = m.colorbar()
        cb.set_label(units)
        #m.scatter(lon, lat, c='y')


        basename = os.path.basename(FILE_NAME)
        plt.title('{0}\n{1}'.format(basename, long_name), fontsize=11)
        fig = plt.gcf()
        plt.show()
        #pngfile = "{0}.py.png".format(basename)
        #fig.savefig(pngfile)

    return [data, lat, lon]
    

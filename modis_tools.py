"""
This was adapted from the MOD17 python example found here:

http://hdfeos.org/zoo/index_openLPDAAC_Examples.php#MOD
"""

import os
import re


import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
import numpy as np

def list_datasets(FILE_NAME):
    from pyhdf.SD import SD, SDC

    # Open file.
    hdf = SD(FILE_NAME, SDC.READ)
    # List available SDS datasets.
    print(hdf.datasets())


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

    # In basemap, the sinusoidal projection is global, so we won't use it.
    # Instead we'll convert the grid back to lat/lons.
    sinu = pyproj.Proj("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext")
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
        m.drawcoastlines(linewidth=0.5)
        m.drawparallels(np.arange(-10, 5, 5), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-70, -55, 5), labels=[0, 0, 0, 1])
        m.pcolormesh(lon, lat, data, latlon=True)
        cb = m.colorbar()
        cb.set_label(units)

        basename = os.path.basename(FILE_NAME)
        plt.title('{0}\n{1}'.format(basename, long_name), fontsize=11)
        fig = plt.gcf()
        plt.show()
        #pngfile = "{0}.py.png".format(basename)
        #fig.savefig(pngfile)

    return [data, lat, lon]
    

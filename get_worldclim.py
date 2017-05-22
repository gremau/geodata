import pandas as pd
import geofiles as gf
from time import strptime
import pdb
import os


# Function for extracting 30 year normal PRISM precip data
def get_30s_avgmonthly(wc_path, lat, lon, datatype='prec'):
    """
    Extract 800m 30 year normal PRISM data for a list of lat/lon coordinates.

    Args:
        wc_path (str) : path to directory containing prism zip files
        lat (float array): latitude coordinates 
        lon (float array : longitude coordinates
        datamonth (str)  : Monthly normal to return (eg. 'Jan') or 'annual'
        datatype (str)   : abbreviation for the PRISM data type
                           (ppt, tmax, tmin)
    Returns:
        A pandas dataframe with lat/lon coordinates and PRISM normal data

    Note that this is looking for 800m "all" normal files, which are
    downloaded here:
    http://www.prism.oregonstate.edu/normals/
    
    by using the "Download All Normals Data(.bil)" option. The downloaded files
    contain both the annual and monthly normal datasets in separate .bil files.
    To get monthly files provide the proper 3 letter month code to datamonth.
    """
    mocols = ['{0:02}_{1}_{2}'.format(x+1, datatype, 'WCn') for x in range(12)]
    df = pd.DataFrame(index=lat.index, columns=['lat','lon'] + mocols)

    for i in range(12):
        month = "02".format(i+1)
        fname = r'/vsizip/' + os.path.join(wc_path,
            r'wc2.0_30s_{0}.zip'.format(datatype), 
            r'wc2.0_30s_{0}_{1}.tif'.format(datatype, month))
        pdb.set_trace()
        ds = gf.gdalFile(fname)
        for j in range(len(lat)):
            df.loc[j, ['lat', 'lon']] = [lat[j], lon[j]]
            df.loc[j, mocols[i]] = ds.extract_coord_val(
                lat[j], lon[j])
    return df


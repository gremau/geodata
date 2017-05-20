import pandas as pd
import geofiles as gf
from time import strptime
import pdb
import os

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
        bil_ds = gdalFile(bil_file)
        for j in range(len(pnts.index)):
            pt_val = bil_ds.extract_coord_val(pnts.lat[j], pnts.lon[j])
            df.iloc[i, j] = pt_val
    return df

# Function for extracting daily PRISM data from provisional files
def getDailyPrismProvis(year, month, metdata, data_path, bil_name, coords_file):
    # Read in site coordinates, get date range and create a DataFrame
    #to fill
    pnts = pd.read_csv(coords_file)
    drange =  pd.date_range('{0}-1-1'.format(year),
            '{0}-12-31'.format(year), freq='D')
    drange = drange[drange.month==month]
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
        bil_file = (r'/vsizip/' + data_path + bil_name + '/' +
        r'PRISM_{0}_provisional_4kmD2_{1}{2}{3}_bil.bil'.format(metdata, 
            *ymd_tuple))
        bil_ds = gdalFile(bil_file)
        for j in range(len(pnts.index)):
            pt_val = bil_ds.extract_coord_val(pnts.lat[j], pnts.lon[j])
            df.iloc[i, j] = pt_val
    return df


# Function for extracting monthly PRISM data
# Note that this uses 1981-2015 files and will need to be altered if different
# files are used
def getMonthlyPrism( metdata, data_path, coords_file ):
    # Read in site coordinates, get date range and create a DataFrame
    #to fill
    pnts = pd.read_csv(coords_file)
    drange =  pd.date_range('1-1-1981', '9-30-2015', freq='M')
    df = pd.DataFrame(index=drange, columns=pnts.sitecode)
    for i in range(len(drange)):
        # Create a tuple to fill the file dates in,
        # pad month & day with zeros
        ym_tuple = (str(drange.year[i]),
                str(drange.month[i]).zfill(2))
        # Create GDAL vsizip file path (read directly from zip archive)
        # See https://trac.osgeo.org/gdal/wiki/UserDocs/ReadInZip
        # If for some reason this vsizip interface won't work, extract the
        # archive and remove the 'vsizip' part of pathname
        if metdata=='ppt':
            bil_file = (r'/vsizip/' + data_path +
                    r'PRISM_{0}_stable_4kmM3_198101_201509_bil.zip/'.format(
                        metdata) +
                    r'PRISM_{0}_stable_4kmM3_{1}{2}_bil.bil'.format(
                        metdata, *ym_tuple))
        elif metdata=='tmean':
            bil_file = (r'/vsizip/' + data_path +
                    r'PRISM_{0}_stable_4kmM2_198101_201509_bil.zip/'.format(
                        metdata) +
                    r'PRISM_{0}_stable_4kmM2_{1}{2}_bil.bil'.format(
                        metdata, *ym_tuple))

        bil_ds = gdalFile(bil_file)
        for j in range(len(pnts.index)):
            pt_val = bil_ds.extract_coord_val(pnts.lat[j], pnts.lon[j])
            df.iloc[i, j] = pt_val
    return df

# Function for extracting monthly PRISM data from provisional files
def getMonthlyPrismProvis(year, metdata, data_path, bil_name, coords_file):
    # Read in site coordinates, get date range and create a DataFrame
    #to fill
    pnts = pd.read_csv(coords_file)
    drange =  pd.date_range('{0}-10-01'.format(year),
            '{0}-12-31'.format(year), freq='M')
    #drange = drange[drange.month==month]
    df = pd.DataFrame(index=drange, columns=pnts.sitecode)
    for i in range(len(drange)):
        # Create a tuple to fill the file dates in,
        # pad month & day with zeros
        ym_tuple = (str(drange.year[i]),
                str(drange.month[i]).zfill(2))
        # Create GDAL vsizip file path (read directly from zip archive)
        # See https://trac.osgeo.org/gdal/wiki/UserDocs/ReadInZip
        # If for some reason this vsizip interface won't work, extract the
        # archive and remove the 'vsizip' part of pathname
        bil_file = (r'/vsizip/' + data_path + bil_name + '/' +
        r'PRISM_{0}_provisional_4kmM2_{1}{2}_bil.bil'.format(metdata, 
            *ym_tuple))
        #pdb.set_trace()
        bil_ds = gdalFile(bil_file)
        for j in range(len(pnts.index)):
            pt_val = bil_ds.extract_coord_val(pnts.lat[j], pnts.lon[j])
            df.iloc[i, j] = pt_val
    return df



# Function for extracting 30 year normal PRISM precip data
def get_800m_30ynormal(prism_path, lat, lon, datamonth='annual',
        datatype='ppt'):
    """
    Extract 800m 30 year normal PRISM data for a list of lat/lon coordinates.

    Args:
        prism_path (str) : path to directory containing prism zip files
        lat (float array): latitude coordinates 
        lon (float array : longitude coordinates
        datamonth (str)  : Monthly normal to return (eg. 'Jan') or 'annual'
        datatype (str)   : abbreviation for the PRISM data type
    Returns:
        A pandas dataframe with lat/lon coordinates and PRISM normal data

    Note that this is looking for 800m "all" normal files, which are
    downloaded here:
    
    http://www.prism.oregonstate.edu/normals/
    
    by using the "Download All Normals Data(.bil)" option. The downloaded files
    contain both the annual and monthly normal datasets in separate .bil files.
    To get monthly files provide the proper 3 letter month code to datamonth.
    """
    if datamonth is not 'annual':
        datamonth = str(strptime('Feb', '%b').tm_mon)

    df = pd.DataFrame(index=lat.index,
            columns=['lat','lon', datatype + '_PRISMn'])

    bil_file = r'/vsizip/' + os.path.join(prism_path,
            r'PRISM_{0}_30yr_normal_800mM2_all_bil.zip'.format(datatype), 
            r'PRISM_{0}_30yr_normal_800mM2_{1}_bil.bil'.format(datatype,
                datamonth))
    bil_obj = gf.gdalFile(bil_file)
    for j in range(len(lat)):
        df.loc[j, ['lat', 'lon']] = [lat[j], lon[j]]
        df.loc[j, datatype + '_PRISMn'] = bil_obj.extract_coord_val(
                lat[j], lon[j])
    return df


import shapefile as shp
import csv

def csv2shp( infile, outfile,

    Â'''
    https://gis.stackexchange.com/questions/35593/using-the-python-shape-library-pyshp-how-to-convert-csv-file-to-shp
    '''

#Set up blank lists for data
x,y,id_no,date,target=[],[],[],[],[]

#read data from csv file and store in lists
with open('input.csv', 'rb') as csvfile:
    r = csv.reader(csvfile, delimiter=';')
    for i,row in enumerate(r):
        if i > 0: #skip header
            x.append(float(row[3]))
            y.append(float(row[4]))
            id_no.append(row[0])
            date.append(''.join(row[1].split('-')))#formats the date correctly
            target.append(row[2])

#Set up shapefile writer and create empty fields
w = shp.Writer(shp.POINT)
w.autoBalance = 1 #ensures gemoetry and attributes match
w.field('X','F',10,8)
w.field('Y','F',10,8)
w.field('Date','D')
w.field('Target','C',50)
w.field('ID','N')

#loop through the data and write the shapefile
for j,k in enumerate(x):
    w.point(k,y[j]) #write the geometry
    w.record(k,y[j],date[j], target[j], id_no[j]) #write the attributes

#Save shapefile
w.save(outfile)

from libGMAOMERRA2Data import *

import os
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

import h5py


# Select data


#DATAFIELD_NAME =  'TO3'
#DATAFIELD_UNIT = DATAFIELD_NAME+' (Ozone:Db) '

DATAFIELD_NAME =  'TQV'
DATAFIELD_UNIT = DATAFIELD_NAME+' (PWV:mm) '


os.environ["HDFEOS_ZOO_DIR"] = "/Users/dagoret/DATA/MERRA2/inst1_2d_asm_Nx_M2I1NXASM/2025"
#os.environ["HDFEOS_ZOO_DIR"] =  "/sps/lsst/groups/auxtel/MERRA2/data/inst1_2d_asm_Nx_M2I1NXASM/2025"

# If a certain environment variable is set, look there for the input
# file, otherwise look in the current directory.

hdffile = 'MERRA2_400.inst1_2d_asm_Nx.20250101.nc4'


FILE_NAME= hdffile

base_filename=os.path.basename(FILE_NAME).split('.hdf')[0]
p = re.compile('[.]')
root_filename=p.sub('_',base_filename)
rootimg_dir=os.path.join('test_images',root_filename)

try:
    FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
except KeyError:
    pass


# Access to data
#
(data3D,unit,longname) = GetGeoRefData(FILE_NAME,DATAFIELD_NAME)

data= data3D[0,:,:] ## the first index is the time, for example take the first one

(lat,un_lat,nm_lat) = Get1DData(FILE_NAME,'lat')
latitude = lat[:]
(lon,un_lon,nm_lon) = Get1DData(FILE_NAME,'lon')
longitude = lon[:]
(tim,un_tim,nm_tim) = Get1DData(FILE_NAME,'time')
thetime=tim[:]# Open the file


# Plots maps

X,Y=np.meshgrid(longitude,latitude)
longitude=X
latitude=Y


PlotData(longitude,latitude,data,12,6,title=base_filename,labelz=longname,
        longs=all_longs,
        lats=all_lats,
        tags=all_tags)



PlotGeoData(longitude,latitude,data,12,8,title=base_filename,labelz=longname,
            longs=all_longs,
            lats=all_lats,
            tags=all_tags)


# Select area
LongMin=-90
LongMax=-30
LatMin=-55
LatMax=15

PlotGeoData2(longitude,latitude,data,LatMin,LatMax,LongMin,LongMax,10,14,title=base_filename,labelz=longname,
            longs=all_longs,
            lats=all_lats,
            tags=all_tags)



# Select area
LongMin=-80
LongMax=-60
LatMin=-40
LatMax=-20

PlotGeoData2(longitude,latitude,data,LatMin,LatMax,LongMin,LongMax,14,8,title=base_filename,labelz=longname,
            longs=all_longs,
            lats=all_lats,
            tags=all_tags)




  # Select area
LongMin=-75
LongMax=-65
LatMin=-35
LatMax=-25

PlotGeoData2(longitude,latitude,data,LatMin,LatMax,LongMin,LongMax,14,8,title=base_filename,labelz=longname,
            longs=all_longs,
            lats=all_lats,
            tags=all_tags)


# Select area
LongMin=-73
LongMax=-67
LatMin=-33
LatMax=-27

PlotGeoData2(longitude,latitude,data,LatMin,LatMax,LongMin,LongMax,14,8,title=base_filename,labelz=longname,
            longs=all_longs,
            lats=all_lats,
            tags=all_tags)


# Select area
LongMin=-130
LongMax=-100
LatMin=20
LatMax=50

PlotGeoData2(longitude,latitude,data,LatMin,LatMax,LongMin,LongMax,14,8,title=base_filename,labelz=longname,
            longs=all_longs,
            lats=all_lats,
            tags=all_tags)


# Select area
LongMin=-120
LongMax=-110
LatMin=30
LatMax=40

PlotGeoData2(longitude,latitude,data,LatMin,LatMax,LongMin,LongMax,14,8,title=base_filename,labelz=longname,
            longs=all_longs,
            lats=all_lats,
            tags=all_tags)


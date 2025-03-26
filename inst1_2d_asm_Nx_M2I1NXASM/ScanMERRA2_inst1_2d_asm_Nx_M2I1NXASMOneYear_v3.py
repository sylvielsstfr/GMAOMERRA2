#!/usr/bin/env python
# coding: utf-8

# # Scan MERRA-2 atmospheric properties during one Year
# ----------------------------------------------------------------------------------
#
# - author: Sylvie Dagoret-Campagne
# - creation January 12 2017
# - last update February 17 2023 on imac
# - last update January 2024/01/28 at CC with kernel ``conda_jax0325_py310``
# - last update January 2024/10/15 at CC with kernel ``conda_jax0325_py310``
# - at CCIN2P3 anaconda3_py39_auxtel kernel (kernel not working anymore)
# - loal Python 3 kernel
#
# Link:
#
# http://disc.sci.gsfc.nasa.gov/datareleases/merra_2_data_release
#
# ### purpose:
#
# Scan One month of MERRA-2 predictions of the dataset inst1_2d_asm_Nx_M2I1NXASM.
# Extract the relevant atmospheric variables.
# Build the correcponding time series and dataset in pandas.
# Plot the variables. Save the pandas dataset into a file.
# Convert the pandas dataset into an astropy fits table and save into a fits file as well.

# ## 1) python libraries
# ---------------------------

# Set up matplotlib and use a nicer set of plot parameters
#get_ipython().run_line_magic('config', 'InlineBackend.rc = {}')
import matplotlib
import matplotlib as mpl
matplotlib.rc_file("templates/matplotlibrc")
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')


import datetime

from matplotlib.dates import MonthLocator, WeekdayLocator,DateFormatter
from matplotlib.dates import MONDAY


mondays = WeekdayLocator(MONDAY)
months = MonthLocator(range(1, 13), bymonthday=1, interval=1)
monthsFmt = DateFormatter("%b '%y")



import os

def default_proj_lib():
    proj_lib = os.getenv('PROJ_LIB')
    if proj_lib not in (None, 'PROJ_LIB'):
        return proj_lib
    try:
        import conda
    except ImportError:
        conda = None
    if conda is not None or os.getenv('CONDA_PREFIX') is None:
        conda_file_dir = conda.__file__
        conda_dir = conda_file_dir.split('lib')[0]
        proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
        if os.path.exists(proj_lib):
            return proj_lib
        return None
    return None




theprojlib=default_proj_lib()
print("The PROJ_LIB={}".format(theprojlib))



import os
import re
import numpy as np
#from mpl_toolkits.basemap import Basemap
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord

from astropy.table import Table

import h5py
import libGMAOMERRA2Data as merra2  # My own library


############################################################################
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(f):
        os.makedirs(f)
#########################################################################


# ## 2)  Configuration
# -------------------------



#YEARNUM="2023"
YEARNUM="2024"


# In[15]:


# SELECT OBSERVATORY
OBS_NAME='lsst'


# where are the HDF files
#HDFEOS_ZOO_DIR="/Volumes/LaCie2/DATA/MERRA-2/inst1_2d_asm_Nx_M2I1NXASM"
# For Simulation of CTIO atmosphere in May-Jun 2017
#HDFEOS_ZOO_DIR="/sps/lsst/data/AtmosphericCalibration/MERRA-2/May-Jun-2017/subset_M2I1NXASM_V5.12.4_20180424_201411"
# Data 2018
#HDFEOS_ZOO_DIR="/sps/lsst/data/MERRA2/M2I1NXASM.5.12.4/"+YEARNUM
#HDFEOS_ZOO_DIR="/Users/dagoret/DATA/MERRA-2/M2I1NXASM_V5.12.4/"+YEARNUM
#HDFEOS_ZOO_DIR=os.path.join("/Users/dagoret/DATA/MERRA-2/M2I1NXASM_V5.12.4",YEARNUM)
HDFEOS_ZOO_DIR=os.path.join("/sps/lsst/groups/auxtel/MERRA2/data/inst1_2d_asm_Nx_M2I1NXASM",YEARNUM)
#HDFEOS_ZOO_DIR=os.path.join("/Users/dagoret/DATA/merra-2import/M2I1NXASM.5.12.4",YEARNUM)
#HDFEOS_ZOO_DIR=os.path.join("/Users/sylvie/DATA/MERRA2/data/inst1_2d_asm_Nx_M2I1NXASM",YEARNUM)


# In[17]:


path=HDFEOS_ZOO_DIR


# ### Here I describe the content of the input files



DATA_TAG=['PS','T10M','T2M','TO3','TOX','TQI','TQL','TQV','TS','U10M','U2M','U50M','V10M','V2M','V50M' ]


# In[19]:


DATA_TITLE=['Surface Pressure',
            '10 meter air temperature',
            '2 meter air temperature',
            'total column ozone',
            'total column odd oxygen',
            'total precipitable ice water',
            'total precipitable liquid water',
            'total precipitable water vapor',
            'Surface Temperature skin',
            '10 meter eastward wind',
            '2 meter eastward wind',
            '50 meter eastward wind',
            '10 meter northward wind',
            '2 meter northward wind',
            '50 meter northward wind'
           ]

NB_DATAFIELDS=len(DATA_TAG)


# ### List of output files

# The selected data field
DATA_NAME =  'inst1_2d_asm_Nx_M2I1NXASM'   #


pandas_filename='MERRA2_'+YEARNUM+'_'+DATA_NAME+'_'+OBS_NAME+'_'+'AllYear'+'.csv'
fits_filename='MERRA2_'+YEARNUM+'_'+DATA_NAME+'_'+OBS_NAME+'_'+'AllYear' +'.fits'
hdf5_filename='MERRA2_'+YEARNUM+'_'+DATA_NAME+'_'+OBS_NAME+'_'+'AllYear'+'.h5'


# ### Select where in the world


# Select observatory
loc=merra2.observatory_location(OBS_NAME)


# ### 2.2) Getting the list of the files
# ------------------------------


nc4_files = [f for f in os.listdir(path) if f.endswith('.nc4')]

print(nc4_files[:5])



keysel_filename='^MERRA2_400.inst1_2d_asm_Nx.'+YEARNUM+'.*'


print('Selection key' ,keysel_filename)

nc4_files2 = []
for file in nc4_files:
    if re.findall(keysel_filename,file):
        nc4_files2.append(file)

nc4_files2=np.array(nc4_files2)



# ### 2.4) Sort files by increasing time

nc4_files=np.sort(nc4_files2)


# ### 2.5) Build the full filename before reading


NBFILES=len(nc4_files)
full_nc4files=[]

for file in nc4_files:
    fname = os.path.join(path, file)
    full_nc4files.append(fname)


# ## 3)  Extract data and write them into pandas dataset and time series
# --------------------------------------------------------------------------------------

# In[ ]:


ts0=[]  # intermediate data series
ts1=[]
ts2=[]
ts3=[]  # intermediate data series
ts4=[]
ts5=[]
ts6=[]  # intermediate data series
ts7=[]
ts8=[]
ts9=[]  # intermediate data series
ts10=[]
ts11=[]
ts12=[]
ts13=[]
ts14=[]

df_inst1_2d_asm_Nx=[] # final pandas dataset for all atmospheric quantities

for file in full_nc4files: # loop on data file of each day of the month

    print(file)

    #Retrieve 1D parameters longitude, latitude, time
    (m_lat,m_un_lat,m_nm_lat) = merra2.Get1DData(file,'lat') # latitude (array, unit, name)
    m_latitude = m_lat[:]
    (m_lon,m_un_lon,m_nm_lon) = merra2.Get1DData(file,'lon') # longitude(array, unit, name)
    m_longitude = m_lon[:]
    (m_tim,m_un_tim,m_nm_tim)= merra2.Get1DData(file,'time') # time (array, unit, name)
    m_time=m_tim[:]

    # with python3 obliged to transform byte string into a string
    m_un_tim2=m_un_tim.decode("utf-8")

    NbDataPerFile=m_time.shape[0] # number of data sample per file
    #start_time = re.findall("^minutes since[ ]([0-9.].+[0-9.].+[0-9.].+)[ ]00:00:00$",m_un_tim) # extract start time
    start_time = re.findall("^minutes since[ ]([0-9.].+[0-9.].+[0-9.].+)",m_un_tim2) # extract start time

    #print 'start_time = ', start_time
    time_rng = pd.date_range(start_time[0], periods=NbDataPerFile, freq='H') # one data per hour



    m_X,m_Y=np.meshgrid(m_longitude,m_latitude) # build meash-grid in longitude and latitude
    (sel_long, sel_lat)=merra2.GetBinIndex(m_X,m_Y,loc[0],loc[1]) # get bin in longitude and latitude for the site


    # loop
    for index in range(NB_DATAFIELDS):
        (m_data,m_unit,m_longname)=merra2.GetGeoRefData(file,DATA_TAG[index]) # 3D array : time x longitude x latitude
        dt=m_data[:,sel_lat,sel_long]
        if index==0:
            ts0 = pd.Series(dt, index=time_rng)
        elif index==1:
            ts1 = pd.Series(dt, index=time_rng)
        elif index==2:
            ts2 = pd.Series(dt, index=time_rng)
        elif index==3:
            ts3 = pd.Series(dt, index=time_rng)

        elif index==4:
            ts4 = pd.Series(dt, index=time_rng)
        elif index==5:
            ts5 = pd.Series(dt, index=time_rng)
        elif index==6:
            ts6 = pd.Series(dt, index=time_rng)

        elif index==7:
            ts7 = pd.Series(dt, index=time_rng)
        elif index==8:
            ts8 = pd.Series(dt, index=time_rng)
        elif index==9:
            ts9 = pd.Series(dt, index=time_rng)

        elif index==10:
            ts10 = pd.Series(dt, index=time_rng)
        elif index==11:
            ts11 = pd.Series(dt, index=time_rng)
        elif index==12:
            ts12 = pd.Series(dt, index=time_rng)

        elif index==13:
            ts13 = pd.Series(dt, index=time_rng)
        elif index==14:
            ts14 = pd.Series(dt, index=time_rng)


        #clf_timeseries.append(ts)
        # Create the dataframe
    df = pd.DataFrame({DATA_TAG[0]: ts0,
                       DATA_TAG[1]: ts1,
                       DATA_TAG[2]: ts2,
                       DATA_TAG[3]: ts3,

                       DATA_TAG[4]: ts4,
                       DATA_TAG[5]: ts5,
                       DATA_TAG[6]: ts6,

                       DATA_TAG[7]: ts7,
                       DATA_TAG[8]: ts8,
                       DATA_TAG[9]: ts9,

                       DATA_TAG[10]: ts10,
                       DATA_TAG[11]: ts11,
                       DATA_TAG[12]: ts12,
                       DATA_TAG[13]: ts13,
                       DATA_TAG[14]: ts14 }, index=time_rng)
    df_inst1_2d_asm_Nx.append(df)



# ### Concatenation

df_inst1_2d_asm_Nx=pd.concat(df_inst1_2d_asm_Nx)


print(df_inst1_2d_asm_Nx.info())


# ## 4) Plot the time dependence of time series and dataset

# ####  Fill information on years and months from the filenames


# ## 5) Output


df_inst1_2d_asm_Nx.index.name='time'
print(df_inst1_2d_asm_Nx.describe())

# ## 5)  Save dataset  in file pandas (csv)
# ----------------------------------------


dataset=df_inst1_2d_asm_Nx
dataset.index.name='time'

print(dataset.describe())
print(dataset.head())


# ### write pandas

dataset.to_csv(pandas_filename)




# ### Check


saved_dataset=pd.read_csv(pandas_filename)



# ## 6) Convert dataset into a table and then save in a fits file
# --------------------------------------------------------------------------

table = Table.from_pandas(saved_dataset)
table.write(fits_filename,format='fits',overwrite=True)


#





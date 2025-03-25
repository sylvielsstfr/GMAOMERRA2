#!/usr/bin/env python
# coding: utf-8

# # Scan MERRA-2 atmospheric properties during one year
# ----------------------------------------------------------------------------------
# 
# - author: Sylvie Dagoret-Campagne
# - creation January 12 2017
# - update January 12 2016
# - update April 25th 2018
# - update June 5th 2019
# - update December 2022
# - update 2024/10/15
# - last update January 2024 28 at CC with kernel ``conda_jax0325_py310``
# - update 2025/03/25 with anaconda3-py311 (2025) at CC
# 
# 
# Link:
# 
# http://disc.sci.gsfc.nasa.gov/datareleases/merra_2_data_release
# 
# ### purpose:
# 
# Scan One year of MERRA-2 predictions of the dataset tavg1_2d_aer_Nx_M2T1NXAER. 
# Extract the relevant atmospheric variables.
# Build the correcponding time series and dataset in pandas.
# Plot the variables. Save the pandas dataset into a file.
# Convert the pandas dataset into an astropy fits table and save into a fits file as well.
# 
# 

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

# ### SELECT YEAR

YEARNUM="2024"


# SELECT OBSERVATORY
OBS_NAME='lsst'

#where are the HDF files
#HDFEOS_ZOO_DIR="/Volumes/DAGORETBACK/MERRA-2/inst1_2d_asm_Nx_M2I1NXASM/2016"
#HDFEOS_ZOO_DIR="/Volumes/LaCie2/DATA/MERRA-2/tavg1_2d_aer_Nx_M2T1NXAER"
#HDFEOS_ZOO_DIR="/sps/lsst/data/AtmosphericCalibration/MERRA-2/May-Jun-2017/subset_M2T1NXAER_V5.12.4_20180424_195815"
#HDFEOS_ZOO_DIR="/Users/dagoret/DATA/MERRA-2/M2T1NXAER.5.12.4"
#HDFEOS_ZOO_DIR=os.path.join("/Users/dagoret/DATA/MERRA-2/M2T1NXAER.5.12.4",YEARNUM)
#HDFEOS_ZOO_DIR=os.path.join("/Users/dagoret/DATA/merra-2import/M2T1NXAER.5.12.4",YEARNUM)
HDFEOS_ZOO_DIR=os.path.join("/sps/lsst/groups/auxtel/MERRA2/data/tavg1_2d_aer_Nx_M2T1NXAER",YEARNUM)

path=HDFEOS_ZOO_DIR


# ### Here I describe the content of the input files

DATA_TAG=['TOTANGSTR','TOTEXTTAU','TOTSCATAU']

DATA_TITLE=['Total Aerosol Angstrom parameter 470-870 nm',
            'Total Aerosol Extinction AOT 550 nm',
            'Total Aerosol Scattering AOT 550 nm'
           ]

NB_DATAFIELDS=len(DATA_TAG)

# The selected data field
DATA_NAME =  'tavg1_2d_aer_Nx_M2T1NXAER'   # 

pandas_filename='MERRA2_'+YEARNUM+'_'+DATA_NAME+'_'+OBS_NAME+'_'+'AllYear'+'.csv'
fits_filename='MERRA2_'+YEARNUM+'_'+DATA_NAME+'_'+OBS_NAME+'_'+'AllYear' +'.fits'



# ### Select where in the world

# Select observatory
loc=merra2.observatory_location(OBS_NAME)


# ### 2.2) Getting the list of the files
# ------------------------------

nc4_files = [f for f in os.listdir(path) if f.endswith('.nc4')]  

print(nc4_files[:5])


# ### 2.3) Select files of a given month

keysel_filename='^MERRA2_400.tavg1_2d_aer_Nx.'+YEARNUM+'.*'

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

ts0=[]  # intermediate data series
ts1=[]
ts2=[]


df_tavg1_2d_aer_Nx=[] # final pandas dataset for all atmospheric quantities

for file in full_nc4files: # loop on data file of each day of the month
    
    #Retrieve 1D parameters longitude, latitude, time
    try:
        (m_lat,m_un_lat,m_nm_lat) = merra2.Get1DData(file,'lat') # latitude (array, unit, name)
        m_latitude = m_lat[:]
        (m_lon,m_un_lon,m_nm_lon) = merra2.Get1DData(file,'lon') # longitude(array, unit, name)
        m_longitude = m_lon[:]
        (m_tim,m_un_tim,m_nm_tim)= merra2.Get1DData(file,'time') # time (array, unit, name)
        m_time=m_tim[:]
    
    except Exception as inst:
        print(type(inst))    # the exception type
        print(inst.args)     # arguments stored in .args
        print(inst)  
        print("SKIP")
        continue
    
    # with python3 obliged to transform byte string into a string    
    m_un_tim2=m_un_tim.decode("utf-8")  
    
    NbDataPerFile=m_time.shape[0] # number of data sample per file
    #start_time = re.findall("^minutes since[ ]([0-9.].+[0-9.].+[0-9.].+)[ ]00:00:00$",m_un_tim) # extract start time
    start_time = re.findall("^minutes since[ ]([0-9.].+[0-9.].+[0-9.].+)",m_un_tim2) # extract start time
    
    #print 'start_time = ', start_time
    time_rng = pd.date_range(start_time[0], periods=NbDataPerFile, freq='H') # one data per hour
    
    #print('---------------------------------------------')
    #print('NbDataPerFile = ', NbDataPerFile)
    #print('start_time = ', start_time)
    #print('time_rng   = ', time_rng[:5])
    
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

            
        #clf_timeseries.append(ts)
        # Create the dataframe
    df = pd.DataFrame({DATA_TAG[0]: ts0, 
                       DATA_TAG[1]: ts1,
                       DATA_TAG[2]: ts2 }, index=time_rng)
    df_tavg1_2d_aer_Nx.append(df)  
    
# ### Concatenation

df_tavg1_2d_aer_Nx=pd.concat(df_tavg1_2d_aer_Nx)


print(df_tavg1_2d_aer_Nx.info())



# ## 5) Output

df_tavg1_2d_aer_Nx.index.name='time'
print(df_tavg1_2d_aer_Nx.describe())

# ## 5)  Save dataset  in file pandas (csv)
# ----------------------------------------


dataset=df_tavg1_2d_aer_Nx

dataset.index.name='time'


print(dataset.describe())

print(dataset.head())


dataset.to_csv(pandas_filename)
saved_dataset=pd.read_csv(pandas_filename)


print(saved_dataset.head())


# ## 6) Convert dataset into a table and then save in a fits file
# --------------------------------------------------------------------------

table = Table.from_pandas(saved_dataset)

table.write(fits_filename,format='fits',overwrite=True)


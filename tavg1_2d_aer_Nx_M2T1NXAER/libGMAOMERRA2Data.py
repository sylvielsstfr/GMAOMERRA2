"""
------------------------------------------------------------------------------
libGMAOLERRA2Data.py
=================

Author : Sylvie Dagoret-Campagne
Date   : November 25 2016



-----------------------------------------------------------------------------
"""

import os
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

import h5py

 #LSST site
Longitude_lsst = -70.7366833333333 # deg
Latitude_lsst = -30.240741666666672 #deg
Altitude_lsst = 2749.999999999238 #m
 #CTIO Site
Longitude_ctio = -70.815 # deg
Latitude_ctio = -30.165277777777778 #deg
Altitude_ctio = 2214.9999999993697 #m
# Cerro Paranal
 
Longitude_paranal = -70.40300000000002 #deg
Latitude_paranal  = -24.625199999999996 #deg
Altitude_paranal = 2635.0000000009704 #m
# Observatoire de Haute Provence
Longitude_ohp=5.71222222222
Latitude_ohp=43.9316666667
Altitude_ohp=650.    
    

#--------------------------------------------------------------------------
def ensure_dir(f):
    '''
    ensure_dir : check if the directory f exist. If not, it is created
    '''
    d = os.path.dirname(f)
    if not os.path.exists(f):
        os.makedirs(f)
#-----------------------------------------------------------------------------    

def loc_ctio():
    return(Longitude_ctio,Latitude_ctio,Altitude_ctio)
    
def loc_lsst():
    return(Longitude_lsst,Latitude_lsst,Altitude_lsst)
    
def loc_ohp():
    return(Longitude_ohp,Latitude_ohp,Altitude_ohp)
    
def loc_none():
    return(0,0,0)
#--------------------------------------------------------------------    
    
#--------------------------------------------------------------------    
def observatory_location(obs):
    if obs== 'ctio':
        loc=loc_ctio()
    elif obs=='lsst':
        loc=loc_lsst()
    elif obs=='ohp':
        loc=loc_ohp()
    else:
        loc=loc_none()
    return loc
#--------------------------------------------------------------------
                                                                                                                 

#---------------------------------------------------------------------------------
def GetGeoRefData(file,datafield):
    '''
    GetGeoRefData(file,datafield) : read the data labeled datafield in file
    =================================================================
    
    Retrieve data from a HDF file    
    
    input : 
    ------
        file : name of input file
        datafield : label of required data field
    output:
    ------
        data3D : output array
    
    '''
    
    nc4f= h5py.File(file, mode='r') 
    data = nc4f[datafield][:,:,:]
    units = nc4f[datafield].attrs['units']
    long_name = nc4f[datafield].attrs['long_name']
    _FillValue = nc4f[datafield].attrs['_FillValue']
    data[data == _FillValue] = np.nan
    data = np.ma.masked_where(np.isnan(data), data)
    #nc4f.close()
    
    return data,units,long_name
#-----------------------------------------------------------------------------    


#---------------------------------------------------------------------------------
def Get1DData(file,datafield):
    '''
    Get1DData(file,datafield) : read the data labeled datafield in file
    =================================================================
    
    Retrieve data from a HDF file    
    
    input : 
    ------
        file : name of input file
        datafield : label of required data field
    output:
    ------
        data3D : output array
    
    '''
    
    nc4f= h5py.File(file, mode='r') 
    data = nc4f[datafield]
    units = nc4f[datafield].attrs['units']
    long_name = nc4f[datafield].attrs['long_name']
    #nc4f.close()
    
    
    return data,units,long_name
#-----------------------------------------------------------------------------   


#--------------------------------------------------------------------------------
def AreaSelect(X,Y,data,LongMin,LongMax,LatMin,LatMax):
    '''
    AreaSelect(X,Y,data,LongMin,LongMax,LatMin,LatMax)
    ==================================================
    
    Select an area    
    
    input:
    ------
        X,Y  : Longitude and lattitude 2D array 
        data : Data array
        LongMin,LongMax,LatMin,LatMax : Longitude and Latitude boundaries
    
    output:
    -------
        Xsel,Ysel : Longitude and lattitude 2D array of selected zone
        extracted_data : data array selected
    
    '''
    flags_long=np.logical_and(X>=LongMin, X<=LongMax)   # flags in X where are the selected longitudes
    flags_lat=np.logical_and(Y>=LatMin, Y<=LatMax)      # flags in Y where are the selected longitudes
    flags_longlat=np.logical_and(flags_long,flags_lat)  # flags where the region is selected in the long-lat matrix

    (selected_lat_indexes,selected_long_indexes)=np.where(flags_longlat==True) # list of indexes
    selected_long=longitude[:,selected_long_indexes] # all selected longitudes
    selected_lat=latitude[selected_lat_indexes,:]    # all selected latitudes

    min_long_index=np.min(selected_long_indexes)
    max_long_index=np.max(selected_long_indexes)

    min_lat_index=np.min(selected_lat_indexes)
    max_lat_index=np.max(selected_lat_indexes)

    # output    
    extracted_data=data[min_lat_index:max_lat_index,min_long_index:max_long_index] # extract the data
    Xsel=X[min_lat_index:max_lat_index,min_long_index:max_long_index] # extract the Long
    Ysel=Y[min_lat_index:max_lat_index,min_long_index:max_long_index] # extract the lat
    
    return Xsel,Ysel,extracted_data
#---------------------------------------------------------------------------------    
    
#---------------------------------------------------------------------------------    
def SelectBin(X,Y,data,Long0,Lat0,DLong=0.625,DLat=0.498614958449):
    '''
    SelectBin(X,Y,data,Long0,Lat0,DLong=1.0,DLat=1.0)
    =================================================
    
    Select one bin    
    
    input:
    -----
        X,Y        : Longitude and Latitude 2D array
        data       : 2D array of data 
        Long0,Lat0 : The Longitude and Latitude of the bin
        DLong,DLat : The angular bin width
    
    output:
    ------
        sel_min_long_index,sel_min_lat_index : index of selected bin 
        extracted_data                       : data of the selected bin
    
    '''
    sel_flags_long=np.logical_and(X>=Long0-float(DLong)/2., X<=Long0+float(DLong)/2.)   # flags in X where are the selected longitudes
    sel_flags_lat=np.logical_and(Y>=Lat0-float(DLat)/2., Y<=Lat0+float(DLat)/2.)      # flags in Y where are the selected longitudes
    sel_flags_longlat=np.logical_and(sel_flags_long,sel_flags_lat)  # flags where the region is selected in the long-lat matrix

    (selected_lat_indexes,selected_long_indexes)=np.where(sel_flags_longlat==True) # list of indexes


    selected_X=X[:,selected_long_indexes] # all selected longitudes
    selected_Y=Y[selected_lat_indexes,:] 

    sel_min_long_index=np.min(selected_long_indexes)
    sel_max_long_index=np.max(selected_long_indexes)

    sel_min_lat_index=np.min(selected_lat_indexes)
    sel_max_lat_index=np.max(selected_lat_indexes)

    extracted_data=data[sel_min_lat_index:sel_max_lat_index+1,sel_min_long_index:sel_max_long_index+1] # extract the data
    
    return sel_min_long_index,sel_min_lat_index,extracted_data[0][0]    
#---------------------------------------------------------------------------------    
    
  #---------------------------------------------------------------------------------    
def GetBinIndex(X,Y,Long0,Lat0,DLong=0.625,DLat=0.498614958449):
    '''
    GetBinIndex(X,Y,Long0,Lat0,DLong=1.0,DLat=1.0)
    =================================================
    
    Select one bin    
    
    input:
    -----
        X,Y        : Longitude and Latitude 2D array
        Long0,Lat0 : The Longitude and Latitude of the bin
        DLong,DLat : The angular bin width
    
    output:
    ------
        sel_min_long_index,sel_min_lat_index : index of selected bin 
        extracted_data                       : data of the selected bin
    
    '''
    sel_flags_long=np.logical_and(X>=Long0-float(DLong)/2., X<=Long0+float(DLong)/2.)   # flags in X where are the selected longitudes
    sel_flags_lat=np.logical_and(Y>=Lat0-float(DLat)/2., Y<=Lat0+float(DLat)/2.)      # flags in Y where are the selected longitudes
    sel_flags_longlat=np.logical_and(sel_flags_long,sel_flags_lat)  # flags where the region is selected in the long-lat matrix

    (selected_lat_indexes,selected_long_indexes)=np.where(sel_flags_longlat==True) # list of indexes


    selected_X=X[:,selected_long_indexes] # all selected longitudes
    selected_Y=Y[selected_lat_indexes,:] 

    sel_min_long_index=np.min(selected_long_indexes)
    sel_max_long_index=np.max(selected_long_indexes)

    sel_min_lat_index=np.min(selected_lat_indexes)
    sel_max_lat_index=np.max(selected_lat_indexes)

    #extracted_data=data[sel_min_lat_index:sel_max_lat_index+1,sel_min_long_index:sel_max_long_index+1] # extract the data
    
    return sel_min_long_index,sel_min_lat_index   
#---------------------------------------------------------------------------------   
    
#---------------------------------------------------------------------------------
def PlotData(X,Y,data,sizex=8,sizey=8,labelx='longitude',labely='latitude',labelz='Unit',title=''):
    '''
    PlotData(X,Y,data,sizex=8,sizey=8,labelx='longitude',labely='latitude',labelz='Unit',title='')
    ==============================================================================================

    Plot in matplotlib the 2D array of data

    input:
    ------
        X,Y   : 2D array of lontitude and latitude
        data  : Data array
        sizex,sizey   :  size of figure
        labelx,labely,labelz,title : labels of axis and title
        
    output:
        figure in matplotlib
    
    '''


    fig=plt.figure(figsize=(sizex,sizey))
    im = plt.pcolormesh(X,Y,data)
    cbar=plt.colorbar(im, orientation='vertical')
    cbar.set_label(labelz)
    Xmin=X.min()
    Xmax=X.max()
    Ymin=Y.min()
    Ymax=Y.max()
    plt.axis([Xmin, Xmax,Ymin,Ymax])
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(title)
    #plt.tight_layout()
    plt.show()
#------------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
if __name__ == "__main__":
    
    DATAFIELD_NAME =  'TO3'
    DATAFIELD_UNIT = DATAFIELD_NAME+' (Ozone:Db) '
    
    #os.environ["HDFEOS_ZOO_DIR"] = "/Users/dagoret-campagnesylvie/MacOSX/LSST/MyWork/GitHub/GMAOMERRA2data/inst1_2d_asm_Nx_M2I1NXASM"
    
    os.environ["HDFEOS_ZOO_DIR"] = "/sps/lsst/data/AtmosphericCalibration/MERRA-2/May-Jun-2017/subset_M2I1NXASM_V5.12.4_20180424_201411"

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    #hdffile = 'MERRA2_400.inst1_2d_asm_Nx.20161001.nc4'
    hdffile = 'MERRA2_400.inst1_2d_asm_Nx.20170603.nc4'
    
   
    FILE_NAME= hdffile
    
    base_filename=os.path.basename(FILE_NAME).split('.hdf')[0]
    p = re.compile('[.]')
    root_filename=p.sub('_',base_filename)    
    rootimg_dir=os.path.join('test_images',root_filename)
    
    try:
        FILE_NAME = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    
    
    (data3D,unit,longname)=GetGeoRefData(FILE_NAME,DATAFIELD_NAME)
    data= data3D[0,:,:] ## Ozone has no additional dimensions
    (lat,un_lat,nm_lat) = Get1DData(FILE_NAME,'lat')
    latitude = lat[:]
    (lon,un_lon,nm_lon) = Get1DData(FILE_NAME,'lon')
    longitude = lon[:]
    (tim,un_tim,nm_tim)=Get1DData(FILE_NAME,'time')
    thetime=tim[:]
    
    
    X,Y=np.meshgrid(longitude,latitude)
    longitude=X
    latitude=Y


    # Plot world data
    #-------------------------
    PlotData(longitude,latitude,data,7,3.5,title=base_filename,labelz=longname)
    
    # Select area
    LongMin=-100
    LongMax=-30
    LatMin=-55
    LatMax=15  
    
    X=longitude
    Y=latitude
    
    (Xsel,Ysel,extracted_data) = AreaSelect(X,Y,data,LongMin,LongMax,LatMin,LatMax)
    # plot area
    #------------
    PlotData(Xsel,Ysel,extracted_data,6,6,title=base_filename,labelz=longname)
    
    

    
    #LSST site
    Longitude_lsst = -70.7366833333333 # deg
    Latitude_lsst = -30.240741666666672 #deg
    Altitude_lsst = 2749.999999999238 #m
    #CTIO Site
    Longitude_ctio = -70.815 # deg
    Latitude_ctio = -30.165277777777778 #deg
    Altitude_ctio = 2214.9999999993697 #m
    # Cerro Paranal
    Longitude_paranal = -70.40300000000002 #deg
    Latitude_paranal  = -24.625199999999996 #deg
    Altitude_paranal = 2635.0000000009704 #m
    # Observatoire de Haute Provence
    Longitude_ohp=5.71222222222
    Latitude_ohp=43.9316666667
    Altitude_ohp=650.    
    
    
    
    # Select one bin
    #-----------------
    (ctio_min_long_index, ctio_min_lat_index,extrdata)=SelectBin(X,Y,data,Longitude_ctio,Latitude_ctio)
        
    print('ctio_min_lat_index=',ctio_min_lat_index)
    print('ctio_min_long_index=',ctio_min_long_index)
    print('ctio_data = ',extrdata,DATAFIELD_UNIT)
    
    
    (lsst_min_long_index, lsst_min_lat_index,extrdata)=SelectBin(X,Y,data,Longitude_lsst,Latitude_lsst)
        
    print('lsst_min_lat_index=',lsst_min_lat_index)
    print('lsst_min_long_index=',lsst_min_long_index)
    print('lsst_data = ',extrdata,DATAFIELD_UNIT)
    
    
    (ohp_min_long_index, ohp_min_lat_index,extrdata)=SelectBin(X,Y,data,Longitude_ohp,Latitude_ohp)
    
    print('ohp_min_lat_index=',ohp_min_lat_index)
    print('ohp_min_long_index=',ohp_min_long_index)
    print('ohp_data = ',extrdata,DATAFIELD_UNIT)
    
    
    
#---------------------------------------------------------------------------------
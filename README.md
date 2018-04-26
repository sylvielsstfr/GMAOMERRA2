# GMAOMERRA2
Access to MERRA2 runs, to extract atmospheric quantities relevant for CTIO/LSST atmospheric attenuation

- author : Sylvie Dagoret-Campagne
- date  : April 26th 2018
- affiliation : LAL/IN2P3/CNRS


## inst1\_2d\_asm\_Nx\_M2I1NXASM :
Notebooks and scripts to decode MERRA-2 data for the dataset **inst1\_2d\_asm\_Nx\_M2I1NXASM**
This dataset contains:

- pressure
- pwv
- ozone

Subsets are extracted as *.csv* files.


## tavg1\_2d\_aer\_Nx\_M2T1NXAER :
Notebooks and scripts to decode MERRA-2 data for the dataset **tavg1\_2d\_aer\_Nx\_M2T1NXAER**
This dataset contains data on aerosols.
Subsets are extracted as *.csv* files.

## tavg1\_2d\_rad\_NX\_M2T1NXRAD :
Notebooks and scripts to decode MERRA-2 data for the dataset **tavg1\_2d\_rad\_Nx\_M2T1NXRAD**
This dataset contains information on cloud coverage and attenuation.
Subsets are extracted as *.csv* files.

## tavg1\_2d\_csp\_Nx\_M2T1NXCSP :
Decode old dataset **tavg1\_2d\_csp\_Nx\_M2T1NXCSP** on clouds. It is not used for CTIO 2017 data.
Subsets are extracted as *.csv* files.

## MergeDataSets
The above previous tools produced subsets are extracted as *.csv* files indexed by time.
These scripts/notebooks merge the data from the different dataset in a single dataset file in *.csv*.

## atmtranspsim
Does the libratran simulation of atmospheric transparency from the merged dataset and produce a results in a single fits file. 

## atmanimation
Does the animation of the atmospheric transparency profile from the *.fits* file produced previously 

## dataacess	
Memo on how to get the MERRA-2 "*.nc4*" files.		 

 
import numpy as np

BaseURL="https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2"
DataSet="M2T1NXAER.5.12.4"
Year="2016"

DataSetName_root="MERRA2_400.tavg1_2d_aer_Nx.20180801.nc4"
DataSetName="MERRA2_400.tavg1_2d_aer_Nx."

All_Months=["01","02","03","04","05","06","07","08","09","10","11","12"]
All_Days=["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"]

for month in All_Months:
	for day in All_Days:
		url=BaseURL+'/'+DataSet+'/'+Year+'/'+month+'/'+DataSetName+Year+month+day+".nc4"
		print url


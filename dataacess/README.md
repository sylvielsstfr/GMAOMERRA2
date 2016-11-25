
#Top page at GES DISC
=====================

http://disc.sci.gsfc.nasa.gov/uui/datasets?keywords=%22MERRA-2%22


product for my test:

inst1_2d_asm_Nx_M2I1NXASM

0.5 degree x 0.625 degree
hourly


## Top page for the product
---------------------------
http://disc.sci.gsfc.nasa.gov/uui/datasets/M2I1NXASM_V5.12.4/summary?keywords=%22MERRA-2%22


http://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/


http://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/2016/10/

194 MB/per file

 
## List files

wget -q -nH -nd http://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/2016/10/ | cut -f4 -d\"

curl -s http://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/2016/10/ | grep MERRA2_400 | cut -f4 -d\"



wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies http://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/2016/10/


## Load the file
wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies  http://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/2016/10//MERRA2_400.inst1_2d_asm_Nx.20161031.nc4
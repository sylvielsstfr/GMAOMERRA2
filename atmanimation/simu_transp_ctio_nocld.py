#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 17:34:54 2018

@author: dagoret
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import os
from astropy.io import fits
import numpy as np
import pandas as pd

from astropy.io import fits

#  column 0 : count number
#  column 1 : aerosol value
#  column 2 : pwv value
#  column 3 : ozone value
#  column 6 : data start 
#
index_atm_count=0
index_atm_aer=1
index_atm_pwv=2
index_atm_oz=3
index_atm_ps=4
index_atm_cloud=5
index_atm_data=6


file_simu='MERRA2_2017_M2I1NXASM_M2T1NXAER_M2T1NXRAD_ctio_atmsim.fits'
hdu = fits.open(file_simu)
data=hdu[0].data

file_merra2='MERRA2_2017_M2I1NXASM_M2T1NXAER_M2T1NXRAD_ctio_AllYear.csv'
df_merra2=pd.read_csv(file_merra2,index_col=0)
df_merra2.index.name='time'

file_logbook_ctio='ctiofulllogbook_jun2017_v4.csv'
df_ctio=pd.read_csv(file_logbook_ctio,sep=';')
df_ctio=df_ctio.reindex(columns=['date','P','T','RH','airmass','seeing','exposure','object','filter','disperser','focus','W','file']).set_index('date').sort_index()

transp=data[1:,index_atm_data:]
WL=data[0,index_atm_data:]


transp_sum=np.sum(transp,axis=1)
selected_indexes=np.where(transp_sum>0)[0]

NB_Frames=len(selected_indexes)

fig, ax = plt.subplots()
fig.set_tight_layout(True)
ax.grid(True)
ax.set_title('Atmospheric transparency CTIO Obs May-Jun 2017')
ax.set_ylabel('atm. transparency')

transp_curv, = ax.plot(WL, transp[selected_indexes[0],:], 'r-', linewidth=2)

# Query the figure's on-screen size and DPI. Note that when saving the figure to
# a file, we need to provide a DPI for that separately.
print('fig size: {0} DPI, size in inches {1}'.format(
    fig.get_dpi(), fig.get_size_inches()))


def update(i):
    
    sel_index=selected_indexes[i]
    label = 'wavelength (nm) '+str(i)+') date :: ',df_merra2.index.get_values()[sel_index]
    #print(label)
    tau_cloud=data[sel_index+1,index_atm_cloud]
    att_cloud=np.exp(-tau_cloud)
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    #transp_curv.set_ydata(att_cloud*transp[sel_index,:])
    transp_curv.set_ydata(transp[sel_index,:])
    ax.set_xlabel(label)
    return transp_curv, ax


if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(fig, update, frames=np.arange(0, NB_Frames), interval=200)
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        anim.save('sim.gif', dpi=80, writer='imagemagick')
    else:
        # plt.show() will just loop the animation forever.
        plt.show()

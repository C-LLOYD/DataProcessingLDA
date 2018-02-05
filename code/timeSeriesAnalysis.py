###########################################################################################
##############		TIME SERIES ANALYSIS SCRIPT		###############################
####
####		This script imports time series dataFrames and plots variables 
####
####
####
###########################################################################################
####
##		Initialise python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from movingAverageFilterReynoldsStresses import movingWeightedAverage as mwa
#
###########################################################################################
timeSeries = '.000003.pkl'

SA0 		= pd.read_pickle(str('../data/processedData/dataQualityTests/8Hz/angle0/seeding/timeSeries/raw/raw_8hz_seeding_angle0'+ timeSeries))
NSA0	 	= pd.read_pickle(str('../data/processedData/dataQualityTests/8Hz/angle0/noSeeding/timeSeries/raw/raw_8hz_noSeeding_angle0'+ timeSeries))
SA27 		= pd.read_pickle(str('../data/processedData/dataQualityTests/8Hz/angle2p7/seeding/timeSeries/raw/raw_8hz_seeding_angle2p7'+ timeSeries))
NSA27 	= pd.read_pickle(str('../data/processedData/dataQualityTests/8Hz/angle2p7/noSeeding/timeSeries/raw/raw_8hz_noSeeding_angle2p7'+ timeSeries))

#SA0 		= pd.read_pickle(str('../data/processedData/dataQualityTests/8Hz/angle0/seeding/timeSeries/w50_MA_high/w50_MA_high_8hz_seeding_angle0'+ timeSeries))
#NSA0	 	= pd.read_pickle(str('../data/processedData/dataQualityTests/8Hz/angle0/noSeeding/timeSeries/w50_MA_high/w50_MA_high_8hz_noSeeding_angle0'+ timeSeries))
#SA27 		= pd.read_pickle(str('../data/processedData/dataQualityTests/8Hz/angle2p7/seeding/timeSeries/w50_MA_high/w50_MA_high_8hz_seeding_angle2p7'+ timeSeries))
#NSA27 	= pd.read_pickle(str('../data/processedData/dataQualityTests/8Hz/angle2p7/noSeeding/timeSeries/w50_MA_high/w50_MA_high_8hz_noSeeding_angle2p7'+ timeSeries))

#../data/processedData/dataQualityTests/8Hz/angle0/noSeeding/timeSeries/raw/raw_8hz_noSeeding_angle0.000001.pkl
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'16'})
plt.rc('text', usetex=True)
fig = plt.figure()
ax1 = fig.add_subplot(4,1,1)
plt.plot(SA0.timeStamp,SA0.Ux,'xk')
plt.minorticks_on()
plt.grid(True, which='minor',alpha=0.6)
plt.grid(True, which='major',linewidth=0.9)
plt.xlim(0,40)
#plt.ylim(-0.05,0.05)
ax1.locator_params(nbins=4, axis='y')
cur_axes = plt.gca()
cur_axes.axes.xaxis.set_ticklabels([])
#ax.set_xlabel(xLabel,fontsize='30')
#plt.ylabel(r'$U_y$ (m/s)',fontsize='30')
plt.tight_layout()
###########################
ax2 = fig.add_subplot(4,1,2)
plt.plot(NSA0.timeStamp,NSA0.Ux,'xk')
plt.minorticks_on()
plt.grid(True, which='minor',alpha=0.6)
plt.grid(True, which='major',linewidth=0.9)
plt.xlim(0,40)
#plt.ylim(-0.05,0.05)
ax2.locator_params(nbins=4, axis='y')
cur_axes = plt.gca()
cur_axes.axes.xaxis.set_ticklabels([])
#ax.set_xlabel(xLabel,fontsize='30')
plt.ylabel(r'$U_x$ (m/s)',fontsize='20')
plt.tight_layout()
###########################
ax3 = fig.add_subplot(4,1,3)
plt.plot(SA27.timeStamp,SA27.Ux,'xk')
#plt.plot(fil.timeStamp,UyMean)
plt.minorticks_on()
plt.grid(True, which='minor',alpha=0.6)
plt.grid(True, which='major',linewidth=0.9)
plt.xlim(0,40)
#plt.ylim(-0.05,0.05)
ax3.locator_params(nbins=4, axis='y')
cur_axes = plt.gca()
cur_axes.axes.xaxis.set_ticklabels([])
#plt.ylabel(r'$U_y$ (m/s)',fontsize='30')
plt.tight_layout()
###########################
ax3 = fig.add_subplot(4,1,4)
plt.plot(NSA27.timeStamp,NSA27.Ux,'xk')
#plt.plot(fil.timeStamp,UyMean)
plt.minorticks_on()
plt.grid(True, which='minor',alpha=0.6)
plt.grid(True, which='major',linewidth=0.9)
plt.xlim(0,40)
#plt.ylim(-0.05,0.05)
ax3.locator_params(nbins=4, axis='y')
plt.xlabel(r'$t$ (s)',fontsize='20')
#plt.ylabel(r'$U_y$ (m/s)',fontsize='30')
plt.tight_layout()
plt.subplots_adjust(wspace=0.15, hspace=0.15)
#plt.show()
plt.savefig('../data/processedData/figures/dataQualityTests/Ux_timeSeries.png')
plt.close()


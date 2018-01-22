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
dec 		= pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/timeSeries/raw/raw_4Hz_x400.000041.pkl')
decFil 	= pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/timeSeries/high/high_4Hz_x400.000041.pkl')
mar 		= pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/1703XX_4Hz_x400/timeSeries/raw/raw_4Hz_x400.000142.pkl')
#fil = pd.read_pickle('../data/processedData/smoothPlate/16Hz/x400/170816_16Hz_x400/timeSeries/raw/raw_16Hz_x400.000004.pkl')
#fil = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/timeSeries/raw/raw_4Hz_x400.000014.pkl')
#raw = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/timeSeries/raw/raw_4Hz_x400.000049.pkl')
#print(fil,raw)
DecAv = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_raw.pkl')
MarAv = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/1703XX_4Hz_x400/4Hz_x400_averaged_raw.pkl')

print(DecAv)
#print(MarAv)
#print(MarAv.loc[MarAv['z']<11])
#UyMean = mwa(25,decFil.Uy,decFil.resTime)
mar = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/1703XX_4Hz_x400/timeSeries/raw/raw_4Hz_x400.000022.pkl')
dec = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/timeSeries/raw/raw_4Hz_x400.000014.pkl')
#print(temp)
#print(temp2)
#plt.plot(temp.timeStamp,temp.Uy,'-x')
#plt.plot(temp2.timeStamp,temp2.Uy,'-x')
#plt.show()


plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'16'})
plt.rc('text', usetex=True)
fig = plt.figure()
ax1 = fig.add_subplot(3,1,1)
plt.plot(dec.timeStamp,dec.Uy,'xk')
plt.minorticks_on()
plt.grid(True, which='minor',alpha=0.6)
plt.grid(True, which='major',linewidth=0.9)
plt.xlim(5,10)
plt.ylim(-0.2,0.2)
ax1.locator_params(nbins=4, axis='y')
cur_axes = plt.gca()
cur_axes.axes.xaxis.set_ticklabels([])
#ax.set_xlabel(xLabel,fontsize='30')
#plt.ylabel(r'$U_y$ (m/s)',fontsize='30')
plt.tight_layout()
###########################
ax2 = fig.add_subplot(3,1,2)
plt.plot(dec.timeStamp,dec.Uy,'xk')
plt.minorticks_on()
plt.grid(True, which='minor',alpha=0.6)
plt.grid(True, which='major',linewidth=0.9)
plt.xlim(5,10)
plt.ylim(-0.004,0.004)
ax2.locator_params(nbins=4, axis='y')
cur_axes = plt.gca()
cur_axes.axes.xaxis.set_ticklabels([])
#ax.set_xlabel(xLabel,fontsize='30')
plt.ylabel(r'$U_y$ (m/s)',fontsize='20')
plt.tight_layout()
###########################
ax3 = fig.add_subplot(3,1,3)
plt.plot(mar.timeStamp,mar.Uy,'xk')
#plt.plot(fil.timeStamp,UyMean)
plt.minorticks_on()
plt.grid(True, which='minor',alpha=0.6)
plt.grid(True, which='major',linewidth=0.9)
plt.xlim(5,10)
plt.ylim(-0.004,0.004)
ax3.locator_params(nbins=4, axis='y')
plt.xlabel(r'$t$ (s)',fontsize='20')
#plt.ylabel(r'$U_y$ (m/s)',fontsize='30')
plt.tight_layout()
plt.subplots_adjust(wspace=0.15, hspace=0.15)
#plt.show()
plt.savefig('../data/processedData/figures/filterDependence/timeSeriesAnalysis_z=0p5.png')
plt.close()


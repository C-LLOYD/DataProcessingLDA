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
###
###	timeSeriesPlotter function:	Reads in dataframe and savename, plots 10 seconds of Uy and Ux data, and saves it.
def setPlotParams(axes):
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
	plt.xlim(0,40)
#	plt.ylim(-0.05,0.05)
	axes.locator_params(nbins=4, axis='y')
	
def timeSeriesPlotter(df,saveName):
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'16'})
	plt.rc('text', usetex=True)
	fig = plt.figure()
	ax1 = fig.add_subplot(2,1,1)
	plt.plot(df.timeStamp,df.Ux,'xk')
	setPlotParams(ax1)
	cur_axes = plt.gca()
	cur_axes.axes.xaxis.set_ticklabels([])
#	ax.set_xlabel(r'$t$ (s)',fontsize='30')
	plt.ylabel(r'$U_1$ (m/s)',fontsize='30')
	plt.tight_layout()
###########################
	ax2 = fig.add_subplot(2,1,2)
	plt.plot(df.timeStamp,df.Uy,'xk')
	setPlotParams(ax2)
#plt.ylim(-0.05,0.05)
	ax2.set_xlabel(r'$t$ (s)',fontsize='30')
	plt.ylabel(r'$U_2$ (m/s)',fontsize='30')
	plt.tight_layout()
	plt.subplots_adjust(wspace=0.15, hspace=0.15)
#	plt.show()
	plt.savefig(saveName)
	plt.close()		



timeSeries = [
			'.000001.pkl',
			'.000002.pkl',
			'.000003.pkl',
		]		

fileNamePre = [
			'../data/processedData/dataQualityTests/8Hz/rotation0/angle0/timeSeries/raw/raw_8hz_rotation0_nodip',
			'../data/processedData/dataQualityTests/8Hz/rotation0/angle2p7/timeSeries/raw/raw_8hz_rotation0_dip',
			'../data/processedData/dataQualityTests/8Hz/rotation90/angle0/timeSeries/raw/raw_8hz_rotation90_nodip',
			'../data/processedData/dataQualityTests/8Hz/rotation90/angle2p7/timeSeries/raw/raw_8hz_rotation90_dip',
			'../data/processedData/dataQualityTests/8Hz/rotation45/angle0/timeSeries/raw/raw_8hz_rotation45_nodip',
			'../data/processedData/dataQualityTests/8Hz/rotation45/angle2p7/timeSeries/raw/raw_8hz_rotation45_dip',
		]

saveNames = [
			'../data/processedData/figures/dataQualityTests/rotation0_nodip.png',
			'../data/processedData/figures/dataQualityTests/rotation0_dip.png',
			'../data/processedData/figures/dataQualityTests/rotation45_nodip.png',
			'../data/processedData/figures/dataQualityTests/rotation45_dip.png',
			'../data/processedData/figures/dataQualityTests/rotation90_nodip.png',
			'../data/processedData/figures/dataQualityTests/rotation90_dip.png',
		]

for f in range(len(fileNamePre)):
	for strEnd in range(len(timeSeries)):
		fileName = str(fileNamePre[f] + timeSeries[strEnd])
		print(fileName)
		data = pd.read_pickle(fileName)
		timeSeriesPlotter(data,saveNames[f])



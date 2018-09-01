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
#	plt.xlim(0,40)
#	plt.ylim(-0.05,0.05)
	axes.locator_params(nbins=4, axis='y')
	
def timeSeriesPlotter(df,dfMean,TS,saveName):
	uPrime = df.Ux - dfMean.UxMean.loc[dfMean.fileName == str(TS + '.txt')].as_matrix()
	vPrime = df.Uy - dfMean.UyMean.loc[dfMean.fileName == str(TS + '.txt')].as_matrix()
	uvNorm = uPrime*vPrime/(dfMean.uxRMS.loc[dfMean.fileName == str(TS + '.txt')].as_matrix()
						*dfMean.uyRMS.loc[dfMean.fileName == str(TS + '.txt')].as_matrix()
					)
#	print(uPrime)
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'16'})
	plt.rc('text', usetex=True)
	fig = plt.figure()
#	ax1 = fig.add_subplot(2,1,1)
	plt.plot(df.timeStamp,uvNorm,'-xk')
#	setPlotParams(ax1)
#	cur_axes = plt.gca()
#	cur_axes.axes.xaxis.set_ticklabels([])
#	ax.set_xlabel(r'$t$ (s)',fontsize='30')
#	plt.ylabel(r'$U_1$ (m/s)',fontsize='30')
#	plt.tight_layout()
###########################
#	ax2 = fig.add_subplot(2,1,2)
#	plt.plot(df.timeStamp,vPrime,'-xk')
#	setPlotParams(ax2)
#plt.ylim(-0.05,0.05)
#	ax2.set_xlabel(r'$t$ (s)',fontsize='30')
#	plt.ylabel(r'$U_2$ (m/s)',fontsize='30')
#	plt.tight_layout()
#	plt.subplots_adjust(wspace=0.15, hspace=0.15)
	plt.show()
#	plt.savefig(saveName)
	plt.close()		



timeSeries = [
##			'.000001.pkl',
##			'.000002.pkl',
			'4Hz_x400.000010',
		]		

fileNamePre = [
#	'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/timeSeries/lowFil_0_rotation/lowFil_0_rotation_',
	'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/timeSeries/lowFil/lowFil_',
#	'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/timeSeries/raw/raw_',
#	'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/timeSeries/raw/raw_',
]

saveNames = [
	'../temp.png',
]

dfMean = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl')
print(dfMean)
#dfMean = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl')



for f in range(len(fileNamePre)):
	for ts in range(len(timeSeries)):
		fileName = str(fileNamePre[f] + timeSeries[ts] + '.pkl')
		data = pd.read_pickle(fileName)
		print(data)
		timeSeriesPlotter(data,dfMean,timeSeries[ts],saveNames[f])



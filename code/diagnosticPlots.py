###########################################################################################
##############		PROFILE ANALYSIS SCRIPT		#####################################
####
####		This script imports profile dataFrames and plots variables 
####
####
####
###########################################################################################
####
##		Initialise python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#
###########################################################################################
#
##	import data
fileNameFlatPlate  = [
##	2 min Av:
#		'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
##	5 min Av:
		'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
##	200 sec Av:
		'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',		
		]

fileNameRibletted  = [
##	2 min Av:
#		'../data/processedData/riblettedDenticles/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/riblettedDenticles/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/riblettedDenticles/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
##	5 min Av:
		'../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
##	200 sec Av:
		'../data/processedData/riblettedDenticles/8Hz/x400/180808_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/riblettedDenticles/12Hz/x400/180815_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/riblettedDenticles/16Hz/x400/180823_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
		]

fileNameSmooth = [
##	5 min Av:
		'../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
##	200 sec Av:
		'../data/processedData/smoothDenticles/8Hz/x400/180809_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothDenticles/12Hz/x400/180816_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothDenticles/16Hz/x400/180824_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
]

TSlabels =  ['flat plate',
		 'ribletted denticles',
		 'smooth denticles',
		]

markers = ['o','v','s','*','+','^','.','<','>','p','x',]

###########################################################################################################
###########			Plotter function
def plotter(X,Y,dataLabels,markers,xLabel,yLabel,writeName):
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	for i in range(len(X)):
#		ax.semilogx(X[i],Y[i],label = dataLabels[i], marker=markers[i], linestyle='')#, color = 'k')
		ax.plot(X[i],Y[i],label = dataLabels[i], marker=markers[i], linestyle='')#, color = 'k')
	ax.set_xlabel(xLabel,fontsize='30')
	ax.set_ylabel(yLabel,fontsize='30')#,rotation=0,labelpad=45)
	plt.tight_layout()
#	h.set_rotation(0)
	if legend == True:
		box = ax.get_position()
#		ax.set_position([box.x0, box.y0*1.2, box.width * 0.75, box.height])
		ax.legend(loc='upper left',numpoints=1,prop={'size':18})#, bbox_to_anchor=(1, 0.5),numpoints=1)	
#		plt.legend(loc = 'best')#handles=[plot,],loc=4,prop={'size':18},ncol=1)
#	if np.mean(Y[i]) > 0:
#		ax.set_ylim([0, 1.05*np.max(Y[i])])	
#	else:
#		ax.set_ylim([1.05*np.min(Y[i]),1.05*np.max(Y[i])])
	
	ax.set_xlim([0,1])		
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
#	writeName = str('../data/processedData/figures/diagnosticPlots/'++varNames[var]+'.png')
	plt.savefig(writeName)
#	plt.show()
	plt.close()
#
###########################################################################################################
##
##		Comparing 8Hz smooth plate profiles at diff samping rates
dataFlatPlates = [
	pd.read_pickle(fileNameFlatPlate[0]),
	pd.read_pickle(fileNameFlatPlate[1]),
	pd.read_pickle(fileNameFlatPlate[2]),
	pd.read_pickle(fileNameFlatPlate[3]),	
]

dataRib = [
	pd.read_pickle(fileNameRibletted[0]),
	pd.read_pickle(fileNameRibletted[1]),
	pd.read_pickle(fileNameRibletted[2]),
	pd.read_pickle(fileNameRibletted[3])
]

dataSmo = [
	pd.read_pickle(fileNameSmooth[0]),
	pd.read_pickle(fileNameSmooth[1]),
	pd.read_pickle(fileNameSmooth[2]),
	pd.read_pickle(fileNameSmooth[3])
]

var = [
	"uxRMS",
	"uyRMS",
	"uv"
]
power = [
	1,
	1,
	2,
]
ylabels = [r"$\sqrt{\overline{u'u'}}/U_\infty$",r"$\sqrt{\overline{v'v'}}/U_\infty$",r"$\overline{u'v'}/(U_\infty)^2$"]

writePath = str('../data/processedData/figures/diagnosticPlots/')
writeNamesPre = ['4Hz','8Hz','12Hz','16Hz']
writeNames = [
	str('uRMS.png'),
	str('vRMS.png'),
	str('uv.png')
]
print(dataFlatPlates[0])
legend = False
for d in range(len(dataFlatPlates)):
	for i in range(len(var)):
		plotter(
			[
			dataFlatPlates[d]["UxMean"].loc[dataFlatPlates[d]["fil"]]/np.mean(dataFlatPlates[d]["UxMean"][-4:]),
			dataRib[d]["UxMean"].loc[dataRib[d]["fil"]]/np.mean(dataRib[d]["UxMean"][-4:]),
			dataSmo[d]["UxMean"].loc[dataSmo[d]["fil"]]/np.mean(dataSmo[d]["UxMean"][-4:])
			],
			[
			dataFlatPlates[d][var[i]].loc[dataFlatPlates[d]["fil"]]/np.mean(dataFlatPlates[d]["UxMean"][-4:])**power[i],
			dataRib[d][var[i]].loc[dataRib[d]["fil"]]/np.mean(dataRib[d]["UxMean"][-4:])**power[i],
			dataSmo[d][var[i]].loc[dataSmo[d]["fil"]]/np.mean(dataSmo[d]["UxMean"][-4:])**power[i]
			],
			[
				'Flat',
				'Rib',
				'Smo'
			],
			[
				markers[0],
				markers[1],
				markers[2]
			],
			r'$\overline{U}/U_\infty$', ylabels[i], str(writePath + writeNamesPre[d] +  writeNames[i])
		)
		legend = False




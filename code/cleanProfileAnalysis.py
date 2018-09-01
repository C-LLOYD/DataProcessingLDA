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
		'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
##	5 min Av:
		'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
##	200 sec Av:
		'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',		
		]

fileNameRibletted  = [
##	2 min Av:
		'../data/processedData/riblettedDenticles/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/riblettedDenticles/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/riblettedDenticles/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
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

markers = ['x','o','v','s','*','+','^','.','<','>','p']

###########################################################################################################
###########			Plotter function
def plotter(X,Y,dataLabels,markers,xLabel,yLabel,writeName):
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	for i in range(len(X)):
		ax.semilogx(X[i],Y[i],label = dataLabels[i], marker=markers[i], linestyle='-')#, color = 'k')
#		ax.plot(X[i],Y[i],label = dataLabels[i], marker=markers[i], linestyle='-')#, color = 'k')
	ax.set_xlabel(xLabel,fontsize='30')
	ax.set_ylabel(yLabel,fontsize='30')#,rotation=0,labelpad=45)
	plt.tight_layout()
#	h.set_rotation(0)
	if legend == True:
		box = ax.get_position()
#		ax.set_position([box.x0, box.y0*1.2, box.width * 0.75, box.height])
		ax.legend(loc='upper left',numpoints=1,prop={'size':18})#, bbox_to_anchor=(1, 0.5),numpoints=1)	
#		plt.legend(loc = 'best')#handles=[plot,],loc=4,prop={'size':18},ncol=1)
#	ax.set_xlim([0.1, 180])
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
	#writeName = str('../plots/freeStream_'+varNames[var]+'.png')
#	plt.savefig(writeName)
	plt.show()
	plt.close()
#
###########################################################################################################
#print(pd.read_pickle(fileName8Hz[1]))

dataFlatPlate = [
## 2 min Av:
#	pd.read_pickle(fileNameFlatPlate[0]),	#.loc[pd.read_pickle(fileNameFlatPlate[0]).z>0.5],
#	pd.read_pickle(fileNameFlatPlate[1]),	#.loc[pd.read_pickle(fileNameFlatPlate[1]).z>0.46],
#	pd.read_pickle(fileNameFlatPlate[2]),	#.loc[pd.read_pickle(fileNameFlatPlate[2]).z>1.12],
## 5 min Av:
	pd.read_pickle(fileNameFlatPlate[3]),
## 200 sec Av:
#	pd.read_pickle(fileNameFlatPlate[4]),
]

dataRibletted = [	
## 2 min Av:
#	pd.read_pickle(fileNameRibletted[0]),#.loc[pd.read_pickle(fileNameRibletted[0]).fMean<27],
#	pd.read_pickle(fileNameRibletted[1]),#.loc[(pd.read_pickle(fileNameRibletted[1]).fMean<30) & (pd.read_pickle(fileNameRibletted[1]).z>2.0)],
#	pd.read_pickle(fileNameRibletted[2]),#.loc[pd.read_pickle(fileNameRibletted[2]).z>2.26],
## 5 min Av:
	pd.read_pickle(fileNameRibletted[3]),
]

dataSmooth = [
## 5 min Av:
	pd.read_pickle(fileNameSmooth[0]),
]


#print(pd.read_pickle(fileName8Hz[1]).loc[pd.read_pickle(fileName8Hz[1]).z>2.15])

#
#
#var = ["fMean"]
#ylabels = ['fMean']
var = ["UxMean","uxRMS","uyRMS","uv"]
ylabels = [r'$\overline{U}$ (m/s)',r"$\sqrt{\overline{u'u'}}$ (m/s)",r"$\sqrt{\overline{v'v'}}$ (m/s)",r"$\overline{u'v'}$, (m\textsuperscript{2}/s\textsuperscript{2})"]

writePath = str('../data/processedData/figures/1803XX_profiles/')
writeNamesPre = ['4Hz']#,'8Hz','16Hz']
writeNames = [
	str('_Uprofile_lowFil.png'),
	str('_uRMSprofile_lowFil.png'),
	str('_vRMSprofile_lowFil.png'),
	str('_uvProfile_lowFil.png')
]

offsetRib = [0]
offsetSmo = [0]
temp = False
if temp == True:
	for k in range(len(dataFlatPlate)):
		legend = True
#	print(str(writeNamesPre[k] + ' Flat Plate'))
#	print(dataFlatPlate[k])
#	print(str(writeNamesPre[k] + ' Ribletted Denticles'))
#	print(dataRibletted[k])
		for i in range(len(var)):
			for j in range(len(offsetRib)):
				plotter(
#					[	dataRibletted[0]["z"],
#						dataRibletted[1]["z"],
					[	dataFlatPlate[k]["z"],
						dataRibletted[k]["z"] - offsetRib[j],
						dataSmooth[k]["z"] - offsetSmo[j],

					],
#					[	dataRibletted[0][var[i]],
#						dataRibletted[1][var[i]],
					[	dataFlatPlate[k][var[i]],
						dataRibletted[k][var[i]],
						dataSmooth[k][var[i]],
					],
					['Flat Plate','Ribletted Denticles', 'Smooth Denticles'],
					[markers[0],markers[1],markers[2]],
					r'$y$ (mm)', ylabels[i], str(writePath + writeNamesPre[k] + writeNames[i])
				)
				legend = False

###########################################################################################################
##
##		Comparing 8Hz smooth plate profiles at diff samping rates
dataFlatPlates = [
#	pd.read_pickle(fileNameFlatPlate[3]),
#	pd.read_pickle(fileNameFlatPlate[4]),
#	pd.read_pickle(fileNameFlatPlate[5]),
	pd.read_pickle(fileNameFlatPlate[2]),
	pd.read_pickle(fileNameFlatPlate[6]),	
]

dataRib = [
	pd.read_pickle(fileNameRibletted[3]),
	pd.read_pickle(fileNameRibletted[4]),
	pd.read_pickle(fileNameRibletted[5]),
	pd.read_pickle(fileNameRibletted[6])
]

dataSmo = [
	pd.read_pickle(fileNameSmooth[0]),
	pd.read_pickle(fileNameSmooth[1]),
	pd.read_pickle(fileNameSmooth[2]),
	pd.read_pickle(fileNameSmooth[3])
]

writePath = str('../data/processedData/figures/180807_profiles/')
writeNamesPre = ['4Hz']#,'8Hz','16Hz']
writeNames = [
	str('Uprofile_lowFil.png'),
	str('uRMSprofile_lowFil.png'),
	str('vRMSprofile_lowFil.png'),
	str('uvProfile_lowFil.png')
]
legend = True
for i in range(len(var)):
	plotter(
#				[	dataRibletted[0]["z"],
#					dataRibletted[1]["z"],
		[
#			dataFlatPlates[0]["z"]*dataFlatPlates[0]["fil"],
#			dataRib[0]["z"],
#			dataSmo[0]["z"],
#			dataFlatPlates[0]["z"],
			dataFlatPlates[1]["z"],
			dataRib[3]["z"],
			dataSmo[3]["z"]
		],
		[
#			dataFlatPlates[0][var[i]]*dataFlatPlates[0]["fil"],
#			dataRib[0][var[i]],
#			dataSmo[0][var[i]],
#			dataFlatPlates[0][var[i]],
			dataFlatPlates[1][var[i]],
			dataRib[3][var[i]],
			dataSmo[3][var[i]]
		],
		[
#			'Rib',
#			'Smo',
			'Flat',
			'Rib',
			'Smo'
		],
		[
			markers[0],
			markers[1],
			markers[2]
		],
		r'$y$ (mm)', ylabels[i], str(writePath + writeNames[i])
	)
	legend = False



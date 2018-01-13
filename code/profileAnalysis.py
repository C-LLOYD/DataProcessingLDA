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
testProfiles = True
findMeanSpikeFrac = True
#
##	import data
fileName4Hz = ['../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_raw.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_basicMin.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_basicMed.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_min.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_low.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_med.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_high.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_w50_MA_min.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_w50_MA_low.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_w50_MA_med.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_w50_MA_high.pkl']

fileName8Hz = ['../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_raw.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_basicMin.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_basicMed.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_min.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_low.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_med.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_high.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_w50_MA_min.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_w50_MA_low.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_w50_MA_med.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/8Hz_x400_averaged_w50_MA_high.pkl']

nameEnd = ['4Hz','8Hz']

###########################################################################################################
###########			Plotter function
def plotter(X,Y,dataLabels,markers,xLabel,yLabel,writeName):
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	for i in range(len(X)):
		ax.semilogx(X[i],Y[i],label = dataLabels[i], marker=markers[i], linestyle='-')#, color = 'k')
	ax.set_xlabel(xLabel,fontsize='30')
	ax.set_ylabel(yLabel,fontsize='30')#,rotation=0,labelpad=45)
	plt.tight_layout()
#	h.set_rotation(0)
	if legend == True:
		box = ax.get_position()
#		ax.set_position([box.x0, box.y0*1.2, box.width * 0.75, box.height])
		ax.legend(loc='center left',numpoints=1,prop={'size':18})#, bbox_to_anchor=(1, 0.5),numpoints=1)	
#		plt.legend(loc = 'best')#handles=[plot,],loc=4,prop={'size':18},ncol=1)
#	ax.set_xlim([0.1, 180])
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
	#writeName = str('../plots/freeStream_'+varNames[var]+'.png')
	plt.savefig(writeName)
#	plt.show()
	plt.close()
#
###########################################################################################################

if findMeanSpikeFrac == True:
	for j in range(len(nameEnd)):
		if nameEnd[j] == '4Hz':
			fileName = fileName4Hz
		elif nameEnd[j] == '8Hz':
			fileName = fileName8Hz
		data = [	pd.read_pickle(fileName[1]),
				pd.read_pickle(fileName[2]),
				pd.read_pickle(fileName[3]),
				pd.read_pickle(fileName[4]),
				pd.read_pickle(fileName[5]),
				pd.read_pickle(fileName[6]),
				pd.read_pickle(fileName[7]),
				pd.read_pickle(fileName[8]),
				pd.read_pickle(fileName[9]),
				pd.read_pickle(fileName[10])
			  ]
		meanSpikeFrac = []
		for i in range(len(data)):
			meanSpikeFrac.append(np.mean(data[i]["spikeFrac"].loc[data[i].z > 1]))
		print(meanSpikeFrac)

if testProfiles == True:
	for j in range(len(nameEnd)):
		if nameEnd[j] == '4Hz':
			fileName = fileName4Hz
		elif nameEnd[j] == '8Hz':
			fileName = fileName8Hz
#
		data = [	pd.read_pickle(fileName[0]),
				pd.read_pickle(fileName[1]),
				pd.read_pickle(fileName[2]),
				pd.read_pickle(fileName[3]),
				pd.read_pickle(fileName[4]),
				pd.read_pickle(fileName[5]),
				pd.read_pickle(fileName[6]),
				pd.read_pickle(fileName[7]),
				pd.read_pickle(fileName[8]),
				pd.read_pickle(fileName[9]),
				pd.read_pickle(fileName[10])
			  ]
#
		markers = ['x','o','v','s','*','+','^','.','<','>','p']
		labels =  ['raw','MA, min', 'MA, med', 'MARS, min', 'MARS, low', 'MARS, med', 'MARS, high', 'MA W50, min', 'MA W50, low', 'MA W50, med', 'MA W50, high']
#
##	1.	test effect of filtering method
##			-	raw
##			-	basic, minimal filtering
##			-	RSM,	 minimal filtering
##			-	basic, medium filtering
##			-	RSM,	 medium filtering
		var = ["UxMean","uxRMS","uyRMS","uv"]
		ylabels = [r'$\mu_u$ (m/s)',r'$\sigma_u$ (m/s)',r'$\sigma_v$ (m/s)',r'$\gamma_{uv}$, (m\textsuperscript{2}/s\textsuperscript{2})']
		writeNames = [	str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterMethod_Uprofile.png'),
					str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterMethod_uRMSprofile.png'),
					str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterMethod_vRMSprofile.png'),
					str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterMethod_uvProfile.png')
				]
		legend = True
		for i in range(len(var)):
			plotter(	[data[0]["z"],data[1]["z"],data[2]["z"],data[3]["z"],data[5]["z"]], [data[0][var[i]],data[1][var[i]],data[2][var[i]],data[3][var[i]],data[5][var[i]]],
					[labels[0],labels[1],labels[2],labels[3],labels[5]], [markers[0],markers[1],markers[2],markers[3],markers[5]],
					r'$y$ (mm)', ylabels[i],writeNames[i])
			legend = False
#
##	1.	test effect of filtering strength
##			-	raw
##			-	RSM,	 minimal filtering
##			-	RSM,	 low filtering
##			-	RSM,	 medium filtering
##			-	RSM,	 high filtering
		writeNames = [	str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterStrength_Uprofile.png'),
					str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterStrength_uRMSprofile.png'),
					str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterStrength_vRMSprofile.png'),
					str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterStrength_uvProfile.png')
				]
		legend = True
		for i in range(len(var)):
			plotter(	[data[0]["z"],data[3]["z"],data[4]["z"],data[5]["z"],data[6]["z"]], [data[0][var[i]],data[3][var[i]],data[4][var[i]],data[5][var[i]],data[6][var[i]]],
					[labels[0],labels[3],labels[4],labels[5],labels[6]], [markers[0],markers[3],markers[4],markers[5],markers[6]],
					r'$y$ (mm)', ylabels[i],writeNames[i])
			legend = False
##	1.	test effect of filtering window
##			-	raw
##			-	MA_w50,	 minimal filtering
##			-	MA_w50,	 low filtering
##			-	MA_w50,	 medium filtering
##			-	MA_w50,	 high filtering
		writeNames = [	str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterWindow_Uprofile.png'),
					str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterWindow_uRMSprofile.png'),
					str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterWindow_vRMSprofile.png'),
					str('../data/processedData/figures/filterDependence/' + nameEnd[j] + 'filterWindow_uvProfile.png')
				]
		print(writeNames)
		legend = True
		for i in range(len(var)):
			plotter(	[data[0]["z"],data[7]["z"],data[8]["z"],data[9]["z"],data[10]["z"]], [data[0][var[i]],data[7][var[i]],data[8][var[i]],data[9][var[i]],data[10][var[i]]],
					[labels[0],labels[7],labels[8],labels[9],labels[10]], [markers[0],markers[7],markers[8],markers[9],markers[10]],
					r'$y$ (mm)', ylabels[i],writeNames[i])
#			legend = False


#
###########################################################################################################

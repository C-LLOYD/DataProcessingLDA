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

fileNameRibletted  = [
		'../data/processedData/riblettedDenticles/4Hz/x400/180918_4Hz_spatialDep/4Hz_spatialDep_averaged_lowFil.pkl',
		'../data/processedData/riblettedDenticles/12Hz/x400/180919_12Hz_spatialDep/12Hz_spatialDep_averaged_lowFil.pkl',
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
	plt.savefig(writeName)
#	plt.show()
	plt.close()
#
###########################################################################################################
#print(pd.read_pickle(fileName8Hz[1]))

data = [
	pd.read_pickle(fileNameRibletted[0]),
	pd.read_pickle(fileNameRibletted[1]),
]

writeName = [
	'../data/processedData/figures/spatialDependence/4Hz_',
	'../data/processedData/figures/spatialDependence/12Hz_'
]

loc_x1 = [400.0,400.0,401.0]
loc_x2 = [0.0, 1.0, 0.0]

legend = True
for d in range(len(data)):
	d1 = data[d].loc[(data[d]["x1"]==400.0) & (data[d]["x2"]==0.0)]
	d2 = data[d].loc[(data[d]["x1"]==400.0) & (data[d]["x2"]==1.0)]
	d3 = data[d].loc[(data[d]["x1"]==401.0) & (data[d]["x2"]==0.0)]
	var = ["UxMean","uxRMS","uyRMS","uv"]
	ylabels = [
		r'$\overline{U}$ (m/s)',
		r"$\sqrt{\overline{u'u'}}$ (m/s)",
		r"$\sqrt{\overline{v'v'}}$ (m/s)",
		r"$\overline{u'v'}$, (m\textsuperscript{2}/s\textsuperscript{2})"
	]
	for i in range(len(var)):
		twriteName = str(writeName[d] + var[i] + '.png')
		plotter(
			[d1["z"], d2["z"], d3["z"]],
			[d1[var[i]],d2[var[i]],d3[var[i]]],
			["loc1","loc2","loc3"],
			[markers[0],markers[1],markers[2]],
			r'$y$ (mm)', ylabels[i], twriteName
			
		)
		legend = False





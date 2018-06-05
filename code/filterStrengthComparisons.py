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
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#
##	Script separates data sets by flow rate and plate type
flatPlate16HzFileNames = [
	'../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_raw.pkl',
	'../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
	'../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_w50_MA_high.pkl',
	'../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_w200_MA_high.pkl',
]

flatPlate16HzDataNames = [
	'raw',
	'low fil',
	'w = 50, high fil',
	'w = 200, high fil',
]

flatPlate8HzFileNames = [
	'../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_raw.pkl',
	'../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
	'../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_w50_MA_high.pkl',
	'../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_w200_MA_high.pkl',
]

flatPlate8HzDataNames = [
	'raw',
	'low fil',
	'w = 50, high fil',
	'w = 200, high fil',
]

flatPlate4HzFileNames = [
	'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_raw.pkl',
	'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
	'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_w50_MA_high.pkl',
	'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_w200_MA_high.pkl',
]

flatPlate4HzDataNames = [
	'raw',
	'low fil',
	'w = 50, high fil',
	'w = 200, high fil',
]

markers = ['x','^','s','*']

###########################################################################################################
var = ["UxMean","uxRMS","uyRMS","uv"]
ylabels = [r'$U$ (m/s)',r'$\sqrt{<u'u'>}$ (m/s)',r"$\sqrt{<v'v'>}$ (m/s)",r"$<u'v'>$, (m\textsuperscript{2}/s\textsuperscript{2})"]

data = [0]*4
for k in range(len(flatPlate4HzFileNames)):
	data[k] = pd.read_pickle(flatPlate16HzFileNames[k])
	print(data[k])

legend = True
for i in range(len(var)):
	plotter(
		[data[0]["z"],data[1]["z"],data[2]["z"],data[3]["z"]],
		[data[0][var[i]],data[1][var[i]],data[2][var[i]],data[3][var[i]]],
		[flatPlate16HzDataNames[0],flatPlate16HzDataNames[1],flatPlate16HzDataNames[2],flatPlate16HzDataNames[3]],
		[markers[0],markers[1],markers[2],markers[3]],
		r'$y$ (mm)', ylabels[i], 'noWrite'
	)
	legend = False








##

###########################################################################################
##############		RAW DATA PROCESSING SCRIPT	#####################################
####
####		This script processes raw data of Dec2017.
####
####		Inputs txt files from LDA, converts to dataFrames (and saves), filters data
####		(and saves), then averages and compiles all data into a summary dataFrame.
####
####
###########################################################################################
####
##		Initialise python
import numpy as np
import pandas as pd
import glob
from processingFunctions import txtToDataFrame
from processingFunctions import timeAverage
from processingFunctions import transform
from processingFunctions import spatialFilter
from FilterFunctions import Filter
from skinFrictionFunctions import viscousSublayerEstimation as visFunc
from skinFrictionFunctions import logLawEstimation as logFunc
from skinFrictionFunctions import clauserEstimation as clauserFunc
from skinFrictionFunctions import choiEstimation as choiFunc
#
##########################################################################################
##
##	1.	Loop over data: Put this in at the end
#
##	2.	import txt file: Give a hard coded name for now
writeData = False
saveFil = False
writeSpikeFrac = True
writeSpatialFilter = False
writeSkinFrictionEstimations = False
writeDenticleSkinFriction = False
#
rawPathNames = [
#	'../data/rawData/smoothPlate/16Hz/x400/1703XX_16Hz_x400/*.txt',
	'../data/rawData/smoothPlate/4Hz/x400/180530_4Hz_x400/*.txt',
	'../data/rawData/smoothDenticles/4Hz/x400/180601_4Hz_x400/*.txt',
	'../data/rawData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/*.txt',
]
#
#"../data/rawData/smoothPlate/4Hz/x400/1703XX_4Hz_x400/*.txt"
#"../data/rawData/smoothPlate/16Hz/x400/170816_16Hz_x400/*.txt",
#		"../data/rawData/smoothPlate/8Hz/x400/171214_8Hz_x400/*.txt"
filterType = [	#'movingAverageFilter',
			'movingAverageFilter',
			'movingAverageFilter',
			'movingAverageFilter',
			'movingAverageFilter',
			'movingAverageFilter',
			'movingAverageFilter',
#			'movingAverageFilterReynoldsStresses',
#			'movingAverageFilterReynoldsStresses',
#			'movingAverageFilterReynoldsStresses',
#			'movingAverageFilterReynoldsStresses',
#			'movingAverageFilterReynoldsStresses'
		]
probeRotationAngle = [45.0,45.0]#,45.0,45.0,45.0,45.0,45.0]


filLoops = 	[0,1]#,0,0,0,0,0]#[1,2,1,2]# 0, 1, 2, 1, 2]
NstdDev = 	[0,2]#,0,0,0,0,0]#[4,4,2,2]# 0, 4, 4, 2, 2]
avWindow = [50,50]#,50,50,50,50,50]#,200]#,50]#[50, 50, 50, 50]

saveNames = ['raw','lowFil']#'w50_MA_high','w200_MA_high']#['raw','lowFil']#,'w50_MA_high']#['w50_MA_min', 'w50_MA_low','w50_MA_med','w50_MA_high']#['basicMin','basicMed']#'raw','min','low','med','high']
#
##	First loop converts txt files to data frames, then filters the time series data,
##	then compiles all data sets into an averaged data frame
##
if writeData == True:
	for j in range(len(rawPathNames)):
		for i in range(len(filLoops)):
			rawPath = rawPathNames[j]
			data = []
			print(rawPath)
			for fileName in glob.glob(rawPath):
				print(fileName)				
				tData = txtToDataFrame(fileName)
#
##	3.	filter data
#					print(rawSaveName)
				if isinstance(tData,pd.DataFrame):
					ttData = Filter(tData,filterType[i],'mean',avWindow[i],'none',filLoops[i],NstdDev[i])
					tttData = transform(ttData,probeRotationAngle[i])
#					print(filSaveName)
					if saveFil == True:
						tempSavePath = fileName.split('/')
						filSaveNameEnd = str('/' + saveNames[i] + '/' + saveNames[i] + '_' + tempSavePath[-1].split('.txt')[0] + '.pkl')
						tempSavePath[-1] = 'timeSeries'	
						tempSavePath[tempSavePath.index('rawData')] = 'processedData'			
						filSaveName = str('/'.join(tempSavePath)+filSaveNameEnd)
						tttData.to_pickle(filSaveName)
#	
##	4.	apply averaging and append a final data series
					dataNew  = timeAverage(tttData)
					dataNew['fileName'] = fileName.split("/")[-1]
					if not isinstance(data,pd.DataFrame):
						data = dataNew
					else:
						data=data.append(dataNew)
			dataSorted = data.sort_values(by=['z'])
			print(dataSorted)
			tempSavePath = fileName.split('/')
			tempSavePath[-1] = str(tempSavePath[-1].split('.')[0] + '_averaged_' + saveNames[i] + '.pkl')
			tempSavePath[tempSavePath.index('rawData')] = 'processedData'
			globSaveName = '/'.join(tempSavePath)				
			print(globSaveName)
			dataSorted.to_pickle(globSaveName)
#
#
##################################################################################################################################
masterPath = [
	"../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/",
	"../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/",
	"../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/",
#	"../data/processedData/riblettedDenticles/16Hz/x400/180426_16Hz_x400/",
#	"../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/",
#	"../data/processedData/riblettedDenticles/8Hz/x400/180322_8Hz_x400/",
#	"../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/",
#	"../data/processedData/riblettedDenticles/4Hz/x400/180427_4Hz_x400/",
#	"../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/"
]
		
filterType = ['lowFil']#['w50_MA_min', 'w50_MA_low','w50_MA_med','w50_MA_high']#[	'basicMin','basicMed']#'min','low','med','high'	]

avDataFileName = [
	'4Hz_x400_averaged',
	'4Hz_x400_averaged',
	'4Hz_x400_averaged',
]

if writeSpikeFrac == True:
	for i in range(len(masterPath)):
		rawSeriesNames = glob.glob(str(masterPath[i] + 'timeSeries/raw/*'))
		for j in range(len(filterType)):
			AvData = pd.read_pickle(str(masterPath[i] + avDataFileName[i] + '_' + filterType[j] + '.pkl'))
			AvData["spikeFrac"] = ""
			AvData = AvData.sort_values(by="fileName",ascending=1)
			AvData = AvData.set_index("fileName")
			AvDataFileStruc = AvData.index[0].split('.')
			filterSeriesNames = rawSeriesNames[:]
			for k in range(len(filterSeriesNames)):
				filterSeriesNames[k] = rawSeriesNames[k].replace('raw',filterType[j])
				rawTimeSeries = pd.read_pickle(rawSeriesNames[k])
				filterTimeSeries = pd.read_pickle(filterSeriesNames[k])
				spikeFracEntry = 1-float(len(filterTimeSeries.Ux))/float(len(rawTimeSeries.Ux))
				AvLoc = AvDataFileStruc[:]
				AvLoc[1] = filterSeriesNames[k].split('/')[-1].split('.')[1]
				AvData.set_value(	'.'.join(AvLoc),	"spikeFrac",	spikeFracEntry)
			AvData = AvData.reset_index()
			print(AvData)
			AvData.to_pickle(str(masterPath[i] + avDataFileName[i] + '_' + filterType[j] + '.pkl'))
#
#
##################################################################################################################################
#
#
dataNames = [
	'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_',
	'../data/processedData/riblettedDenticles/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_',
	'../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_',
	'../data/processedData/riblettedDenticles/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_',
	'../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_',
	'../data/processedData/riblettedDenticles/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_',
]
nameEnd = [
	'raw.pkl',
	'lowFil.pkl',
	'w50_MA_high.pkl',
	'w200_MA_high.pkl',
]

#data = pd.read_pickle('../data/processedData/riblettedDenticles/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_raw.pkl')
#data = pd.read_pickle('../data/processedData/riblettedDenticles/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_raw.pkl')
#data = pd.read_pickle('../data/processedData/riblettedDenticles/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_raw.pkl')
data = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_raw.pkl')
if writeSpatialFilter == True:
	for i in range(len(dataNames)):
		for j in range(len(nameEnd)):
			data = pd.read_pickle(str(dataNames[i] + nameEnd[j]))
			newData = spatialFilter(data)
			newData.to_pickle(str(dataNames[i] + nameEnd[j]))

#print(newdata)

if writeSkinFrictionEstimations == True:
	flatPlateDataNames = [
		'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_raw.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_raw.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
	] 
#
##	Set constants
	kappa = 0.42
	nu = 1e-6
#
## 	Loop through file names, read in data, and save uTau estimates
	for d in range(len(flatPlateDataNames)):
		df = pd.read_pickle(flatPlateDataNames[d])
		U = df.UxMean[df.fil].as_matrix()
		yTilde = df.z[df.fil].as_matrix()
		[a,b,e] = visFunc(U,yTilde,nu)	#Estimate offset from wall and uTau from viscous layer
		uTauVis = np.sqrt(b*nu*1e3)		#1e3 due to y measured in mm
		y = (yTilde + a/b)*1e-3
		yConsistent = pd.Series(df.z.as_matrix() + a/b)
		[a,b,e] = logFunc(U,y,uTauVis,kappa,nu)	#get another uTau estimate from loglayer
		uTauLog = b*kappa		
		[uTauClauser, deltaStar, theta, H] = clauserFunc(U,y,kappa,nu)	#get another uTau estimate 
#		dfNew = df.drop(['uTau','theta','deltaStar','H','y'],axis=1)#.reset_index().copy()
		dfNew = df.copy()
		dfNew["y"], dfNew["uTau"], dfNew["deltaStar"], dfNew["theta"], dfNew['H'] = \
			np.nan, np.nan, np.nan, np.nan, np.nan
		dfNew["y"] = yConsistent
		dfNew["uTau"] = pd.Series([uTauVis,uTauLog,uTauClauser])
		dfNew["deltaStar"] = pd.Series([deltaStar])
		dfNew["theta"] = pd.Series([theta])
		dfNew["H"] = pd.Series([H])
		dfNew.to_pickle(flatPlateDataNames[d])

if writeDenticleSkinFriction == True:
	flatPlateData = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl')
	riblettedDenticleData = pd.read_pickle('../data/processedData/riblettedDenticles/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl')
	choiFunc(flatPlateData["UxMean"][flatPlateData["fil"]],
		flatPlateData["y"][flatPlateData["fil"]]*1e-3,	
		flatPlateData["uTau"][0],
		flatPlateData["deltaStar"][0],
		'',''	
		)















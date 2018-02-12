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
from FilterFunctions import Filter
#
##########################################################################################
##
##	1.	Loop over data: Put this in at the end
#
##	2.	import txt file: Give a hard coded name for now
writeData = True
saveFil = True
writeSpikeFrac = False
#
rawPathNames = 	[
				'../data/rawData/dataQualityTests/8Hz/profiles/rotation0/*.txt',
				'../data/rawData/dataQualityTests/8Hz/profiles/rotation45/*.txt',
			]
#
#"../data/rawData/smoothPlate/4Hz/x400/1703XX_4Hz_x400/*.txt"
#"../data/rawData/smoothPlate/16Hz/x400/170816_16Hz_x400/*.txt",
#		"../data/rawData/smoothPlate/8Hz/x400/171214_8Hz_x400/*.txt"
filterType = [	#'movingAverageFilter',
			'movingAverageFilter',
			'movingAverageFilter',
#			'movingAverageFilter',
#			'movingAverageFilterReynoldsStresses',
#			'movingAverageFilterReynoldsStresses',
#			'movingAverageFilterReynoldsStresses',
#			'movingAverageFilterReynoldsStresses',
#			'movingAverageFilterReynoldsStresses'
		]
filLoops = 	[0,2]#[1,2,1,2]# 0, 1, 2, 1, 2]
NstdDev = 	[0,2]#[4,4,2,2]# 0, 4, 4, 2, 2]
avWindow = [50,50]#,50]#[50, 50, 50, 50]
probeRotationAngle = [0,45]
saveNames = ['raw','w50_MA_high']#,'w50_MA_high']#['w50_MA_min', 'w50_MA_low','w50_MA_med','w50_MA_high']#['basicMin','basicMed']#'raw','min','low','med','high']
#
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
				if saveFil == True:
					tempSavePath = fileName.split('/')
					filSaveNameEnd = str('/' + saveNames[i] + '/' + saveNames[i] + '_' + tempSavePath[-1].split('.txt')[0] + '.pkl')
					tempSavePath[-1] = 'timeSeries'	
					tempSavePath[tempSavePath.index('rawData')] = 'processedData'			
					filSaveName = str('/'.join(tempSavePath)+filSaveNameEnd)
#					print(rawSaveName)
				if isinstance(tData,pd.DataFrame):
					ttData = Filter(tData,filterType[i],'mean',avWindow[i],'none',filLoops[i],NstdDev[i])
					if probeRotationAngle[j] == 45:
						tttData = transform(ttData)
					else:
						tttData = ttData
#					print(filSaveName)
					tttData.to_pickle(filSaveName)
#	
##	4.	apply averaging and append a final data series
					dataNew  = timeAverage(tttData)
					dataNew['fileName'] = fileName.split("/")[-1]
					if not isinstance(data,pd.DataFrame):
						data = dataNew
					else:
						data=data.append(dataNew)
			print(data)
			tempSavePath = fileName.split('/')
			tempSavePath[-1] = str(tempSavePath[-1].split('.')[0] + '_averaged_' + saveNames[i] + '.pkl')
			tempSavePath[tempSavePath.index('rawData')] = 'processedData'
			globSaveName = '/'.join(tempSavePath)				
			print(globSaveName)
			data.to_pickle(globSaveName)
#
#
##################################################################################################################################
masterPath = [	"../data/processedData/smoothPlate/4Hz/x400/1703XX_4Hz_x400/",
		#	"../data/processedData/smoothPlate/8Hz/x400/171214_8Hz_x400/"	]
		]
filterType = ['high']#['w50_MA_min', 'w50_MA_low','w50_MA_med','w50_MA_high']#[	'basicMin','basicMed']#'min','low','med','high'	]

avDataFileName = [	'4Hz_x400_averaged']#,'8Hz_x400_averaged'	]

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

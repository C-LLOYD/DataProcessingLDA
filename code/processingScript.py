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
from FilterFunctions import Filter
#
##########################################################################################
##
##	1.	Loop over data: Put this in at the end
#
##	2.	import txt file: Give a hard coded name for now
writeData = True
saveFil = True
#
rawPathNames = 	[#"../data/rawData/smoothPlate/4Hz/x400/171211_4Hz_x400/*.txt",
		"../data/rawData/smoothPlate/8Hz/x400/171214_8Hz_x400/*.txt"]

filLoops = 	[0, 1, 2, 1, 2]
NstdDev = 	[0, 4, 4, 2, 2]
saveNames = ['raw','min','low','med','high']
#
if writeData == True:
	for j in range(len(rawPathNames)):
		for i in range(len(filLoops)):
			rawPath = rawPathNames[j]
			data = []
			print(rawPath)
			for fileName in glob.glob(rawPath):
				print(fileName)				
				td = txtToDataFrame(fileName)
#
##	3.	filter data
				if saveFil == True:
					tempSavePath = fileName.split('/')
					filSaveNameEnd = str('/' + saveNames[i] + '/' + saveNames[i] + '_' + tempSavePath[-1].split('.txt')[0] + '.pkl')
					tempSavePath[-1] = 'timeSeries'	
					tempSavePath[tempSavePath.index('rawData')] = 'processedData'			
					filSaveName = str('/'.join(tempSavePath)+filSaveNameEnd)
#					print(rawSaveName)
				if isinstance(td,pd.DataFrame):
					tempData = Filter(td,'movingAverageFilterReynoldsStresses','mean',200,'none',filLoops[i],NstdDev[i])
#					tempData = td
					tempData.to_pickle(filSaveName)
#	
##	4.	apply averaging and append a final data series
					dataNew  = timeAverage(tempData)
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
			data.to_pickle(globSaveName)
#
#
##################################################################################################################################

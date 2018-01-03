###########################################################################################
##############		RAW DATA PROCESSING SCRIPT	###################################
####
####		This script processes raw data from the profile tests of Mar2017.
####
####
###########################################################################################
####
##		Initialise python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from processingFunctions import txtToDataFrame
from processingFunctions import timeAverage
from processingFunctions import plotter
from processingFunctions import findDimensionlessParameters as FDP
from FilterFunctions import Filter
#
##########################################################################################
##
##	1.	loop over all data sets ... use wild cards to identify txt files
##	2.	import txt data and save as temp dataFrame
##	3.	filter this data frame using phase space filter
##	4.	apply averaging and save a new data frame with the columns:
##				PumpSpeed,	Z,	UxMean,	UyMean,	uxRMS,	uyRMS,	uv	(X is constant for each data set - store this in name)
##	5.	remove temp dataframe from storage and read in the next file
##	6.	append the new dataframe with additional rows for each txt file
##	7.	when all data is read, save the data frame and run plotting scripts
##	8.	plot all 5 variables against X, for each pump speed and Z positions (might only be one z position so check)
##
#
##	1.	Loop over data: Put this in at the end
#
##	2.	import txt file: Give a hard coded name for now
writeData = True
saveRaw = True
saveFil = True
testProfiles = False
#
rawPath = 	["../data/rawData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400.0000*"]#.txt"]
#rawPath = 	["../data/rawData/smoothPlate/8Hz/x400/171214_8Hz_x400/*.txt"]
#
#writePath = 	["../Data/processedData/dataFrames/4hz_300mm_profiles.pkl",
#			"../Data/processedData/dataFrames/4hz_400mm_profiles.pkl",
#			"../Data/processedData/dataFrames/8hz_300mm_profiles.pkl",
#			"../Data/processedData/dataFrames/8hz_400mm_profiles.pkl",
#			"../Data/processedData/dataFrames/16hz_300mm_profiles.pkl",
#			"../Data/processedData/dataFrames/16hz_400mm_profiles.pkl"]
#
##
#dataPath = 	["../Data/processedData/dataFrames/4hz_300mm_profiles.pkl",
#			"../Data/processedData/dataFrames/4hz_400mm_profiles.pkl",
#			"../Data/processedData/dataFrames/8hz_300mm_profiles.pkl",
#			"../Data/processedData/dataFrames/8hz_400mm_profiles.pkl",
#			"../Data/processedData/dataFrames/16hz_300mm_profiles.pkl",
#			"../Data/processedData/dataFrames/16hz_400mm_profiles.pkl"]
#
#	Then add new fileName to end?

if writeData == True:
	for j in range(len(rawPath)):
		data = []
#		print(rawPath[j])
		for fileName in glob.glob(rawPath[j]):
			print(fileName)
			if saveRaw == True:
				tempSavePath = fileName.split('/')
				rawSaveNameEnd = str('/raw/' + tempSavePath[-1].split('.txt')[0] + '.pkl')
				tempSavePath[-1] = 'timeSeries'
				tempSavePath[tempSavePath.index('rawData')] = 'processedData'
				rawSaveName = str('/'.join(tempSavePath)+rawSaveNameEnd)				
			td = txtToDataFrame(fileName)
			td.to_pickle(rawSaveName)
#
##	3.	filter data
			if saveFil == True:
				tempSavePath = fileName.split('/')
				filSaveNameEnd = str('/fil/fil.' + tempSavePath[-1].split('.txt')[0] + '.pkl')
				tempSavePath[-1] = 'timeSeries'	
				tempSavePath[tempSavePath.index('rawData')] = 'processedData'			
				filSaveName = str('/'.join(tempSavePath)+filSaveNameEnd)
#				print(rawSaveName)
			if isinstance(td,pd.DataFrame):
				tempData = Filter(td,'movingAverageFilterReynoldsStresses','mean',200,'none')
#				tempData = td
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
		tempSavePath[-1] = str(tempSavePath[-1].split('.')[0]+'_averaged.pkl')
		tempSavePath[tempSavePath.index('rawData')] = 'processedData'
		globSaveName = '/'.join(tempSavePath)				
		data.to_pickle(globSaveName)

if testProfiles == True:
	data = pd.read_pickle('temp.pkl')
	
	plt.semilogx(data.z,data.UxMean,linestyle='-',marker='x')
	plt.show()
	plt.close()
	plt.plot(data.z,data.uxRMS,linestyle='-',marker='x')
	plt.show()
	plt.close()
	plt.plot(data.z,data.uyRMS,linestyle='-',marker='x')
	plt.show()
	plt.close()
	plt.plot(data.z,data.uv,linestyle='-',marker='x')
	plt.show()
	plt.close()
#
#
##################################################################################################################################

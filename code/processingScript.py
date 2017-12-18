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
#
rawPath = 	["../data/rawData/smoothPlate/4Hz/x400/171211_4Hz_x400/*.txt"]
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

if writeData == True:
	for j in range(len(rawPath)):
		data = []
#		print(rawPath[j])
		for fileName in glob.glob(rawPath[j]):
			td = txtToDataFrame(fileName)
#			print(	td.NXYZ[1],
#					td.NXYZ[3], 
#					len(td.sampleNumber),
#					round(np.mean(td.Ux),4),
#					round(np.std(td.Ux)/np.mean(td.Ux),4),
#					round(np.ptp(td.Ux)/np.mean(td.Ux),4)				
#				 )
#			plt.scatter(td.sampleNumber.loc[td.timeStamp < 15],td.Ux.loc[td.timeStamp < 15])
#			plt.show()
#			plt.close()
#			print(len(td.)
#
##	3.	filter data
			if isinstance(tempData,pd.DataFrame):
				tempData = Filter(tempData,'movingAverageFilter','mean',10,'none','none')
#
##	4.	apply averaging and append a final data series
#				dataNew  = timeAverage(tempData)
#				if not isinstance(data,pd.DataFrame):
#					data = dataNew
#				else:
#
#		print(data)
#		data.to_pickle(writePath[j])
#
#
##################################################################################################################################

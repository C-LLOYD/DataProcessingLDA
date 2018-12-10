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
from processingFunctions import timeAverageS
from processingFunctions import timeAverageV
from processingFunctions import transform
from processingFunctions import spatialFilter
from FilterFunctions import Filter
from skinFrictionFunctions import viscousSublayerEstimation as visFunc
from skinFrictionFunctions import logLawEstimation as logFunc
from skinFrictionFunctions import clauserEstimation as clauserFunc
from skinFrictionFunctions import compositeProfile
from skinFrictionFunctions import compositeFit
from skinFrictionFunctions import compositeLeastSquares
from skinFrictionFunctions import minFunc
from skinFrictionFunctions import choiEstimation as choiFunc
from skinFrictionFunctions import wakeFitting
from skinFrictionFunctions import logLawFitting
#
##########################################################################################
##
####	Grouped switches:
####	Raw data from txt files to data frames:
#	1. profile data
writeData = False
saveFil = False
#	2. long time series data (takes a long time to process!!)
writeStatConvergenceSeries = False
#
####	Appending to existing df:
#	1. manual spatial filtering:
writeSpatialFilter = False
#	2. add error bars, calculated from longer TS
writeFlatPlateErrorBars = False
writeDenticlePlateErrorBars = False
#	3. calculate parameters and variables from the profiles
writeSkinFrictionEstimations = False
writeDenticleSkinFriction = False
#
rawPathNames = [
	'../data/rawData/smoothPlate/16Hz/x400/180822_16Hz_x400/*.txt',
	'../data/rawData/riblettedDenticles/16Hz/x400/180823_16Hz_x400/*.txt',
	'../data/rawData/smoothDenticles/16Hz/x400/180824_16Hz_x400/*.txt',
##
	'../data/rawData/smoothPlate/12Hz/x400/180814_12Hz_x400/*.txt',
	'../data/rawData/riblettedDenticles/12Hz/x400/180815_12Hz_x400/*.txt',
	'../data/rawData/smoothDenticles/12Hz/x400/180816_12Hz_x400/*.txt',
#	'../data/rawData/riblettedDenticles/12Hz/x400/180919_12Hz_spatialDep/*.txt',
##
	'../data/rawData/smoothPlate/8Hz/x400/180807_8Hz_x400/*.txt',
	'../data/rawData/riblettedDenticles/8Hz/x400/180808_8Hz_x400/*.txt',
	'../data/rawData/smoothDenticles/8Hz/x400/180809_8Hz_x400/*.txt',
##
	'../data/rawData/smoothPlate/4Hz/x400/180530_4Hz_x400/*.txt',
	'../data/rawData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/*.txt',
	'../data/rawData/smoothDenticles/4Hz/x400/180601_4Hz_x400/*.txt',
#	'../data/rawData/riblettedDenticles/4Hz/x400/180918_4Hz_spatialDep/*.txt',
]
probeRotationAngle = [45.0,45.0]#,45.0,45.0,45.0,45.0,45.0]

filLoops = 	[2,0]#,0,0,0,0,0]#[1,2,1,2]# 0, 1, 2, 1, 2]
NstdDev = 	[3,0]#,0,0,0,0,0]#[4,4,2,2]# 0, 4, 4, 2, 2]
avWindow = [15,0]#,50,50,50,50,50]#,200]#,50]#[50, 50, 50, 50]

saveNames = ['lowFil','raw']

##	1.	Loop over data:
if writeData == True:
	for j in range(len(rawPathNames)):
		for i in range(len(filLoops)):
			rawPath = rawPathNames[j]
			data = []
			print(rawPath)
##	2.	loop through and import all txt files
			for fileName in glob.glob(rawPath):
				print(fileName)
				tData = txtToDataFrame(fileName)
##	3.	filter data and transform data
				if isinstance(tData,pd.DataFrame):
					[ttData,spikeFrac] = Filter(tData,'movingAverageFilter','mean',avWindow[i],'none',filLoops[i],NstdDev[i])
					tttData = transform(ttData,probeRotationAngle[i])
					if saveFil == True:
						tempSavePath = fileName.split('/')
						filSaveNameEnd = str('/' + saveNames[i] + '/' + saveNames[i] + '_' + tempSavePath[-1].split('.txt')[0] + '.pkl')
						tempSavePath[-1] = 'timeSeries'	
						tempSavePath[tempSavePath.index('rawData')] = 'processedData'			
						filSaveName = str('/'.join(tempSavePath)+filSaveNameEnd)
						tttData.to_pickle(filSaveName)
#	
##	4.	apply averaging and append a final data series
					dataNew  = timeAverageS(tttData)
					dataNew['fileName'] = fileName.split("/")[-1]
					dataNew['spikeFrac'] = spikeFrac
					if not isinstance(data,pd.DataFrame):
						data = dataNew
					else:
						data=data.append(dataNew)
			tdataSorted = data.sort_values(by=['z'])
			dataSorted = tdataSorted.reset_index()	
			tempSavePath = fileName.split('/')
			tempSavePath[-1] = str(tempSavePath[-1].split('.')[0] + '_averaged_' + saveNames[i] + '.pkl')
			tempSavePath[tempSavePath.index('rawData')] = 'processedData'
			globSaveName = '/'.join(tempSavePath)				
			print(globSaveName)
			print(dataSorted)
			dataSorted.to_pickle(globSaveName)
#
#
##################################################################################################################################
######
######	Algorithm for writing convergence of stats for longer time series. This is used to estimate errors induced when 
#####		averaging for shorter time windows.
#####
rawPathNames = [
	'../data/rawData/smoothPlate/4Hz/x400/180926_4Hz_x400/*.txt',
	'../data/rawData/smoothPlate/8Hz/x400/180926_8Hz_x400/*.txt',
	'../data/rawData/smoothPlate/12Hz/x400/180926_12Hz_x400/*.txt',
	'../data/rawData/smoothPlate/16Hz/x400/180926_16Hz_x400/*.txt'
]

avTimes = [
	300.0,
	200.0,
	200.0,
	200.0,
]
##	if raw data wanted, append the below lists.
##
probeRotationAngle = [45.0]
filLoops = 	[2]
NstdDev = 	[3]
avWindow = [15]
saveNames = ['lowFil_errors']
#
if writeStatConvergenceSeries == True:
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
					[ttData,spikeFrac] = Filter(tData,'movingAverageFilter','mean',avWindow[i],'none',filLoops[i],NstdDev[i])
					tttData = transform(ttData,probeRotationAngle[i])
#					print(filSaveName)

#	next - do time averaging on time series such that we add columns to ts rather than create a single entry for each point
					ttttData = timeAverageV(tttData,avTimes[j])

					if saveFil == True:
						tempSavePath = fileName.split('/')
						filSaveNameEnd = str('/' + saveNames[i] + '/' + saveNames[i] + '_' + tempSavePath[-1].split('.txt')[0] + '.pkl')
						tempSavePath[-1] = 'timeSeries'	
						tempSavePath[tempSavePath.index('rawData')] = 'processedData'			
						filSaveName = str('/'.join(tempSavePath)+filSaveNameEnd)
						ttttData.to_pickle(filSaveName)
#	Now we must group all time series data into profiles
					print(ttttData['UxMean'].iloc[-1])
					dataNew = pd.DataFrame({
						'fileName': fileName.split("/")[-1],
						'spikeFrac': spikeFrac,
						'x1': pd.Series(float(ttttData.NXYZ[1])),
						'x2': pd.Series(float(ttttData.NXYZ[2])),
						'z': pd.Series(float(ttttData.NXYZ[3])),
						'UxMean': ttttData['UxMean'].iloc[-1],
						'UyMean': ttttData['UyMean'].iloc[-1],
						'uxRMS': ttttData['uxRMS'].iloc[-1],
						'uyRMS': ttttData['uyRMS'].iloc[-1],
						'uv': ttttData['uv'].iloc[-1],
						str('e' + str(int(avTimes[j])) + '_UxMean'): ttttData[str('e' + str(int(avTimes[j])) + '_UxMean')].iloc[0],
						str('e' + str(int(avTimes[j])) + '_UyMean'): ttttData[str('e' + str(int(avTimes[j])) + '_UyMean')].iloc[0],
						str('e' + str(int(avTimes[j])) + '_uxRMS'): ttttData[str('e' + str(int(avTimes[j])) + '_uxRMS')].iloc[0],
						str('e' + str(int(avTimes[j])) + '_uyRMS'): ttttData[str('e' + str(int(avTimes[j])) + '_uyRMS')].iloc[0],
						str('e' + str(int(avTimes[j])) + '_uv'): ttttData[str('e' + str(int(avTimes[j])) + '_uv')].iloc[0],
					})
#
					if not isinstance(data,pd.DataFrame):
						data = dataNew
					else:
						data=data.append(dataNew)
#				print(data)
##	
##	4.	apply averaging and append a final data series
#					dataNew  = timeAverage(tttData)
#					dataNew['fileName'] = fileName.split("/")[-1]
#					if not isinstance(data,pd.DataFrame):
#						data = dataNew
#					else:
#						data=data.append(dataNew)
#			print(data)
			tdataSorted = data.sort_values(by=['z'])
			dataSorted = tdataSorted.reset_index()
			print(dataSorted)
			tempSavePath = fileName.split('/')
			tempSavePath[-1] = str(tempSavePath[-1].split('.')[0] + '_averaged_' + saveNames[i] + '.pkl')
			tempSavePath[tempSavePath.index('rawData')] = 'processedData'
			globSaveName = '/'.join(tempSavePath)				
			print(globSaveName)
			dataSorted.to_pickle(globSaveName)
##################################################################################################################################
####	New method - we manually find spikes/noisey data based on v'v' profiles (most noisey)
if writeSpatialFilter == True:
	data = [
		'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl', 
		'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/riblettedDenticles/8Hz/x400/180808_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothDenticles/8Hz/x400/180809_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/riblettedDenticles/12Hz/x400/180815_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothDenticles/12Hz/x400/180816_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/riblettedDenticles/16Hz/x400/180823_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothDenticles/16Hz/x400/180824_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
	]
	loc1 = [
		[0,11],	#flat plate 4Hz
		[0,4],	#flat plate 8Hz
		[0,11],	#flat plate 12Hz
		[0,11],	#flat plate 16Hz
		[0,19],	#rib plate 4Hz
		[0,19],	#smo plate 4Hz
		[0,17],	#rib plate 8Hz
		[0,7],	#smo plate 8Hz
		[0,10],	#rib plate 12Hz
		[0,9],	#smo plate 12Hz
		[0,11],	#rib plate 16Hz
		[0,6],	#smo plate 16Hz
	]
	loc2 = [
		[],	#flat plate 4Hz
		[14],	#flat plate 8Hz
		[17],	#flat plate 12Hz
		[],	#flat plate 16Hz
		[23],	#rib plate 4Hz
		[],	#smo plate 4Hz
		[],	#rib plate 8Hz
		[15],	#smo plate 8Hz
		[17],	#rib plate 12Hz
		[],	#smo plate 12Hz
		[14],	#rib plate 16Hz
		[],	#smo plate 16Hz
	]
	loc3 = [
		[],	#flat plate 4Hz
		[],	#flat plate 8Hz
		[18],	#flat plate 12Hz
		[],	#flat plate 16Hz
		[],	#rib plate 4Hz
		[],	#smo plate 4Hz
		[],	#rib plate 8Hz
		[18],	#smo plate 8Hz
		[],	#rib plate 12Hz
		[],	#smo plate 12Hz
		[],	#rib plate 16Hz
		[],	#smo plate 16Hz
	]

	for d in range(len(data)):
		import matplotlib.pyplot as plt
		df = pd.read_pickle(data[d])
		print(data[d])
		df['fil'] = True
		df['fil'][loc1[d][0]:loc1[d][1]] = False
		df['fil'][loc2[d]] = False
		df['fil'][loc3[d]] = False
#		print(df)
###	show effect of filtering on v'v' - most robust for finding spikes
#		plt.semilogx(df.z,df.uyRMS)
#		plt.scatter(df.z[~df.fil],df.uyRMS[~df.fil],color='r')
#		plt.scatter(df.z[df.fil],df.uyRMS[df.fil])
#		plt.show()
###	show effect of filtering on all other variables
#		plt.semilogx(df.z[df.fil],df.UxMean[df.fil])
#		plt.scatter(df.z[~df.fil],df.UxMean[~df.fil],color='r')
#		plt.show()
#		plt.semilogx(df.z[df.fil],df.uxRMS[df.fil])
#		plt.scatter(df.z[~df.fil],df.uxRMS[~df.fil],color='r')
#		plt.show()
#		plt.semilogx(df.z[df.fil],df.uyRMS[df.fil])
#		plt.scatter(df.z[~df.fil],df.uyRMS[~df.fil],color='r')
#		plt.show()
#		plt.semilogx(df.z[df.fil],df.uv[df.fil])
#		plt.scatter(df.z[~df.fil],df.uv[~df.fil],color='r')
#		plt.show()
		print(df)
		df.to_pickle(data[d])


#####
##################################################################################################################################
#####
#####		functions to write error bars into profiles
if writeFlatPlateErrorBars == True:
	flatPlateData = [
		'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
	]
	errorData = [
		'../data/processedData/smoothPlate/4Hz/x400/180926_4Hz_x400/180926_sampleTimeExperiments_averaged_lowFil_errors.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/180926_8Hz_x400/180926_sampleTimeExperiments_averaged_lowFil_errors.pkl',
		'../data/processedData/smoothPlate/12Hz/x400/180926_12Hz_x400/180926_sampleTimeExperiments_averaged_lowFil_errors.pkl',
		'../data/processedData/smoothPlate/16Hz/x400/180926_16Hz_x400/180926_sampleTimeExperiments_averaged_lowFil_errors.pkl',
	]
	avTime = [
		'e300',
		'e200',
		'e200',
		'e200',
	]
	var = [
		'UxMean',
		'uxRMS',
		'uyRMS',
		'uv',
	]
# Both flatPlate and denticle plates use different methods to find offset
# Flat plate assumes profiles are identical, but shifted in z

	for d in range(len(flatPlateData)):
		dp = pd.read_pickle(flatPlateData[d])
		de = pd.read_pickle(errorData[d])
	# calculating shift in z:
	# 	prime denote the error matrix, else the profile matrix
	#	We locate minimum velocity of error matrix to define zPrime and Uprime
		Uprime = de.UxMean.min()#.as_matrix()
		zPrime = de.z.loc[de.UxMean == de.UxMean.min()].as_matrix()
		U = dp.UxMean.loc[dp.fil].as_matrix()
	#	Find bounding U and z values of the profile matrix using prime values
		tindex = U < Uprime 
		minU = U[tindex][-1]
		maxU = U[~tindex][0]
		minZ = dp.z.loc[dp.fil].loc[dp.UxMean.loc[dp.fil] == minU].as_matrix()[0]
		maxZ = dp.z.loc[dp.fil].loc[dp.UxMean.loc[dp.fil] == maxU].as_matrix()[0]
	#	Linear interpolation to find deltaZ and then shift error matrix accordingly
		deltaZ = (maxZ-minZ)*(Uprime - minU)/(maxU-minU) + minZ - zPrime
		de["z"] = de["z"] + deltaZ
##	Can check deltaZ through plotting against profile data
#		import matplotlib.pyplot as plt
#		plt.semilogx(dp.z,dp.UxMean)
#		plt.scatter(de.z,de.UxMean)
#		plt.show()
#		print(dp.z[0])
	# Now for each variable we use linear interpolation to expand error matrix to length of profile
		for v in range(len(var)):
			varName = str(avTime[d] + '_' + str(var[v]))
			dp[varName] = 0.0
			for i in range(len(dp.z.loc[dp.fil])):
				tz = dp.z.loc[dp.fil].as_matrix()[i]
				tindex = de.z < tz
	#			cases correspond to start and end of profile where linear interpolation doesn't work. 
	#			We assume error  constant here
				if sum(tindex) == 0:
					dp[varName].loc[dp.z == tz] = de[varName].loc[de.z == de.z.min()].as_matrix()
				elif sum(~tindex) == 0:
					dp[varName].loc[dp.z == tz] = de[varName].loc[de.z == de.z.max()].as_matrix()
				else:
					maxZ = de.z[~tindex].as_matrix()[0]
					minZ = de.z[tindex].as_matrix()[-1]
					maxE = de[varName].loc[de.z == maxZ].as_matrix()
					minE = de[varName].loc[de.z == minZ].as_matrix()
					dp[varName].loc[dp.z == tz] = minE + (maxE-minE)*(tz-minZ)/(maxZ-minZ)
#		import matplotlib.pyplot as plt
#		plt.semilogx(dp.z,dp.e300_UxMean)
#		plt.scatter(de.z,de.e300_UxMean)
#		plt.show()			
		print(dp)
		dp.to_pickle(flatPlateData[d])

if writeDenticlePlateErrorBars == True:
	denticleData = [
		[
			'../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
			'../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		],
##
		[
			'../data/processedData/riblettedDenticles/8Hz/x400/180808_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
			'../data/processedData/smoothDenticles/8Hz/x400/180809_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		],
##
		[
			'../data/processedData/riblettedDenticles/12Hz/x400/180815_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
			'../data/processedData/smoothDenticles/12Hz/x400/180816_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		],
##
		[
			'../data/processedData/riblettedDenticles/16Hz/x400/180823_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
			'../data/processedData/smoothDenticles/16Hz/x400/180824_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
		],
	]
	errorData = [
		'../data/processedData/smoothPlate/4Hz/x400/180926_4Hz_x400/180926_sampleTimeExperiments_averaged_lowFil_errors.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/180926_8Hz_x400/180926_sampleTimeExperiments_averaged_lowFil_errors.pkl',
		'../data/processedData/smoothPlate/12Hz/x400/180926_12Hz_x400/180926_sampleTimeExperiments_averaged_lowFil_errors.pkl',
		'../data/processedData/smoothPlate/16Hz/x400/180926_16Hz_x400/180926_sampleTimeExperiments_averaged_lowFil_errors.pkl',
	]
	avTime = [
		'e300',
		'e200',
		'e200',
		'e200',
	]
	var = [
		'UxMean',
		'uxRMS',
		'uyRMS',
		'uv',
	]
# Flat plate error profiles are ~1mm away from first point in the grid
# Assume this is the same for denticles i.e z_error_min = zmin + 1 (shift whole profile)
	for rate in range(len(denticleData)):
		for plate in range(len(denticleData[0])):
			print(denticleData[rate][plate])
			dp = pd.read_pickle(denticleData[rate][plate])
			de = pd.read_pickle(errorData[rate])

			U = dp.UxMean.loc[dp.fil].as_matrix()/dp.UxMean.loc[dp.fil].max()
			found = False
			i = 0
			while ~found:
				Uprime = de.UxMean.iloc[i]/de.UxMean.max()
				zPrime = de.z.iloc[i]
		#	Find bounding U and z values of the profile matrix using prime values
				tindex = U < Uprime
				if sum(tindex) > 0:	
					found = True
					break
				else:
					i += 1
#			print(tindex)
#			print(zPrime) 
#			print(U)
#			print(Uprime)
			minU = U[tindex][-1]
			maxU = U[~tindex][0]
			minZ = dp.z.loc[dp.fil].loc[dp.UxMean.loc[dp.fil]/dp.UxMean.loc[dp.fil].max() == minU].as_matrix()[0]
			maxZ = dp.z.loc[dp.fil].loc[dp.UxMean.loc[dp.fil]/dp.UxMean.loc[dp.fil].max() == maxU].as_matrix()[0]
		#	Linear interpolation to find deltaZ and then shift error matrix accordingly
			deltaZ = (maxZ-minZ)*(Uprime - minU)/(maxU-minU) + minZ - zPrime
			de["z"] = de["z"] + deltaZ
#			import matplotlib.pyplot as plt
#			plt.semilogx(dp.z,dp.UxMean/dp.UxMean.max())
#			plt.scatter(de.z,de.UxMean/de.UxMean.max())
#			plt.show()
	# Now for each variable we use linear interpolation to expand error matrix to length of profile
			for v in range(len(var)):
				varName = str(avTime[rate] + '_' + str(var[v]))
				dp[varName] = 0.0
				for i in range(len(dp.z.loc[dp.fil])):
					tz = dp.z.loc[dp.fil].as_matrix()[i]
					tindex = de.z < tz
		#			cases correspond to start and end of profile where linear interpolation doesn't work. 
		#			We assume error  constant here
					if sum(tindex) == 0:
						dp[varName].loc[dp.z == tz] = de[varName].loc[de.z == de.z.min()].as_matrix()
					elif sum(~tindex) == 0:
						dp[varName].loc[dp.z == tz] = de[varName].loc[de.z == de.z.max()].as_matrix()
					else:
						maxZ = de.z[~tindex].as_matrix()[0]
						minZ = de.z[tindex].as_matrix()[-1]
						maxE = de[varName].loc[de.z == maxZ].as_matrix()
						minE = de[varName].loc[de.z == minZ].as_matrix()
						dp[varName].loc[dp.z == tz] = minE + (maxE-minE)*(tz-minZ)/(maxZ-minZ)
#			import matplotlib.pyplot as plt
#			plt.semilogx(dp.z,dp.e300_UxMean)
#			plt.scatter(de.z,de.e300_UxMean)
#			plt.show()
			print(dp)
			dp.to_pickle(denticleData[rate][plate])

##################################################################################################################################
#####		minimisation method to curve fit raw data to find the boundary layer thickness
#####
##		
##
testing = True
if testing == True:
#	dfName = '../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl'
#	dfName = '../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl'
	dfName = '../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl'
#	dfName = '../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl'
#	nu = 1.0171e-6
#	nu = 9.842e-7
	nu = 1.006e-6,
#	nu = 9.807e-7
	df = pd.read_pickle(dfName)
	U = df.UxMean.loc[df.fil].as_matrix()
	yTilde = df.z.loc[df.fil].as_matrix()*1e-3
#	lim = 5e-3
	from flatPlateParameterisationFunctions import compFit
	from flatPlateParameterisationFunctions import compositeProfile as comp
	import scipy.optimize as opt
	x0 = [1e-3,1e-5,0.3,5.0]
	x1 = [40e-3,1e-3,0.5,18.0]
#	opt.show_options(solver='minimize', method='Nelder-Mead')#, disp=True)
#	temp = opt.minimize(compFit, x0, args=(U, yTilde, nu, 1, 'none'), method='Nelder-Mead',tol=1e-3,options={'disp':True,'xtol':1e-3,'ftol':1e-3})
#	temp = opt.minimize(compFit, x0, args=(U, yTilde, nu, 1, 'none'), method='Nelder-Mead',tol=1e-3,options={'adaptive':True,'disp':True,'xtol':1e-3,'ftol':1e-3})
#	temp = opt.differential_evolution(compFit, [(3e-3,8e-3),(1e-4,1e-3),(0.3,0.5),(8.0,18.0)],args=(U, yTilde, nu, 1, 'none'))
	import matplotlib.pyplot as plt
	from scipy.optimize._differentialevolution import DifferentialEvolutionSolver as DES
	solver = DES(compFit, [(3e-3,40e-3),(1e-4,1e-3),(0.3,0.5),(8.0,18.0)],args=(U, yTilde, nu, 1, 'none'),mutation=0.5,popsize=200,disp=True)#,callback=True,strategy='rand2bin'),recombination=0.5)
	for i in range(1000):
   		x, e = next(solver)
   		# access the energies
#		print(solver.population_energies)
		uTau = x[0]	
		dy = x[1]
		kappa = x[2]
		a = x[3]
		[delta,wakeCoeffs,UPlus,yPlus] = compFit(x,U,yTilde,nu,1,'none',flag=True)
		Ufit1 = comp(UPlus,yPlus,delta*uTau/nu,kappa,nu,a,1,'none')		
		plt.semilogx(yPlus,UPlus,linestyle='',marker='o')	
		plt.plot(yPlus,Ufit1)
#		plt.ylim([0,20])
		print(i,x,e)
#		plt.plot(i,e,'rx')
#		plt.plot(i,x[0],'kx')
#		plt.plot(i,x[1],'kx')
#		plt.plot(i,x[2],'kx')
#		plt.plot(i,x[3],'kx')
		plt.draw()
		plt.pause(1e-17)
#		time.sleep(0.1)


#	print(temp)
#	[delta,wakeCoeffs,UPlus,yPlus] = compFit(temp.x,U,yTilde,nu,1,'none',flag=True)
#	X0 = [1e-4,1e-4]
#	temp = opt.minimize(genericPolyFit, [X0], args=(yTilde,U,lim), method='Nelder-Mead',tol=1e-16,
#			 			options={'disp':True,'xtol':1e-16,'ftol':1e-16})
#	print(temp)
#	uTau = temp.x[0]
#	dy = temp.x[1]
#	kappa = temp.x[2]
#	a = temp.x[3]

#	from flatPlateParameterisationFunctions import compositeProfile
#	Ufit1 = compositeProfile(UPlus,yPlus,delta*uTau/nu,kappa,nu,a,1,'none')
#	Ufit2 = wakeCoeffs[0]*np.log(yPlus) + wakeCoeffs[1]*yPlus**3 + wakeCoeffs[2]*yPlus**2 + wakeCoeffs[3]*yPlus + wakeCoeffs[4]
#	import matplotlib.pyplot as plt
#	plt.semilogx(yPlus,UPlus,linestyle='',marker='o')
#	plt.plot(yPlus,Ufit1)
#	plt.plot(yPlus,Ufit2)
#	plt.ylim([0,20])
#	plt.xlabel('y+')
#	plt.ylabel('U+')
#	plt.savefig('wakefit2.png')
	
##################################################################################################################################
#####			useful - execfile('processingScript.py')
#####		Skin friction estimates for flat plate
#####
##		This iteration uses a minimisation approach to find the global (and local) minima of a four dimensional solution space
##		This is the most robust method so far - further testing required to properly guage the solution space and ensure we 
##		identify the correct global minima
##
testing = False
if testing == True:
#	dfName = '../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl'
#	dfName = '../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl'
#	dfName = '../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl'
	dfName = '../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl'

#	nu = 1.0171e-6
#	nu = 9.842e-7
#	nu = 1.006e-6,
	nu = 9.807e-7

	df = pd.read_pickle(dfName)
#
##	First step is just looking at linear section of profile
	U = df.UxMean.loc[df.fil].as_matrix()
	yTilde = df.z.loc[df.fil].as_matrix()*1e-3
#
##	Note - can't just look at viscous sublayer for any profiles since there are very few points - too much noise!
#
##	Grid method tests: we use a large grid in 4 dimensions
#	tuTau 	= 	np.linspace(0.0,5e-2,5)
#	tdy    	=	np.linspace(0.0,5e-3,5)
#	tkappa 	=  	np.linspace(0.3,0.5,5)
#	ta 		=	np.linspace(5.0,15.0,5)

#	x0 = [5e-3, 5e-4, 0.4, 10.3]
	
#	x0 = [5e-3,4e-4,0.48,11.3]
#	x0 = [15e-3,4e-4,0.4,11.3]
#	x0 = [20e-3,4e-4,0.25,11.0]

#
##	Initial conditions:
	x0_uTau = np.linspace(18.5e-3,21.5e-3,2)
#	x0_uTau = np.linspace(10e-3,30e-3,5)
	x0_dy = np.linspace(1e-4,1e-3,2)
	x0_kappa = np.linspace(0.3,0.5,2)
	x0_a = np.linspace(9,13,2)

	X0u, X0y, X0k, X0a = np.meshgrid(x0_uTau,x0_dy,x0_kappa,x0_a) 
#
##	Solution matrices: the minimisation function returns optimised X values that minimise the error obtained from the 'testFunc'
	Xu, Xy, Xk, Xa, Xe, index = np.zeros(np.shape(X0u)),np.zeros(np.shape(X0u)),np.zeros(np.shape(X0u)),np.zeros(np.shape(X0u)),np.zeros(np.shape(X0u)),np.zeros(np.shape(X0u))
	import scipy.optimize as opt
	from skinFrictionFunctions import testFunc
	for i in range(len(X0u)):
		for j in range(len(X0y)):
			print(str("completed: " + str(	float(i)/float(len(X0u)) + float(j)/float(len(X0u)*len(X0y))	)	))
			for k in range(len(X0k)):
				for l in range(len(X0a)):
#					temp = opt.minimize(testFunc, [X0u[i,j,k,l],X0y[i,j,k,l],X0k[i,j,k,l],X0a[i,j,k,l]], args=(U, yTilde, nu, 1.0, 'none'), method='L-BFGS-B',tol=1e-8,
#			 			bounds=[[1e-5,50e-3],[0.0,5e-2],[0.2,0.8],[3,18]],options={'eps':1e-8,'disp':False,'gtol':1e-8,'ftol':1e-8})
					temp = opt.minimize(testFunc, [X0u[i,j,k,l],X0y[i,j,k,l],X0k[i,j,k,l],X0a[i,j,k,l]], args=(U, yTilde, nu, 1, 'none'), method='Nelder-Mead',tol=1e-16,
			 			options={'disp':False,'xtol':1e-16,'ftol':1e-16})
					Xu[i,j,k,l] = temp.x[0]
					Xy[i,j,k,l] = temp.x[1]
					Xk[i,j,k,l] = temp.x[2]
					Xa[i,j,k,l] = temp.x[3]
					Xe[i,j,k,l] = temp.fun
					index[i,j,k,l] = temp.success
					if temp.success == False:
						print(temp.success)
#	test2 = testFunc(x0,U, yTilde, nu, 1.0, 'none')
#	print(test2)
#	test = opt.minimize(testFunc, x0, args=(U, yTilde, nu, 1.0, 'none'), method='L-BFGS-B',tol=1e-16,
#			 bounds=[[12e-3,17e-3],[0.0,2e-3],[0.3,0.6],[5,15]],options={'eps':1e-15,'disp':True,'gtol':1e-16,'ftol':1e-16,'maxls':100})
#	test = opt.minimize(testFunc, x0, args=(U, yTilde, nu, 1.0, 'none'), method='L-BFGS-B',tol=1e-16,
#			 bounds=[[1e-4,50e-3],[0.0,5e-3],[0.38,0.38],[10.3,10.3]],options={'eps':1e-15,'disp':False,'gtol':1e-16,'ftol':1e-16})
#	test = opt.basinhopping(testFunc, x0, niter=1, T=0.01, stepsize=0.5, minimizer_kwargs={'args':[U, yTilde, nu, 1.0, 'none'],'method':'L-BFGS-B'})
	#check basinhopping !!


#	print(test)
#	print(test.x)
#	uTau = test.x[0]
#	deltaY = test.x[1]
#	kappa = test.x[2]
#	a = test.x[3]
	locmin = np.unravel_index(Xe.argmin(), Xe.shape)
	uTau = Xu[locmin]
	deltaY = Xy[locmin]
	kappa = Xk[locmin]
	a = Xa[locmin]
	print(uTau,deltaY,kappa,a)
	UPlus = U/uTau
	yPlus = (yTilde-deltaY)*uTau/nu
	Utheory = compositeProfile(UPlus,yPlus,kappa,nu,a,0.04,'channel')
#
	yPlus2 = np.logspace(0,4,10000)
	Utheory2 = compositeProfile(UPlus,yPlus2,kappa,nu,a,1,'none')
	import matplotlib.pyplot as plt
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.semilogx(yPlus,Utheory)
	ax.semilogx(yPlus2,Utheory2)
	ax.scatter(yPlus,UPlus)
	ax.set_xlabel(r"$y^+$",fontsize='30')
	ax.set_ylabel(r"$U^+$",fontsize='30')#,rotation=0,labelpad=45)
	plt.tight_layout()
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
	#writeName = str('../plots/freeStream_'+varNames[var]+'.png')
	#plt.savefig('./16Hz_fit.png')
	plt.show()
	plt.close()

	plt.semilogx(yPlus2,Utheory2 - np.log(yPlus2)/kappa)
	print((Utheory2 - np.log(yPlus2)/kappa))
	plt.show()

	
#	testFunc(x0,uTilde, yTilde, nu, PI, method)



#	yPlus = np.zeros([len(tuTau),len(tdy),len(yTilde)])
#	uPlus = np.zeros([len(tuTau),len(U)])
#	for i in range(len(tuTau)):
#		print('initialising matrices ... ' + str(float(i)/float(len(tuTau))) )
#		uPlus[i][:] = U/tuTau[i]
#		for j in range(len(tdy)):
#			yPlus[i][j][:] = (yTilde - tdy[j])*tuTau[i]/nu




testing = False
if testing == True:
##	Set up variables
#
#	dfName = '../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl'
#	dfName = '../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl'
#	dfName = '../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl'
	dfName = '../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl'

#	nu = 1.0171e-6
#	nu = 9.842e-7
#	nu = 1.006e-6,
	nu = 9.807e-7

	df = pd.read_pickle(dfName)
#
	U = df.UxMean.loc[df.fil].as_matrix()
	yTilde = df.z.loc[df.fil].as_matrix()*1e-3
#
##	Set solution space - for each point in the meshgrid we calculate the rms error between our data and the composite profile
	x_uTau = np.linspace(15e-3,25e-3,80)
	x_dy = np.linspace(2e-4,8e-4,80)
	x_k = np.linspace(0.3,0.5,80)
	x_a = np.linspace(9.3,13.3,80)	
	Xu, Xy, Xk, Xa = np.meshgrid(x_uTau,x_dy,x_k,x_a)
	Xe = np.zeros(np.shape(Xu))
#
##	run testFunc for every point in the grid
	import scipy.optimize as opt
	from skinFrictionFunctions import testFunc
	for i in range(len(x_uTau)):
		for j in range(len(x_dy)):
			print(str("completed: " + str(float(i)/float(len(Xu)) + float(j)/float(len(Xu)*len(Xy)))))
			for k in range(len(x_k)):
				for l in range(len(x_a)):
					Xe[i,j,k,l] = testFunc([Xu[i,j,k,l],Xy[i,j,k,l],Xk[i,j,k,l],Xa[i,j,k,l]],U, yTilde, nu, 1.0, 'none')


	from pyevtk.hl import gridToVTK
	locmin = np.unravel_index(Xe.argmin(), Xe.shape)
	locmax = np.unravel_index(Xe.argmax(), Xe.shape)
	
#	if saveName is not False:
	unorm = (x_uTau-np.min(x_uTau))/(np.max(x_uTau)-np.min(x_uTau))
	ynorm = (x_dy-np.min(x_dy))/(np.max(x_dy)-np.min(x_dy))
	knorm = (x_k-np.min(x_k))/(np.max(x_k)-np.min(x_k))
	anorm = (x_a-np.min(x_a))/(np.max(x_a)-np.min(x_a))
	enorm = (Xe-np.min(Xe))/(np.max(Xe)-np.min(Xe))
#		enorm = (errorGrid - errorGrid[locmin])/(errorGrid[locmax]- errorGrid[locmin])
#		print(enorm)
##		NOTE:- treat
	for i in range(len(unorm)):
		gridToVTK(str('./errorSpace/fine_u_'+str(i)), ynorm,knorm,anorm, pointData={'enorm':enorm[i]})
#	import matplotlib.pyplot as plt
#	plt.contour(enorm[:,:,0,0])
#	plt.show()


##################################################################################################################################
#####					Skin friction estimates for flat plate
##
####			NOTE: - This method is old - it adopts an empirical fit for uTau which is not valid for our boundary layer
####			delete this later.
####
if writeSkinFrictionEstimations == True:
	flatPlateDataNames = [
		'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
	] 
	label = [
		'4Hz',
		'8Hz',
		'12Hz',
		'16Hz',
	]
	nu = [
		1.0171e-6,
		9.842e-7,
		1.006e-6,
		9.807e-7	
		
	]
##
##	Set up grids of kappa, a, and deltaY
	kappaV = np.linspace(0.35,0.45,21)
	deltaYV = np.linspace(0.0,1e-3,101)
	aV = np.linspace(9,12,61)
#	kappaV = [0.384]
#	kappaV = [0.4]
#	aV = [10.3061]
#	aV = [10.3]
	
	PI = 0.4
	
	for d in range(len(flatPlateDataNames)):
		df = pd.read_pickle(flatPlateDataNames[d])
		U = df.UxMean.as_matrix()
		yTilde = df.z.as_matrix()*1e-3	
		[kappa,deltaY,a] = minFunc(U,yTilde,kappaV,deltaYV,aV,nu[d],PI,'none',label[d])


		print(kappa,deltaY,a)
		y = df.z.as_matrix()*1e-3 - deltaY
		U = df.UxMean.as_matrix()
#		y = (yTilde - deltaY)
		[uTau, deltaStar, theta, H] = clauserFunc(U,y,kappa,nu[d])

		print(uTau)
		if 'uTau' in df:
			dfNew = df.drop(['y','nu','uTau','theta','delta', 'deltaStar','H','kappa','a'],axis=1)#.reset_index().copy()
		else:
			dfNew = df.copy()

		dfNew["y"], dfNew["nu"], dfNew["uTau"], dfNew["delta"], dfNew["deltaStar"], dfNew["theta"], dfNew["H"], dfNew["kappa"], dfNew["a"] = \
			np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

		dfNew["y"] = pd.Series(df.z.as_matrix()*1e-3 - deltaY)
		dfNew["nu"] = pd.Series([nu[d]])
		dfNew["uTau"] = pd.Series([uTau])
		dfNew["delta"] = pd.Series([y[U>0.98*np.mean(U[-4:])][0]])
		dfNew["deltaStar"] = pd.Series([deltaStar])
		dfNew["theta"] = pd.Series([theta])
		dfNew["H"] = pd.Series([deltaStar/theta])
		dfNew["kappa"] = pd.Series([kappa])
		dfNew["a"] = pd.Series([a])
#		print(dfNew)
#		dfNew.to_pickle(flatPlateDataNames[d])
##################################################################################################################################
#####
#####		Skin friction estimates for denticle plates
#####
if writeDenticleSkinFriction == True:
	import matplotlib.pyplot as plt
#	dfBase = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl')
#	riblettedDenticleData = pd.read_pickle('../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/4Hz_x400_averaged_lowFil.pkl')
#	print(flatPlateData)
	flatPlateNames = [
#		'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
##
#		'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
###
#		'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
##
#		'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
	]

	denticleNames = [
#		'../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
		'../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
##
#		'../data/processedData/riblettedDenticles/8Hz/x400/180808_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/smoothDenticles/8Hz/x400/180809_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
##
#		'../data/processedData/riblettedDenticles/12Hz/x400/180815_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/smoothDenticles/12Hz/x400/180816_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
##
#		'../data/processedData/riblettedDenticles/16Hz/x400/180823_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
#		'../data/processedData/smoothDenticles/16Hz/x400/180824_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
	]

	errorSaveNames = [
		'test'
	]

	nu = [
#		9.962e-7,
		9.938e-7,
#		9.842e-7,
#		9.938e-7,
#		9.950e-7,
#		9.914e-7,
#		9.842e-7,
#		9.986e-7,
	]
	
	cutoff = [
#		5,
		5,
#		30,
#		30,
#		20,
#		20,
#		30,
#		20,
	]
	uTauVV = [
#		np.linspace(3e-3,6e-3,51),
		np.linspace(4e-3,7e-3,51),
 #		np.linspace(0.9e-2,1.2e-2,51),
# 		np.linspace(0.9e-2,1.2e-2,51),
 #		np.linspace(1.6e-2,2.0e-2,51),
#		np.linspace(1.6e-2,2.0e-2,51),
#	 	np.linspace(1.9e-2,3.0e-2,51),
#	 	np.linspace(1.9e-2,3.0e-2,51),	
	]
	deltaYV = np.linspace(1.0e-3,2.6e-3,51)
	deltaUPlusVV = [
#		np.linspace(-0,0,2),
#		np.linspace(-0,0,2),
#		np.linspace(-0,4,61),
#		np.linspace(-0,4,61),	
#		np.linspace(2,6,61),	
#		np.linspace(2,6,61),	
#		np.linspace(0,0,2),	
#		np.linspace(0,0,2),
	]
	aVV = [
#		np.linspace(9,12,51),
		np.linspace(10,13,51),
#		np.linspace(9,12,51),
#		np.linspace(9,12,51),
#		np.linspace(7,10,51),
#		np.linspace(7,10,51),
#		np.linspace(5,10,51),
#		np.linspace(4,7,31),
	]	
#	uTauV = [5e-3]
#	deltaYV = [2e-3]

	for d in range(len(flatPlateNames)):
		print(denticleNames[d])

		uTauV = uTauVV[d]
#		deltaUPlusV = deltaUPlusVV[d]
		aV = aVV[d]
		dfBase = pd.read_pickle(flatPlateNames[d])
		df = pd.read_pickle(denticleNames[d])
		print(df)

#		deltaY = wakeFitting(dfBase,df,deltaYV)
##		print(dy)
#		[uTau, t, deltaUPlus] = compositeFit(dfBase,df,uTauV,[deltaY,deltaY],deltaUPlusV, nu[d], errorSaveNames[0])
#		print(uTau,deltaUPlus)
#		[uTau, t, deltaUPlus] = logLawFitting(dfBase,df,uTauV,deltaUPlusV,nu[d],[deltaY,deltaY],'test')
#		print(df.z.as_matrix()*1e-3-deltaYV[-1])
#		print(df)
#		print(dfBase)
#		[uTau, deltaY, deltaStar, theta] = choiFunc(dfBase,df,uTauV,deltaYV, nu[d], errorSaveNames[0])
		[uTau, deltaY, a] = compositeFit(dfBase,df[cutoff[d]:],uTauV,deltaYV,aV, nu[d], errorSaveNames[0])
		print(uTau, deltaY, a)
		print(uTau, dfBase.uTau[0])



	
#		plt.semilogx(dfBase.y*dfBase.uTau[0]/dfBase.nu[0], dfBase.UxMean/dfBase.uTau[0],label='flatPlate',marker='o')
#		plt.plot((df.z*1e-3-deltaY)*uTau/nu[d], df.UxMean/uTau,label='denticles',marker='^')
#		plt.plot((df.z*1e-3-deltaY)*uTau/nu[d], df.UxMean/uTau + deltaUPlus,label='denticles',marker='>')
#		plt.show()


		if 'uTau' in df:
			dfNew = df.drop(['y','nu','uTau'],axis=1).reset_index()#.copy()
		else:
			dfNew = df.copy().reset_index()

		dfNew["y"], dfNew["nu"], dfNew["uTau"] = \
			np.nan, np.nan, np.nan
		print(dfNew)
		dfNew["y"] = df.z.as_matrix()*1e-3 - deltaY
#		print(pd.Series(df.z.as_matrix()*1e-3 - deltaY))
		print(dfNew["y"])
		dfNew["nu"] = pd.Series([nu[d]])
		dfNew["uTau"] = pd.Series([uTau])
		print(dfNew)

		dfNew.to_pickle(denticleNames[d])




######	Plotting script for wake portion
#		Xbase = (np.mean(dfBase.UxMean.as_matrix()[-4:])-dfBase.UxMean.as_matrix())/(np.mean(dfBase.UxMean.as_matrix()[-4:])*dfBase.deltaStar[0]/dfBase.delta[0])
#		Ybase = dfBase.y.as_matrix()/dfBase.delta[0]
#		Xbase = dfBase.y.as_matrix()*dfBase.uTau.as_matrix()[0]/(np.mean(dfBase.UxMean.as_matrix()[-4:])*dfBase.deltaStar.as_matrix()[0])	
#		Ybase = (np.mean(dfBase.UxMean.as_matrix()[-4:]) - dfBase.UxMean.as_matrix())/dfBase.uTau.as_matrix()[0]
#		X = (df.z.as_matrix()*1e-3-deltaY)*uTau/(np.mean(df.UxMean.as_matrix()[-4:])*deltaStar)	
#		Y = (np.mean(df.UxMean.as_matrix()[-4:]) - df.UxMean.as_matrix())/uTau
#		plt.plot(Xbase,Ybase,label='base')
#		plt.plot(X,Y)

#	plt.xlim([0,20])
#	plt.ylim([0,0.15])
#		plt.show()















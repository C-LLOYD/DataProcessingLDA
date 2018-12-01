#######
#######	Algorithm for writing spike fraction of filtered series - possibly useful when spatially filtering data
#######
masterPath = [
	'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/',
#	'../data/processedData/riblettedDenticles/16Hz/x400/180823_16Hz_x400/',
#	'../data/processedData/smoothDenticles/16Hz/x400/180824_16Hz_x400/',
##
	'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/',
#	'../data/processedData/riblettedDenticles/12Hz/x400/180815_12Hz_x400/',
#	'../data/processedData/smoothDenticles/12Hz/x400/180816_12Hz_x400/',
#	'../data/processedData/riblettedDenticles/12Hz/x400/180919_12Hz_spatialDep/',
##
	'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/',
#	'../data/processedData/riblettedDenticles/8Hz/x400/180808_8Hz_x400/',
#	'../data/processedData/smoothDenticles/8Hz/x400/180809_8Hz_x400/',
##
	'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/',
#	'../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/',
#	'../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/',
#	'../data/processedData/riblettedDenticles/4Hz/x400/180918_4Hz_spatialDep/',
]
#	"../data/processedData/smoothDenticles/16Hz/x400/180824_16Hz_x400/",
##	"../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/",
#	"../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/",
#	"../data/processedData/riblettedDenticles/16Hz/x400/180426_16Hz_x400/",
#	"../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/",
#	"../data/processedData/riblettedDenticles/8Hz/x400/180322_8Hz_x400/",
#	"../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/",
#	"../data/processedData/riblettedDenticles/4Hz/x400/180427_4Hz_x400/",
#	"../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/"
#]
		
filterType = ['lowFil']#['w50_MA_min', 'w50_MA_low','w50_MA_med','w50_MA_high']#[	'basicMin','basicMed']#'min','low','med','high'	]

avDataFileName = [
#	'16Hz_x400_averaged',
#	'16Hz_x400_averaged',
	'16Hz_x400_averaged',
#	'12Hz_x400_averaged',
#	'12Hz_x400_averaged',
	'12Hz_x400_averaged',
#	'12Hz_spatialDep_averaged',
#	'8Hz_x400_averaged',
#	'8Hz_x400_averaged',
	'8Hz_x400_averaged',
#	'4Hz_x400_averaged',
#	'4Hz_x400_averaged',
	'4Hz_x400_averaged',
#	'4Hz_spatialDep_averaged',
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
#####		Spatial filtering - currently not used
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
if writeSpatialFilterOld == True:
	for i in range(len(dataNames)):
		for j in range(len(nameEnd)):
			data = pd.read_pickle(str(dataNames[i] + nameEnd[j]))
			newData = spatialFilter(data)
			newData.to_pickle(str(dataNames[i] + nameEnd[j]))

#print(newdata)

##################################################################################################################################

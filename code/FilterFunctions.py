#Currently written as a script but will be made into a function at a later date ..
#
#	Script is used to identify and remove spikes in a given data set by using the method
#	of Goring and Nikora (2002)
#
#
##	Initialise python
import numpy as np
import re
import matplotlib.pyplot as mpl
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import erfinv
#
#########################################################################################################################
#
##	Define function
##	This function controls the converging process and outputs the data file at the end
##	It calls the appropriate spike locator function, dependent on whether the user chooses
##	either moving average, global average, or ellipsoid method. If ellipsoid method is used, 
##	need to specify an additional averageMethod in order to choose either MAD or MEAN methods. 
def Filter(data,filterMethod,averageMethod,window,writePaths_figures,Nloops,Nstds):
#
##	Decompose the important components of the dataframe:
	Ux = data.Ux.as_matrix()
	Uy = data.Uy.as_matrix()
	t = data.timeStamp.as_matrix()
	s = data.sampleNumber.as_matrix()
	resT = data.resTime.as_matrix()
	NXYZ = data.NXYZ.as_matrix().astype(np.float)
	NXYZ = NXYZ[~np.isnan(NXYZ)]
#
##
	UyNew = Uy
	UxNew = Ux
	Spikes = np.isnan(UxNew)	#Spikes are initialised based on the number of NANS in UxNew - there should be zero ...
##	Initialise filtered velocity fields
##
##	Define filtering method
	if filterMethod == 'movingAverageFilter':
		print('Moving Average Filter ...')
		from movingAverageFilter import movingAverageFilter as Fil
		fileAppend = 'filtered_moving_average.pkl'
	elif filterMethod == 'phaseSpaceFilter':
		print('Phase Space Filter ...')
		from phaseSpaceFilter import phaseSpaceFilter as Fil
		fileAppend = 'filtered_phase_space.pkl'
	elif filterMethod == 'globalAverageFilter':
		print('Global Average Filter')
		from globalAverageFilter import globalAverageFilter as Fil
		fileAppend = 'filtered_global_average.pkl'
	elif filterMethod == 'movingAverageFilterReynoldsStresses':
		print('Moving Average Filter RSM...')
		from movingAverageFilterReynoldsStresses import movingAverageFilterReynoldsStresses as Fil
		fileAppend = 'filtered_moving_average.pkl'
	else:
		print('No valid filtering method given ...')
#
##	Run the filter twice but no more than this.
	N = 0
	Nmax = Nloops
	while N<Nmax:
		N=N+1
		if filterMethod == 'movingAverageFilterReynoldsStresses':
			Spikes = Fil(UxNew,UyNew,resT,window,Nstds)
		else:
			XSpikes = Fil(UxNew,resT,window,data,averageMethod,writePaths_figures,'Ux',Nstds)
			YSpikes = Fil(UyNew,resT,window,data,averageMethod,writePaths_figures,'Uy',Nstds)
			Spikes = XSpikes + YSpikes
#
		if len(UxNew) == len(UxNew[~Spikes]):#		test==len(UxNew[~Spikes]):
			break
#
##	Spike replacement: 	We have two options, either remove the data or replace
##					with a local average. Currently code simply removes bad data. 
		else:
			UxNew 	= UxNew[~Spikes]
			UyNew 	= UyNew[~Spikes]
			t	= t[~Spikes]
			resT	= resT[~Spikes]
			s	= s[~Spikes]
			print('Number of spikes = '+str(len(Ux)-len(UxNew))+' = '+str((float(len(Ux))-float(len(UxNew)))/(float(len(Ux))/100))+'%')
#			
##	Add variables to existing data frame.
##	Variables first need to be converted to 'pandas.series'
	UxNew = pd.Series(UxNew)
	UyNew = pd.Series(UyNew)
	t = pd.Series(t)
	s = pd.Series(s)
	resT = pd.Series(resT)
	data2 = pd.DataFrame({'NXYZ':pd.Series(NXYZ),'sampleNumber':s,
			'timeStamp':t,'resTime':resT,'Ux':UxNew,'Uy':UyNew})
#
####		NEEDS CHANGING : FLOW RATE IS CURRENTLY HARD CODED INTO THE WRITE PATH!
#	data2.to_pickle(writePath_dataFrames+'x_'+str(int(float(data.NXYZ[1])))+'_z_'+str(int(float(data.NXYZ[3])))+'_data_'+fileAppend)
	return data2;
##
##
##

###########################################################################################
##############		TIME SERIES ANALYSIS SCRIPT		###############################
####
####		This script imports time series dataFrames and plots variables 
####
####
####
###########################################################################################
####
##		Initialise python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from movingAverageFilterReynoldsStresses import movingWeightedAverage as mwa
#
###########################################################################################
fil = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/timeSeries/min/min_4Hz_x400.000049.pkl')
raw = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/timeSeries/raw/raw_4Hz_x400.000049.pkl')
#print(fil,raw)
#filAv = pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/171211_4Hz_x400.pkl')

#	create time series for RS
UxMean = mwa(200,fil.Ux.as_matrix(),fil.resTime.as_matrix())
UyMean = mwa(200,fil.Uy.as_matrix(),fil.resTime.as_matrix())
ruu = (fil.Ux.as_matrix()-UxMean)**2
rvv = (fil.Uy.as_matrix()-UyMean)**2
ruv = (fil.Ux.as_matrix()-UyMean)*(fil.Uy.as_matrix()-UyMean)

plt.subplot(2,3,1)
plt.plot(fil.timeStamp,fil.Ux)
plt.subplot(2,3,2)
plt.plot(fil.timeStamp,fil.Uy)
plt.plot(fil.timeStamp,UyMean)
plt.subplot(2,3,3)
plt.plot(fil.timeStamp,ruu)
plt.subplot(2,3,4)
plt.plot(fil.timeStamp,rvv)
plt.subplot(2,3,5)
plt.plot(fil.timeStamp,ruv)
plt.show()
plt.close()


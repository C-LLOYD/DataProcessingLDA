###########################################################################################
##############		PROFILE ANALYSIS SCRIPT		#####################################
####
####		This script imports profile dataFrames and plots variables 
####
####
####
###########################################################################################
####
##		Initialise python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#
###########################################################################################
testProfiles = True
#
##	import data
fileName = ['../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged.pkl']


if testProfiles == True:
	data = pd.read_pickle(fileName[0])
	print(data)
	plt.subplot(2,2,1)
	plt.semilogx(data.z,data.UxMean,linestyle=' ',marker='x')
	plt.subplot(2,2,2)
	plt.plot(data.z,data.uxRMS,linestyle=' ',marker='x')
	plt.subplot(2,2,3)
	plt.plot(data.z,data.uyRMS,linestyle=' ',marker='x')
	plt.subplot(2,2,4)
	plt.plot(data.z,data.uv,linestyle=' ',marker='x')
	plt.show()
	plt.close()

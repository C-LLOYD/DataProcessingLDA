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
fileName = ['../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_raw.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_min.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_low.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_med.pkl',
		'../data/processedData/smoothPlate/4Hz/x400/171211_4Hz_x400/4Hz_x400_averaged_high.pkl']


if testProfiles == True:
	raw = pd.read_pickle(fileName[0])
	mini = pd.read_pickle(fileName[1])
	low = pd.read_pickle(fileName[2])
	med = pd.read_pickle(fileName[3])
	high = pd.read_pickle(fileName[4])
	print(raw)
	plt.subplot(2,2,1)
	plt.semilogx(raw.z,raw.UxMean,linestyle=' ',marker='x')
	plt.semilogx(mini.z,mini.UxMean,linestyle=' ',marker='o')
	plt.semilogx(low.z,low.UxMean,linestyle=' ',marker='+')
	plt.semilogx(med.z,med.UxMean,linestyle=' ',marker='.')
	plt.semilogx(high.z,high.UxMean,linestyle=' ',marker='*')
	plt.subplot(2,2,2)
	plt.semilogx(raw.z,raw.uxRMS,linestyle=' ',marker='x')
	plt.semilogx(mini.z,mini.uxRMS,linestyle=' ',marker='o')
	plt.semilogx(low.z,low.uxRMS,linestyle=' ',marker='+')
	plt.semilogx(med.z,med.uxRMS,linestyle=' ',marker='.')
	plt.semilogx(high.z,high.uxRMS,linestyle=' ',marker='*')
	plt.subplot(2,2,3)
	plt.semilogx(raw.z,raw.uyRMS,linestyle=' ',marker='x')
	plt.semilogx(mini.z,mini.uyRMS,linestyle=' ',marker='o')
	plt.semilogx(low.z,low.uyRMS,linestyle=' ',marker='+')
	plt.semilogx(med.z,med.uyRMS,linestyle=' ',marker='.')
	plt.semilogx(high.z,high.uyRMS,linestyle=' ',marker='*')
	plt.subplot(2,2,4)
	plt.semilogx(raw.z,raw.uv,linestyle=' ',marker='x')
	plt.semilogx(mini.z,mini.uv,linestyle=' ',marker='o')
	plt.semilogx(low.z,low.uv,linestyle=' ',marker='+')
	plt.semilogx(med.z,med.uv,linestyle=' ',marker='.')
	plt.semilogx(high.z,high.uv,linestyle=' ',marker='*')
	plt.show()
	plt.close()

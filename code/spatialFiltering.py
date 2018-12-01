############################################################
####	script to manually find poor data close to the wall
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
####
data = [
#	'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl', 
#	'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/riblettedDenticles/8Hz/x400/180808_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothDenticles/8Hz/x400/180809_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/riblettedDenticles/12Hz/x400/180815_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothDenticles/12Hz/x400/180816_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/riblettedDenticles/16Hz/x400/180823_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothDenticles/16Hz/x400/180824_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
]
####
for d in range(len(data)):
	df = pd.read_pickle(data[d]).reset_index()
##	loop through length of df.z and slowly remove data from plot
	print(data[d])
	df['fil'] = True
	df['fil'][0:17] = False
#	df['fil'][23] = False
	print(df)
	plt.semilogx(df.z,df.uyRMS)
	plt.scatter(df.z[~df.fil],df.uyRMS[~df.fil],color='r')
	plt.scatter(df.z[df.fil],df.uyRMS[df.fil])
	plt.show()
 #	plt.semilogx(df.z[df.fil],df.UxMean[df.fil])
#	plt.scatter(df.z[~df.fil],df.UxMean[~df.fil],color='r')
#	plt.show()
#	plt.semilogx(df.z[df.fil],df.uxRMS[df.fil])
#	plt.scatter(df.z[~df.fil],df.uxRMS[~df.fil],color='r')
#	plt.show()
#	plt.semilogx(df.z[df.fil],df.uyRMS[df.fil])
#	plt.scatter(df.z[~df.fil],df.uyRMS[~df.fil],color='r')
#	plt.show()
#	plt.semilogx(df.z[df.fil],df.uv[df.fil])
#	plt.scatter(df.z[~df.fil],df.uv[~df.fil],color='r')
#	plt.show()
#	df.to_pickle(data[d])


	

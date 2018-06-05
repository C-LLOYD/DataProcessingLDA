import pandas as pd
import matplotlib.pyplot as plt

data0 = [
	pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil_0_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil_0_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil_0_rotation.pkl'),
]

dataNeg1 = [
	pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil_-1_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil_-1_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil_-1_rotation.pkl'),
]

data1 = [
	pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil_1_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil_1_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil_1_rotation.pkl'),
]

dataNeg5 = [
	pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil_-5_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil_-5_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil_-5_rotation.pkl'),
]

data5 = [
	pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil_5_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil_5_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil_5_rotation.pkl'),
]

data45 = [
	pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil_45_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil_45_rotation.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/16Hz/x400/180426_16Hz_x400/16Hz_x400_averaged_lowFil_45_rotation.pkl'),
]

for i in range(len(data0)):
	plt.semilogx(data0[i].z,data0[i].uxRMS,label='probe1 0')
	plt.semilogx(data0[i].z,data0[i].uyRMS,label='probe2 0')
#	plt.semilogx(dataNeg1[i].z,dataNeg1[i].uxRMS,label='probe1 -1')
#	plt.semilogx(dataNeg1[i].z,dataNeg1[i].uyRMS,label='probe2 -1')
#	plt.semilogx(data1[i].z,data1[i].uxRMS,label='probe1 1')
#	plt.semilogx(data1[i].z,data1[i].uyRMS,label='probe2 1')
#	plt.semilogx(dataNeg5[i].z,dataNeg5[i].uxRMS,label='probe1 -5')
#	plt.semilogx(dataNeg5[i].z,dataNeg5[i].uyRMS,label='probe2 -5')
#	plt.semilogx(data5[i].z,data5[i].uxRMS,label='probe1 5')
#	plt.semilogx(data5[i].z,data5[i].uyRMS,label='probe2 5')
#	plt.semilogx(data45[i].z,data45[i].uxRMS,label='probe1 45')
#	plt.semilogx(data45[i].z,data45[i].uyRMS,label='probe2 45')
	plt.legend()
	plt.show()
#	plt.semilogx(data0[i].z,data0[i].UxMean,label='probe1 0')
#	plt.semilogx(data0[i].z,-1*data0[i].UyMean,label='probe2 0')
#	plt.semilogx(dataNeg1[i].z,dataNeg1[i].UxMean,label='probe1 -1')
#	plt.semilogx(dataNeg1[i].z,dataNeg1[i].UyMean,label='probe2 -1')
#	plt.semilogx(data1[i].z,data1[i].UxMean,label='probe1 1')
#	plt.semilogx(data1[i].z,data1[i].UyMean,label='probe2 1')
#	plt.semilogx(dataNeg5[i].z,dataNeg5[i].UxMean,label='probe1 -5')
#	plt.semilogx(dataNeg5[i].z,dataNeg5[i].UyMean,label='probe2 -5')
#	plt.semilogx(data5[i].z,data5[i].UxMean,label='probe1 5')
#	plt.semilogx(data5[i].z,data5[i].UyMean,label='probe2 5')
#	plt.legend()
#	plt.show()

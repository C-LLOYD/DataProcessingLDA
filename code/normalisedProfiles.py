####################################################################
####	Plots normalised variables using shear velocity
####
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from skinFrictionFunctions import compositeProfile as comp
####
flatPlateNames = [
	'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothPlate/4Hz/x400/180530_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
##
	'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothPlate/8Hz/x400/180807_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
##
	'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothPlate/12Hz/x400/180814_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
##
	'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothPlate/16Hz/x400/180822_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
]

ribNames = [
	'../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
##
	'../data/processedData/riblettedDenticles/8Hz/x400/180808_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothDenticles/8Hz/x400/180809_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
##
	'../data/processedData/riblettedDenticles/12Hz/x400/180815_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothDenticles/12Hz/x400/180816_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
##
	'../data/processedData/riblettedDenticles/16Hz/x400/180823_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
#	'../data/processedData/smoothDenticles/16Hz/x400/180824_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
]

smoNames = [
#	'../data/processedData/riblettedDenticles/4Hz/x400/180531_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
	'../data/processedData/smoothDenticles/4Hz/x400/180601_4Hz_x400/4Hz_x400_averaged_lowFil.pkl',
##
#	'../data/processedData/riblettedDenticles/8Hz/x400/180808_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
	'../data/processedData/smoothDenticles/8Hz/x400/180809_8Hz_x400/8Hz_x400_averaged_lowFil.pkl',
##
#	'../data/processedData/riblettedDenticles/12Hz/x400/180815_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
	'../data/processedData/smoothDenticles/12Hz/x400/180816_12Hz_x400/12Hz_x400_averaged_lowFil.pkl',
##
#	'../data/processedData/riblettedDenticles/16Hz/x400/180823_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
	'../data/processedData/smoothDenticles/16Hz/x400/180824_16Hz_x400/16Hz_x400_averaged_lowFil.pkl',
]

saveName = [
	'../data/processedData/figures/fittingDependence/4Hz_',
	'../data/processedData/figures/fittingDependence/8Hz_',
	'../data/processedData/figures/fittingDependence/12Hz_',
	'../data/processedData/figures/fittingDependence/16Hz_',
]

nuDenticles = [
	9.962e-7,
#	9.938e-7,
	9.842e-7,
#	9.938e-7,
	9.950e-7,
#	9.914e-7,
	9.842e-7,
#	9.986e-7,
	]

############################################################################################################
####	1. plot normalised profiles for flat plates - comparing data to fitting comp. profile
## 

for d in range(len(flatPlateNames)):
	df = pd.read_pickle(flatPlateNames[d])
	print(df.kappa[0], df.uTau[0], df.a[0], (df.z*1e-3-df.y)[0], df.deltaStar[0]*np.mean(df.UxMean[-4:])/df.nu[0])
	uTau = df.uTau[0]
	nu = df.nu[0]
	yPlus = df.y.as_matrix()*uTau/nu
	UPlus = df.UxMean.as_matrix()/uTau
	Ucomp = comp('temp',yPlus,df.kappa[0],nu,df.a[0],1,'none')
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
#np.logspace(-1,3,100)
	ax.semilogx(yPlus,UPlus,'ko')
	ax.semilogx(yPlus,Ucomp,'r',linewidth=2)
	ax.set_xlabel(r'$y^+$',fontsize='30')
	ax.set_ylabel(r'$U^+$',fontsize='30')
	plt.tight_layout()
#	box = ax.get_position()
#	ax.legend(loc='upper left',numpoints=1,prop={'size':18})
	ax.set_ylim([0, 22])
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
#	plt.savefig(saveName[d] + 'UxMean.png')
#	plt.show()
	plt.close()

################################################################################################
#### 2. compare sharkskin to flat plate for first and second order stats
##
## 2a) means

for d in range(len(flatPlateNames)):
	dfFlat = pd.read_pickle(flatPlateNames[d])
	dfRib = pd.read_pickle(ribNames[d])
	dfSmo = pd.read_pickle(smoNames[d])
#	print(dfSmo,dfRib)
##
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
#np.logspace(-1,3,100)
	ax.semilogx(dfFlat.y*dfFlat.uTau[0]/dfFlat.nu[0],dfFlat.UxMean/dfFlat.uTau[0],'ko',label='Flat Plate')
	ax.semilogx(dfSmo.y*dfSmo.uTau[0]/dfSmo.nu[0],dfSmo.UxMean/dfSmo.uTau[0],'r^',label='Smo. Denticles')
	ax.semilogx(dfRib.y*dfRib.uTau[0]/dfRib.nu[0],dfRib.UxMean/dfRib.uTau[0],'g*',label='Rib. Denticles')
	ax.set_xlabel(r'$y^+$',fontsize='30')
	ax.set_ylabel(r'$U^+$',fontsize='30')
	plt.tight_layout()
	box = ax.get_position()
	ax.legend(loc='upper left',numpoints=1,prop={'size':18})
	ax.set_ylim([0, 22])
	ax.set_xlim([1, 5000])
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
	plt.savefig(saveName[d] + 'sharkskin_UxMean.png')
#	plt.show()
	plt.close()

## 2b) rms x

for d in range(len(flatPlateNames)):
	dfFlat = pd.read_pickle(flatPlateNames[d])
	dfRib = pd.read_pickle(ribNames[d])
	dfSmo = pd.read_pickle(smoNames[d])
#	print(dfSmo,dfRib)
##
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
#np.logspace(-1,3,100)
#	ax.semilogx(dfFlat.y*dfFlat.uTau[0]/dfFlat.nu[0],(dfFlat.uxRMS)**2/((dfFlat.uTau[0])**2.0),'ko',label='Flat Plate')
#	ax.semilogx(dfSmo.y*dfSmo.uTau[0]/dfSmo.nu[0],(dfSmo.uxRMS)**2/((dfFlat.uTau[0])**2.0),'r^',label='Smo. Denticles')
#	ax.semilogx(dfRib.y*dfRib.uTau[0]/dfRib.nu[0],(dfRib.uxRMS)**2/((dfFlat.uTau[0])**2.0),'g*',label='Rib. Denticles')
	ax.semilogx(dfFlat.y*dfFlat.uTau[0]/dfFlat.nu[0],dfFlat.uxRMS/np.mean(dfFlat.UxMean[-4:]),'ko',label='Flat Plate')
	ax.semilogx(dfSmo.y*dfSmo.uTau[0]/dfSmo.nu[0],dfSmo.uxRMS/np.mean(dfSmo.UxMean[-4:]),'r^',label='Smo. Denticles')
	ax.semilogx(dfRib.y*dfRib.uTau[0]/dfRib.nu[0],dfRib.uxRMS/np.mean(dfRib.UxMean[-4:]),'g*',label='Rib. Denticles')
	ax.set_xlabel(r'$y^+$',fontsize='30')
	ax.set_ylabel(r"$\sqrt{\overline{u'u'}}/U_\infty$",fontsize='30')
	plt.tight_layout()
	box = ax.get_position()
	ax.legend(loc='lower right',numpoints=1,prop={'size':18})
	ax.set_ylim([0, 0.16])
	ax.set_xlim([1, 5000])
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
	plt.savefig(saveName[d] + 'sharkskin_uxRMS.png')
#	plt.show()
	plt.close()
#
## 2b) rms y

for d in range(len(flatPlateNames)):
	dfFlat = pd.read_pickle(flatPlateNames[d])
	dfRib = pd.read_pickle(ribNames[d])
	dfSmo = pd.read_pickle(smoNames[d])
#	print(dfSmo,dfRib)
##
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
#np.logspace(-1,3,100)
#	ax.semilogx(dfFlat.y*dfFlat.uTau[0]/dfFlat.nu[0],(dfFlat.uyRMS)**2/((dfFlat.uTau[0])**2.0),'ko',label='Flat Plate')
#	ax.semilogx(dfSmo.y*dfSmo.uTau[0]/dfSmo.nu[0],(dfSmo.uyRMS)**2/((dfFlat.uTau[0])**2.0),'r^',label='Smo. Denticles')
#	ax.semilogx(dfRib.y*dfRib.uTau[0]/dfRib.nu[0],(dfRib.uyRMS)**2/((dfFlat.uTau[0])**2.0),'g*',label='Rib. Denticles')
	ax.semilogx(dfFlat.y*dfFlat.uTau[0]/dfFlat.nu[0],(dfFlat.uyRMS)/np.mean(dfFlat.UxMean[-4:]),'ko',label='Flat Plate')
	ax.semilogx(dfSmo.y*dfSmo.uTau[0]/dfSmo.nu[0],(dfSmo.uyRMS)/np.mean(dfSmo.UxMean[-4:]),'r^',label='Smo. Denticles')
	ax.semilogx(dfRib.y*dfRib.uTau[0]/dfRib.nu[0],(dfRib.uyRMS)/np.mean(dfRib.UxMean[-4:]),'g*',label='Rib. Denticles')
	ax.set_xlabel(r'$y^+$',fontsize='30')
	ax.set_ylabel(r"$\sqrt{\overline{v'v'}}/U_\infty$",fontsize='30')
	plt.tight_layout()
	box = ax.get_position()
	ax.legend(loc='lower right',numpoints=1,prop={'size':18})
	ax.set_ylim([0, 0.065])
	ax.set_xlim([1, 5000])
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
	plt.savefig(saveName[d] + 'sharkskin_uyRMS.png')
#	plt.show()
	plt.close()

## 2b) uv

for d in range(len(flatPlateNames)):
	dfFlat = pd.read_pickle(flatPlateNames[d])
	dfRib = pd.read_pickle(ribNames[d])
	dfSmo = pd.read_pickle(smoNames[d])
#	print(dfSmo,dfRib)
##
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
#np.logspace(-1,3,100)
#	ax.semilogx(dfFlat.y*dfFlat.uTau[0]/dfFlat.nu[0],(dfFlat.uv)/((dfFlat.uTau[0])**2.0),'ko',label='Flat Plate')
#	ax.semilogx(dfSmo.y*dfSmo.uTau[0]/dfSmo.nu[0],(dfSmo.uv)/((dfFlat.uTau[0])**2.0),'r^',label='Smo. Denticles')
#	ax.semilogx(dfRib.y*dfRib.uTau[0]/dfRib.nu[0],(dfRib.uv)/((dfFlat.uTau[0])**2.0),'g*',label='Rib. Denticles')
	ax.semilogx(dfFlat.y*dfFlat.uTau[0]/dfFlat.nu[0],dfFlat.uv/(np.mean(dfFlat.UxMean[-4:])**2),'ko',label='Flat Plate')
	ax.semilogx(dfSmo.y*dfSmo.uTau[0]/dfSmo.nu[0],dfSmo.uv/(np.mean(dfSmo.UxMean[-4:])**2),'r^',label='Smo. Denticles')
	ax.semilogx(dfRib.y*dfRib.uTau[0]/dfRib.nu[0],dfRib.uv/(np.mean(dfRib.UxMean[-4:])**2),'g*',label='Rib. Denticles')
	ax.set_xlabel(r'$y^+$',fontsize='30')
	ax.set_ylabel(r"$\overline{u'v'}/U_\infty^2$",fontsize='30')
	plt.tight_layout()
	box = ax.get_position()
	ax.legend(loc='lower right',numpoints=1,prop={'size':18})
	ax.set_ylim([-0.004, 0])
	ax.set_xlim([1, 5000])
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
	plt.savefig(saveName[d] + 'sharkskin_uv.png')
#	plt.show()
	plt.close()


################################################################################################
####	3.	Cf plot 
ut0 = np.zeros(len(flatPlateNames))
utR = np.zeros(len(flatPlateNames))
utS = np.zeros(len(flatPlateNames))
for d in range(len(flatPlateNames)):
	ut0[d] = pd.read_pickle(flatPlateNames[d]).uTau[0]
	utR[d] = pd.read_pickle(ribNames[d]).uTau[0]
	utS[d] = pd.read_pickle(smoNames[d]).uTau[0]

DRR = utR**2/ut0**2
DRS = utS**2/ut0**2

wPlus = 4e-3*ut0/1e-6
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
plt.rc('text', usetex=True)
plt.plot(wPlus,DRR,'o',label='Ribletted')
plt.plot(wPlus,DRS,'*',label='Smooth')
ax.set_xlabel(r'$W^+$',fontsize='30')
ax.set_ylabel(r"$C_f/C_{f0}$",fontsize='30')
plt.tight_layout()
box = ax.get_position()
ax.legend(loc='upper left',numpoints=1,prop={'size':18})
plt.minorticks_on()
plt.grid(True, which='minor',alpha=0.6)
plt.grid(True, which='major',linewidth=0.9)
plt.savefig('../data/processedData/figures/fittingDependence/dragReduction.png')
#plt.show()
plt.close()




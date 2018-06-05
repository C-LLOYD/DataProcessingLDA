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
###########			Plotter function
def plotter(X,Y,dataLabels,xlims,ylims,markers,lineStyle,xLabel,yLabel,writeName):
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':'20'})
	plt.rc('text', usetex=True)
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	for i in range(len(X)):
		ax.semilogx(X[i],Y[i],label = dataLabels[i], marker=markers[i], linestyle=lineStyle[i])#, color = 'k')
	ax.set_xlabel(xLabel,fontsize='30')
	ax.set_ylabel(yLabel,fontsize='30')#,rotation=0,labelpad=45)
	plt.tight_layout()
#	h.set_rotation(0)
	if legend == True:
		box = ax.get_position()
#		ax.set_position([box.x0, box.y0*1.2, box.width * 0.75, box.height])
		ax.legend(loc='upper left',numpoints=1,prop={'size':18})#, bbox_to_anchor=(1, 0.5),numpoints=1)	
#		plt.legend(loc = 'best')#handles=[plot,],loc=4,prop={'size':18},ncol=1)
	ax.set_ylim(ylims)
	ax.set_xlim(xlims)
	plt.minorticks_on()
	plt.grid(True, which='minor',alpha=0.6)
	plt.grid(True, which='major',linewidth=0.9)
	#writeName = str('../plots/freeStream_'+varNames[var]+'.png')
	plt.savefig(writeName)
	#plt.show()
	plt.close()
#
#########################################################################################
data = [
#	pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_raw.pkl'),
	pd.read_pickle('../data/processedData/smoothPlate/4Hz/x400/180427_4Hz_x400/4Hz_x400_averaged_lowFil.pkl'),
#	pd.read_pickle('../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_raw.pkl'),
#	pd.read_pickle('../data/processedData/smoothPlate/8Hz/x400/180322_8Hz_x400/8Hz_x400_averaged_lowFil.pkl')
]
print(data[0])
nu = 1e-6
kappa = 0.41#[0.36,0.38,0.4,0.42]
B = 6#[5,6,7,8]#np.linspace(5,6,7)
legend = True
for d in range(len(data)):
#	for b in range(len(B)):
#		yPlus = data[d]["y"][data[d]["fil"]]*data[d]["uTau"][0]*1e-3/nu
	yPlus = data[d]["y"]*data[d]["uTau"][0]*1e-3/nu
	print(yPlus)
#		yPlus = data[d]["y"][data[d]["fil"]]*np.mean(data[d]["uTau"])*1e-3/nu
#		uPlus = data[d]["UxMean"][data[d]["fil"]]/data[d]["uTau"][0]
	uPlus = data[d]["UxMean"]/data[d]["uTau"][0]
	uuPlus = data[d]["uxRMS"]**2/data[d]["uTau"][0]**2
	vvPlus = data[d]["uyRMS"]**2/data[d]["uTau"][0]**2
	uvPlus = data[d]["uv"]/data[d]["uTau"][0]**2
#		uPlus = data[d]["UxMean"][data[d]["fil"]]/np.mean(data[d]["uTau"])
#	print(B[b])
	plotter(	[yPlus,yPlus[yPlus<20],yPlus],
		[uPlus,yPlus[yPlus<20],np.log(yPlus)/(kappa) + B],
		[r"$U^+ = U/u_\tau$",r"$U^+ = y^+$",r"$U^+ = \frac{1}{\kappa} \text{ln} y^+ + B$"],[1,500],[0,22],
		['x','','',],['','-','-'],r'$y^+$',r'$U^+$','../data/processedData/figures/1803XX_profiles/lawOfTheWallU2.png')
	legend = True
	plotter(	[yPlus,yPlus,yPlus],
		[uuPlus,vvPlus,uvPlus],
		[r"$\overline{u'^2}$",r"$\overline{v'^2}$",r"$\overline{u'v'}$"],[1,500],[-1.5,10],['x','o','^',],['','','']
		,r'$y^+$',r"$\overline{u_i'u_i'}/u_\tau^2$",'../data/processedData/figures/1803XX_profiles/lawOfTheWalluu2.png')
	








###############################################################################################

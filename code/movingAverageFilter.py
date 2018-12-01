#Currently written as a script but will be made into a function at a later date ..
#
#	Script is used to identify and remove spikes in a given data set by using the method
#	of Goring and Nikora (2002)
#
#
##	Initialise python
import numpy as np
import matplotlib.pyplot as plt
#
########################################################################################################################
##
##	MOVING AVERAGE FILTER
##
##	Define the filtering function:
##	Input: velocity and the averaging window (multiple of 2)
##	Output: index of spikes after several loops.
def movingAverageFilter(U,resT,window,data,method,writePaths_figures,VariableName,Nstds):
#	Half the window for consistency
	W = int(window/2)
#
##	We loop through the length of U and take neighbouring points ....
##	For the first and last few points we can't do this! - Just use a slight weighting.
##	Define averaging window:
	t = resT
	Umeans = np.zeros(len(U))
	Umeans[:] = np.NAN
	std = np.zeros(len(U))
	std[:] = np.NAN
#
##	Calculate means
	Umeans[0:W]	= np.divide(sum(U[0:2*W]*t[0:2*W]),sum(t[0:2*W]))
	Umeans[-W:] = np.divide(sum(U[-2*W:]*t[-2*W:]),sum(t[-2*W:]))		
	for i in range(W,len(U)-W+1):
		Umeans[i] = np.divide(sum(U[i-W:i+W]*t[i-W:i+W]),sum(t[i-W:i+W]))
#
##	Subtract means from series
	std[0:W]	= np.sqrt(np.divide(sum(((U[0:2*W] - Umeans[0:2*W])**2)*t[0:2*W]),sum(t[0:2*W])))
	std[-W:]	= np.sqrt(np.divide(sum(((U[-2*W:] - Umeans[-2*W:])**2)*t[-2*W:]),sum(t[-2*W:])))
#
	for i in range(W,len(U)-W+1):
		std[i]  = np.sqrt(np.divide(sum(((U[i-W:i+W] - Umeans[i-W:i+W])**2)*t[i-W:i+W]),sum(t[i-W:i+W])))
#
##	


#
##	Test : set up range such that if Umeans-2*std < U < Umeans+2*std we use Unew = U, Unew[test] = np.nan
	spikes = (U<Umeans-Nstds*std)+(U>Umeans+Nstds*std)
#
##	for visualisation of filter uncomment below:
#	temp = U.copy()
#	temp[~spikes] = np.nan
#	plt.plot(U,linestyle='',marker='.')
#	plt.plot(temp,linestyle='',marker='o')
#	plt.plot(Umeans)
#	plt.plot(Umeans-Nstds*std)
#	plt.plot(Umeans+Nstds*std)
#	plt.show()
#	plt.close()
	return spikes

########################################################################################################################


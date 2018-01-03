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
##
##	Define calculation of moving average mean
def movingWeightedAverage(win,phi,tau):
#	win is half averaging window
#	phi is variable
#	tau is weighting
	Phi = np.zeros(len(phi))
	Phi[:] = np.NAN
	Phi[0:win] = np.divide(sum(phi[0:2*win]*tau[0:2*win]),sum(tau[0:2*win]))
	Phi[-win:] = np.divide(sum(phi[-2*win:]*tau[-2*win:]),sum(tau[-2*win:]))		
	for i in range(win,len(phi)-win+1):
		Phi[i] = np.divide(sum(phi[i-win:i+win]*tau[i-win:i+win]),sum(tau[i-win:i+win]))	
	return Phi	

def movingAverageFilterReynoldsStresses(u,v,resT,window):
#	Half the window for consistency
	W = int(window/2)
	N = np.linspace(1,len(u),len(u))
	U = movingWeightedAverage(W,u,resT)
	V = movingWeightedAverage(W,v,resT)
	ruu = (u-U)**2
	rvv = (v-V)**2
	ruv = (u-U)*(v-V)
	Ruu = movingWeightedAverage(W,ruu,resT)
	Rvv = movingWeightedAverage(W,rvv,resT)
	Ruv = movingWeightedAverage(W,ruv,resT)
	varRuu = movingWeightedAverage(W,(ruu-Ruu)**2,resT)
	varRvv = movingWeightedAverage(W,(rvv-Rvv)**2,resT)
	varRuv = movingWeightedAverage(W,(ruv-Ruv)**2,resT)
	tol = 4
	spikes = (	  (u < U - tol*np.sqrt(Ruu)) + (u > U + tol*np.sqrt(Ruu)) +
			  (v < V - tol*np.sqrt(Rvv)) + (v > V + tol*np.sqrt(Rvv)) +
			  (ruu < Ruu - tol*np.sqrt(varRuu)) + (ruu > Ruu + tol*np.sqrt(varRuu)) +
			  (rvv < Rvv - tol*np.sqrt(varRvv)) + (rvv > Rvv + tol*np.sqrt(varRvv)) +
			  (ruv < Ruv - tol*np.sqrt(varRuv)) + (ruv > Ruv + tol*np.sqrt(varRuv))		)
#
	plot = False
	if plot == True:
		plt.subplot(2,3,1)
		plt.plot(N[~spikes],u[~spikes],color='k')
		plt.plot(N[spikes],u[spikes],color='r',linestyle=' ',marker='.')
		plt.plot(N,U)
		plt.plot(N,U - tol*np.sqrt(Ruu))
		plt.plot(N,U + tol*np.sqrt(Ruu))
		plt.xlabel('N')
		plt.ylabel('u')
		plt.subplot(2,3,4)
		plt.plot(N[~spikes],v[~spikes],color='k')
		plt.plot(N[spikes],v[spikes],color='r',linestyle=' ',marker='.')
		plt.plot(N,V)
		plt.plot(N,V - tol*np.sqrt(Rvv))
		plt.plot(N,V + tol*np.sqrt(Rvv))
		plt.xlabel('N')
		plt.ylabel('v')
		plt.subplot(2,3,2)
		plt.plot(N[~spikes],ruu[~spikes],color='k')
		plt.plot(N[spikes],ruu[spikes],color='r',linestyle=' ',marker='.')
		plt.plot(N,Ruu)
		plt.plot(N,Ruu - tol*np.sqrt(varRuu))
		plt.plot(N,Ruu + tol*np.sqrt(varRuu))
		plt.xlabel('N')
		plt.ylabel('ruu')
		plt.subplot(2,3,5)
		plt.plot(N[~spikes],rvv[~spikes],color='k')
		plt.plot(N[spikes],rvv[spikes],color='r',linestyle=' ',marker='.')
		plt.plot(N,Rvv)
		plt.plot(N,Rvv - tol*np.sqrt(varRvv))
		plt.plot(N,Rvv + tol*np.sqrt(varRvv))
		plt.xlabel('N')
		plt.ylabel('rvv')
		plt.subplot(2,3,3)
		plt.plot(N[~spikes],ruv[~spikes],color='k')
		plt.plot(N[spikes],ruv[spikes],color='r',linestyle=' ',marker='.')
		plt.plot(N,Ruv)
		plt.plot(N,Ruv - tol*np.sqrt(varRuv))
		plt.plot(N,Ruv + tol*np.sqrt(varRuv))
		plt.xlabel('N')
		plt.ylabel('ruv')
		plt.show()
		plt.close()
	return spikes

########################################################################################################################


###############################################################################
########
########		Skin friction functions
##
##		These functions estimate uTau and y-offset from a dataset using 
##		several methods.
##
##		1.	estimated from viscous sub-layer
##
###############################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#
##	Basic linear regression function: returns alpha, beta, and rSqrd, such that
##	X = beta*Y+alpha
def linearReg(X,Y):
	beta = np.sum( (X-np.mean(X)) * (Y-np.mean(Y)) ) / np.sum( (Y - np.mean(Y))**2 ) 
	alpha = np.mean(X) - beta*np.mean(Y)
	rSqrd = np.sum((X - alpha - beta*Y)**2)
	return [alpha,beta,rSqrd]

###
##	use linear regression to evaluate friction velocity from near wall
##	velocity measurements
##
##	U and y are numpy vectors, nu is float
def viscousSublayerEstimation(U,y,nu):
	yVis = y[y<5.0]
	UVis = U[y<5.0]
	[a,b,e] = linearReg(UVis,yVis)
#	print(len(yVis))
	while (e>1e-5):
		if len(yVis)<5:
			break
		yVis = yVis[:-1]
		UVis = UVis[:-1] 
#		print(len(yVis))
		[a,b,e] = linearReg(UVis,yVis)
#		print(np.sqrt(b*nu))
#		plt.plot(yVis,yVis*b + a)
#		plt.plot(yVis,UVis,linestyle='',marker='x')
#		plt.show()
#		print(eNew)
	return [a,b,e]
###############################################################################
##	Assumes offset has been taken from y vector (i.e yTrue)
def logLawEstimation(U,y,uTauEst,kappa,nu):
	yPlusEst = y*uTauEst/nu
#	print(yPlusEst)
	yLog = y[(yPlusEst > 40) & (yPlusEst < 300)]
	ULog = U[(yPlusEst > 40) & (yPlusEst < 300)]
	[a,b,e] = linearReg(ULog,np.log(yLog))
#	while (e>1e-5):
#		if len(yLog)<10:
#			break
#		yLog = yLog[:-1]
#		ULog = ULog[:-1]
#		[a,b,e] = linearReg(ULog,np.log(yLog))
#	print(uTauEst, b*kappa)
#		print(e)
#		plt.semilogx(y,np.log(y)*b + a)
#		plt.semilogx(y,U,linestyle='',marker='x')
#		plt.show()
	return [a,b,e]


################################################################################
##	Clauser method
def clauserEstimation(U,y,kappa,nu):
	U0 = np.mean(U[-4:])
	deltaStar = np.trapz(1-U/U0,y)
	theta = np.trapz(U*(1-U/U0)/U0,y)
#	print(deltaStar,theta,deltaStar/theta)
	ReDeltaStar = U0*deltaStar/nu
	uTau = U0/(np.log(ReDeltaStar)/kappa + 4.9)
#	print(uTau)
	return [uTau,deltaStar,theta, deltaStar/theta]

################################################################################
###		Choi rough surface uTau estimate and deltaY estimate
def choiEstimation(US,yS,uTauS,deltaStarS,UR,yTildeR):
	U0S = np.mean(US[-4:])
	plt.plot((U0S - US)/uTauS,yS*uTauS/(U0S*deltaStarS), marker='x')
	plt.show()
##
##
##		NOTE: Does it matter if we can't predict uTau properly for flat surface?
##			Do we still get a correct relative error? 








################################################################################								

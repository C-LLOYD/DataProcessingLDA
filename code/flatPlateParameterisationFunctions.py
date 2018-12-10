###############################################################################
########
########		Parameterisation for flat plate cases
##
##		These functions estimate uTau, y-offset, kappa, and a
##
###############################################################################
#
##	import standard modules
import numpy as np
import pandas as pd
import scipy.optimize as opt
#
###############################################################################
###	Theoretical velocity profile - composite method of Monkewitz et al (2007)
###	method taken from paper: Criteria for assessing experiments in zero
####					 pressure gradient boundary layers
##
def compositeProfile(UPlus,yPlus,deltaPlus,kappa,nu,a,PI,method):
	A = -1.0*a
	alpha = 0.5*(-1.0/kappa-A)
	beta = np.sqrt(-2.0*alpha*A - alpha**2)
	R = np.sqrt(alpha**2 + beta**2)
#
	inner = 1.0/kappa*np.log((yPlus - A)/(-1.0*A)) + R**2/(A*(4*alpha - A))*(
		(4*alpha + A)*np.log(-1.0*A*np.sqrt((yPlus - alpha)**2+beta**2)/(R*(yPlus - A)))
		+ (alpha/beta)*(4*alpha + 5*A)*(np.arctan((yPlus-alpha)/beta) + np.arctan(alpha/beta))
	) + np.exp(-np.log(yPlus/30)**2)/2.85
#
	if method == 'none':
		U = inner
	else:
		eta = yPlus/deltaPlus
		if method == 'exp':
			a2 = 132.8410
			a3 = -166.2041
			a4 = 71.9114
#			eta[eta>1.0]=1.0
			outer = 	(2.0*PI/kappa)*(1 - 0.5*np.log(eta)/PI)*(
						1 - np.exp(
							-0.25*(5*a2 + 6*a3 + 7*a4)*eta**4 
							+ a2*eta**5 + a3*eta**6 + a4*eta**7
						)
					)/(1 - np.exp(-0.25*(a2 + 2*a3 + 3*a4)))
		elif method == 'sine':
			outer = 2*PI*np.sin(np.pi*eta/2)**2/kappa
		elif method == 'channel':
			c2 = -20.22
			c3 = 17.1
			c4 = -11.17
#			eta[eta>1.0]=1.0
			outer = 	(2.0*PI/kappa)*(1 - 0.5*np.log(eta)/PI)*(
						1 - np.exp(
							eta**2*(
								c2*(eta-1.5)
						          + c3*(eta**4 - 3.0)
							    + c4*(eta**5 - 3.5)
							)
						)
					)/(1 - np.exp(-0.5*(c2 + 4.0*c3 + 5.0*c4)))
		U = inner + outer #+ np.log(yPlus)/kappa + 6.0
		U[eta>1] = U[eta<=1][-1]

	
	return(U)
###############################################################################
####	Least squares minimisation method based on above composite profile definition
####	NOTE: there is scope to add a weighting factor to the rms calculations
def compositeLeastSquares(UPlus,yPlus,deltaPlus,kappa,nu,a,PI,method):
	Utheory = compositeProfile(UPlus,yPlus,deltaPlus,kappa,nu,a,PI,method)
	cutoff = 0.2
	rms = np.sqrt(
		np.mean(
			((Utheory[(yPlus<cutoff*deltaPlus) & (yPlus>0.001)]
		    		-UPlus[(yPlus<cutoff*deltaPlus) & (yPlus>0.001)]
			)**2))
	) 
	nrms = rms/np.mean(UPlus[(yPlus<cutoff*deltaPlus) & (yPlus>0.001)])
	return(nrms)
###############################################################################
#####		compFit:	Input initial conditions and a method
import matplotlib.pyplot as plt
def compFit(x0,uTilde, yTilde, nu, PI, method, flag=False):
	UPlus = uTilde/x0[0]
	yPlus = (yTilde - x0[1])*x0[0]/nu
	temp = opt.minimize(deltaFunc,1.0,args=(yPlus,UPlus), method='Nelder-Mead',tol=1e-3,
			 			options={'disp':False,'xtol':1e-3,'ftol':1e-3})
	coeffs = deltaFunc(temp.x,yPlus,UPlus,flag=True)
#	print(coeffs)
	U99 = 0.99*max(UPlus)
#	print('cont')
	temp2 = opt.minimize(deltaFunc2,200,args=(temp.x,coeffs,U99), method='Nelder-Mead',tol=1e-3,
						options={'disp':False,'xtol':1e-3,'ftol':1e-3})
#	delta = 37e-3
	deltaPlus = temp2.x
#	delta2 = yPlus[UPlus>0.98*np.mean(UPlus[-4:])][0]
#	Ufit = temp.x*np.log(yPlus) + coeffs[0]*yPlus**3 + coeffs[1]*yPlus**2 + coeffs[2]*yPlus + coeffs[3]
#	plt.plot(yPlus,UPlus,linestyle='',marker='.')
#	plt.plot(yPlus,Ufit)
#	plt.scatter(deltaPlus,U99,marker='x',s=40)
#	plt.scatter(delta2,0.98*np.mean(UPlus[-4:]),marker='*',s=40)
#	plt.xlabel('y+')
#	plt.ylabel('U+')
#	plt.savefig('wakefit.png')
#	plt.close()

	delta = deltaPlus*nu/x0[0]
#	deltaPlus = delta*x0[0]/nu
#	print(delta)
	if flag==False:
		return compositeLeastSquares(UPlus,yPlus,deltaPlus,x0[2],nu, x0[3], PI, method)
	else:
		newCoeffs = np.insert(coeffs,0,temp.x[0])
		return [delta,newCoeffs,UPlus,yPlus]
###############################################################################
#####		deltaFunc:	Input some params, returns delta
def deltaFunc(X0,yTilde,UTilde,lim=100,flag=False):
	y = yTilde[yTilde>lim]
	U = UTilde[yTilde>lim]
	coeffs = np.polyfit(y,U-X0*np.log(y),3)
	Ufit = X0*np.log(y) + coeffs[0]*y**3 + coeffs[1]*y**2 + coeffs[2]*y + coeffs[3]
#	print(Ufit)
	nrms = np.sqrt(
		np.mean((U-Ufit)**2)
	)/np.mean(Ufit)
	if flag==False:
		return nrms
	else:
		return coeffs
###############################################################################
#####		deltaFunc2:	Input some params, returns delta
def deltaFunc2(X0,coeff1,coeffs,U99,flag=False):
	eps = np.absolute(U99 - (coeff1*np.log(X0) + coeffs[0]*X0**3 + coeffs[1]*X0**2 + coeffs[2]*X0 + coeffs[3]))
#	print(eps)
	if flag==False:
		return eps
	else:
		return X0
###############################################################################	



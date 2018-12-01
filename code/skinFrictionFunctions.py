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
from pyevtk.hl import gridToVTK
#
##	Basic linear regression function: returns alpha, beta, and rSqrd, such that
##	X = beta*Y+alpha
def linearReg(X,Y):
	beta = np.sum( (X-np.mean(X)) * (Y-np.mean(Y)) ) / np.sum( (Y - np.mean(Y))**2.0 ) 
	alpha = np.mean(X) - beta*np.mean(Y)
	rSqrd = 1.0 - np.sum((X - alpha - beta*Y)**2.0)/np.sum((X - np.mean(X))**2.0)
	return [alpha,beta,rSqrd]

###
##	use linear regression to evaluate friction velocity from near wall
##	velocity measurements
##
##	U and y are numpy vectors, nu is float
def viscousSublayerEstimation(U,y,nu):
	yVis = y[y<5e-3]
	UVis = U[y<5e-3]
	[a,b,r2] = linearReg(UVis,yVis)
#	print(len(yVis))
	yPlus = (yVis + a/b)*np.sqrt(b*nu)/nu
	while (r2<0.99) & (max(yPlus)>6.0):
		if len(yVis)<5:
			print('Viscous sub-layer too thin!')
			return [np.nan,np.nan]
		yVis = yVis[:-1]
		UVis = UVis[:-1] 
#		print(len(yVis))
		[a,b,r2] = linearReg(UVis,yVis)
#		print(np.sqrt(b*nu))
#		plt.plot(yVis,yVis*b + a)
#		plt.plot(yVis,UVis,linestyle='',marker='x')
#		plt.show()
	print(r2)
	uTau = np.sqrt(b*nu)
	dy = a/b
	return [uTau,dy]
###############################################################################
##	Assumes offset has been taken from y vector (i.e yTrue)
def logLawEstimation(U,y,uTauEst,kappa,nu):
	yPlusEst = y*uTauEst/nu
#	print(U)
	tdelta = y[U<0.95*np.mean(U[-4:])][-1]
#	print(y)
#	print(tdelta)
#	print(yPlusEst)
	yLog = y[(yPlusEst > 30) & (y < 0.7*(tdelta))]
	ULog = U[(yPlusEst > 30) & (y < 0.7*(tdelta))]
	[a,b,e] = linearReg(ULog,np.log(yLog))
#	while (e>1e-5):
#		if len(yLog)<10:
#			break
#		yLog = yLog[:-1]
#		ULog = ULog[:-1]
#		[a,b,e] = linearReg(ULog,np.log(yLog))
#	print(uTauEst, b*kappa)
#		print(e)
#	plt.semilogx(y,np.log(y)*b + a)
#	plt.semilogx(y[(yPlusEst > 50) & (y < 0.7*(tdelta))],U[(yPlusEst > 50) & (y < 0.7*(tdelta))])
#	plt.semilogx(y,U,linestyle='',marker='x')
#	plt.show()
	return [a,b,e]


################################################################################
##	Clauser method
def clauserEstimation(U,y,kappa,nu):
	U0 = np.mean(U[-4:])
	deltaStar = np.trapz(1.0-np.insert(U,0,0)/U0, np.insert(y,0,0))
	theta = np.trapz(np.insert(U,0,0)*(1.0-np.insert(U,0,0)/U0)/U0,np.insert(y,0,0))
#	deltaStar = np.trapz(1.0-U/U0,y)
#	theta = np.trapz(U*(1.0-U/U0)/U0,y)
#	print(deltaStar,theta,deltaStar/theta)
	ReDeltaStar = U0*deltaStar/nu
#	print(ReDeltaStar, U0*theta/nu)
	uTau = U0/(
		np.log(ReDeltaStar)/kappa + 3.3 + 182.0*np.log(ReDeltaStar)/ReDeltaStar - 2466.0/ReDeltaStar
	)
#	print(uTau)
	return [uTau,deltaStar,theta, deltaStar/theta]

################################################################################
###	Theoretical velocity profile - composite method of Monkewitz et al (2007)
###	method taken from paper: Criteria for assessing experiments in zero
####					 pressure gradient boundary layers
##
def compositeProfile(UPlus,yPlus,kappa,nu,a,PI,method):
	A = -1.0*a
	alpha = 0.5*(-1.0/kappa-A)
	beta = np.sqrt(-2.0*alpha*A - alpha**2)
	R = np.sqrt(alpha**2 + beta**2)
#
	inner = 1.0/kappa*np.log((yPlus - A)/(-1.0*A)) + R**2/(A*(4*alpha - A))*(
		(4*alpha + A)*np.log(-1.0*A*np.sqrt((yPlus - alpha)**2+beta**2)/(R*(yPlus - A)))
		+ (alpha/beta)*(4*alpha + 5*A)*(np.arctan((yPlus-alpha)/beta) + np.arctan(alpha/beta))
	) + np.exp(-np.log(yPlus/30)**2)/2.85



	if method == 'none':
		U = inner
	else:
		tdelta = yPlus[UPlus>0.98*np.mean(UPlus[-4:])][0]
		eta = yPlus/tdelta
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
####	
####	Least squares minimisation method based on above composite profile definition
def compositeLeastSquares(UPlus,yPlus,kappa,nu,a,PI,method):
	Utheory = compositeProfile(UPlus,yPlus,kappa,nu,a,PI,method)
	tdelta = yPlus[UPlus>0.98*np.mean(UPlus[-4:])][0]*0.2
	ty = np.log(yPlus)
#	w = np.exp(-1.0*(ty - 1.0*np.mean(ty))**2/(2.0*np.std(ty)**2))/np.sqrt(2.0*np.pi*np.std(ty))
#	w = np.array([1]*len(ty))
#	plt.semilogx(yPlus,w)
#	plt.show()
#	print(w)
#	print(len(UPlus),len(yPlus),len(w),len(Utheory))
	rms = np.sqrt(
		np.mean(
			((Utheory[(yPlus<tdelta) & (yPlus>0.001)]
		    		-UPlus[(yPlus<tdelta) & (yPlus>0.001)]
			)**2))
	) 
	nrms = rms/np.mean(UPlus[(yPlus<tdelta) & (yPlus>0.001)])
	return(nrms)

def testFunc(x0,uTilde, yTilde, nu, PI, method):
	UPlus = uTilde/x0[0]
	yPlus = (yTilde - x0[1])*x0[0]/nu
	return compositeLeastSquares(UPlus,yPlus,x0[2],nu, x0[3], PI, method)

####
####		Error matrix minimisation 
####	This function reads in kappa, deltaY and a vectors and creates matrix
####	of errors. Function returns the values of kappa, deltaY and a that give
####	minimum of the full matrix
def minFunc(U,yTilde,kappa,deltaY,a,nu,PI,method,saveName):
	errorGrid = np.zeros([len(kappa),len(deltaY),len(a)])
	for k in range(len(kappa)):
		print(str('Matrix construction ' + str(k*100/len(kappa)) + ' % complete'))
		for dy in range(len(deltaY)):
			y = (yTilde - deltaY[dy])
			[uTau, deltaStar, theta, H] = clauserEstimation(U,y,kappa[k],nu)
			yPlus = y*uTau/nu
			UPlus = U/uTau
			for i in range(len(a)):
				e = compositeLeastSquares(UPlus,yPlus,kappa[k],nu,a[i],PI,'none')
				if np.isnan(e):
					errorGrid[k][dy][i] = 1.0
				else:
					errorGrid[k][dy][i] = e
	#			print(e)

	locmin = np.unravel_index(errorGrid.argmin(), errorGrid.shape)
	locmax = np.unravel_index(errorGrid.argmax(), errorGrid.shape)
	
	if saveName is not False:
		if len(kappa) == 1:
			knorm = np.array([0,1])
		else:
			knorm = (kappa-np.min(kappa))/(np.max(kappa)-np.min(kappa))
		if len(deltaY) == 1:
			dynorm = np.array([0,1])
		else:
			dynorm = (deltaY-np.min(deltaY))/(np.max(deltaY)-np.min(deltaY))
		if len(a) == 1:
			anorm = np.array([0,1])
		else:
			anorm = (a-np.min(a))/(np.max(a)-np.min(a))
		enorm = (errorGrid - errorGrid[locmin])/(errorGrid[locmax]- errorGrid[locmin])
#		print(enorm)
		gridToVTK(str('./' + saveName), knorm,dynorm,anorm, pointData={'enorm':enorm})
	#	gridToVTK(str('./' + saveName), kappa,deltaY,a, pointData={'error':errorGrid})


	print(str('rms error = ' + str(errorGrid[locmin])))
	return([kappa[locmin[0]],deltaY[locmin[1]],a[locmin[2]]])

###############################################################################################
##
##	

###
###
def Btheory(kappa):
	Brange = [2,10]
	eps1 = 1.6*(np.exp(0.1663*Brange[0])-1.0) - kappa*Brange[0]
	eps2 = 1.6*(np.exp(0.1663*Brange[1])-1.0) - kappa*Brange[1]
	Bguess = np.mean(Brange)
	epsNew = 1.6*(np.exp(0.1663*Bguess)-1.0) - kappa*Bguess
	epsNewNorm = epsNew/(kappa*Bguess)
	N = 0
	while (epsNewNorm > 1e-5) & (N<4):
		N = N + 1
		if (epsNew < 0.0) & (eps1 < 0.0):
			print(1)
			BrangeNew = [Bguess,Brange[1]]
			eps1 = epsNew
			Bguess = np.mean(Brange)
			epsNew = 1.6*(np.exp(0.1663*Bguess)-1.0) - kappa*Bguess
		elif (epsNew < 0.0) & (eps1 > 0.0):
			print(2)
			BrangeNew = [Brange[0],Bguess]
			eps2 = epsNew
			Bguess = np.mean(Brange)
			epsNew = 1.6*(np.exp(0.1663*Bguess)-1.0) - kappa*Bguess
		elif (epsNew > 0.0) & (eps1 < 0.0):
			print(3)
			BrangeNew = [Brange[0],Bguess]
			eps1 = epsNew
			Bguess = np.mean(Brange)
			epsNew = 1.6*(np.exp(0.1663*Bguess)-1.0) - kappa*Bguess
		else:
			print(4)
			BrangeNew = [Brange[0],Bguess]
			eps2 = epsNew
			Bguess = np.mean(Brange)
			epsNew = 1.6*(np.exp(0.1663*Bguess)-1.0) - kappa*Bguess
#
		epsNewNorm = epsNew/(kappa*Bguess)
		print(Brange)
		
		print(epsNewNorm)

################################################################################
###		Choi rough surface uTau estimate and deltaY estimate
def choiEstimation(dfBase, df, uTauV, dyV, nu, saveName):
	errorGrid = np.zeros([len(uTauV),len(dyV),2])
	Xbaset = dfBase.y.as_matrix()*dfBase.uTau.as_matrix()[0]/(np.mean(dfBase.UxMean.as_matrix()[-4:])*dfBase.deltaStar.as_matrix()[0])	
	Ybaset = (np.mean(dfBase.UxMean.as_matrix()[-4:]) - dfBase.UxMean.as_matrix())/dfBase.uTau.as_matrix()[0]
	
	Xbase = Xbaset[dfBase.y > 0.2*dfBase.delta[0]]
	Ybase = Ybaset[dfBase.y > 0.2*dfBase.delta[0]]

	# estimate delta params from clauser function
	# these are dependent on dy but not uTau
	# assume kappa same for both dataframes
	for dy in range(len(dyV)):
		y = df.z.as_matrix()*1e-3 - dyV[dy]
		[t, deltaStar, theta, t] = clauserEstimation(df.UxMean.as_matrix(),y,dfBase.kappa.as_matrix()[0],nu)
		delta = y[df.UxMean.as_matrix()>0.98*np.mean(df.UxMean.as_matrix()[-4:])][0]
		# loop through uTau and calculate RMS error
		for uTau in range(len(uTauV)):
			Xt = y*uTauV[uTau]/(np.mean(df.UxMean.as_matrix()[-4:])*deltaStar)
			Yt = (np.mean(df.UxMean.as_matrix()[-4:]) - df.UxMean.as_matrix())/uTauV[uTau]	
			X = Xt[y > 0.2*delta]
			Y = Yt[y > 0.2*delta]
			#	Calculate range of X
#			Xrange = [min(Xbase),max(Xbase)]
#			newY = np.interp(Xbase,X,Y)#,np.nan,np.nan)
#			plt.plot(Xbase,newY,'^')
#			plt.plot(Xbase,Ybase,'o')
#			plt.show()
#
#			nrmsY = np.sqrt(np.mean((Ybase-newY)**2))/np.mean(Ybase)

			coeffs = np.polyfit(Xbase, Ybase, 4)
			Yref = coeffs[0]*X**4 + coeffs[1]*X**3 + coeffs[2]*X**2 + coeffs[3]*X + coeffs[4]
			nrmsY = np.sqrt(np.mean((Yref-Y)**2))/np.mean(Y)
	
			errorGrid[uTau][dy][0] = nrmsY
			errorGrid[uTau][dy][1] = nrmsY

	locmin = np.unravel_index(errorGrid.argmin(), errorGrid.shape)
	locmax = np.unravel_index(errorGrid.argmax(), errorGrid.shape)

	print(np.shape(errorGrid[:,:,0]),np.shape(uTauV),np.shape(dyV))
	plt.colorbar() 
	plt.contour(dyV,uTauV,errorGrid[:,:,0], 140, cmap="RdBu_r")
#	plt.plot(uTauV,errorGrid[:,locmin[1]])
	plt.show()
	
	dynorm = (dyV-np.min(dyV))/(np.max(dyV)-np.min(dyV))
	uTaunorm = (uTauV-np.min(uTauV))/(np.max(uTauV)-np.min(uTauV))
	enorm = (errorGrid - errorGrid[locmin])/(errorGrid[locmax]- errorGrid[locmin])
#		print(enorm)
	gridToVTK(str('./' + saveName), uTaunorm,dynorm,np.array([0,1]), pointData={'enorm':enorm})

	y = df.z.as_matrix()*1e-3 - dyV[locmin[1]]
	[t, deltaStar, theta, t] = clauserEstimation(df.UxMean.as_matrix(),y,dfBase.kappa.as_matrix()[0],nu)

#	print(errorGrid[locmin])
#	print(uTauV[locmin[0]],dyV[locmin[1]])
	return([uTauV[locmin[0]],dyV[locmin[1]],deltaStar,theta])

#####################################################################################################
####	composite profile estimates: comparing denticle surface to flat plate constants
def compositeFit(dfBase,df,uTauV,dYV,aV, nu, saveName):
	errorGrid = np.zeros([len(uTauV),len(dYV),len(aV)])
	for i in range(len(uTauV)):
		for j in range(len(dYV)):
			for k in range(len(aV)):
				tUPlus = df.UxMean.as_matrix()/uTauV[i]
				tyPlus = (df.z.as_matrix()*1e-3 - dYV[j])*uTauV[i]/nu
				UPlus = tUPlus[tyPlus > 0.01]
				yPlus = tyPlus[tyPlus > 0.01]
				if min(UPlus) < 0:
					e = 1
				else:
					e = compositeLeastSquares(UPlus,yPlus,dfBase.kappa.as_matrix()[0],nu,aV[k],0.5,'none')
				errorGrid[i][j][k] = e

	errorGrid[np.isnan(errorGrid)] = 1
	errorGrid[errorGrid>1] = 1
	locmax = np.unravel_index(errorGrid.argmax(), errorGrid.shape)
	locmin = np.unravel_index(errorGrid.argmin(), errorGrid.shape)

	uTaunorm = (uTauV-np.min(uTauV))/(np.max(uTauV)-np.min(uTauV))
	dynorm = (dYV-np.min(dYV))/(np.max(dYV)-np.min(dYV))
	if np.isnan(dynorm[0]):
		dynorm = np.array([0,1])

	anorm = (aV - np.min(aV))/(np.max(aV) - np.min(aV))
	enorm = (errorGrid - errorGrid[locmin])/(errorGrid[locmax]- errorGrid[locmin])
	
#		print(enorm)
	gridToVTK(str('./' + saveName), uTaunorm,dynorm,anorm, pointData={'enorm':enorm})
	#	gridToVTK(str('./' + saveName), kappa,deltaY,a, pointData={'error':errorGrid})

	print(str('rms error = ' + str(errorGrid[locmin])))
	return([uTauV[locmin[0]],dYV[locmin[1]],aV[locmin[2]]])

###########################################################################################################
###	wake fitting - NOT CHOI - we look at mean flow variables
def wakeFitting(dfBase,df,dyV):
	e = np.zeros(len(dyV))
##	original fitting
	Xbaset = dfBase.y.as_matrix()/dfBase.delta.as_matrix()[0]	
	Ybaset = (np.mean(dfBase.UxMean.as_matrix()[-4:]) - dfBase.UxMean.as_matrix())/(np.mean(dfBase.UxMean.as_matrix()[-4:])*dfBase.deltaStar.as_matrix()[0]/dfBase.delta.as_matrix()[0])
##	modified fitting to remove dependence on delta
#	Xbaset = dfBase.y.as_matrix()/dfBase.delta.as_matrix()[0]	
#	Ybaset = (np.mean(dfBase.UxMean.as_matrix()[-4:]) - dfBase.UxMean.as_matrix())/np.mean(dfBase.UxMean.as_matrix()[-4:])
	
	Xbase = Xbaset[dfBase.y > 0.2*dfBase.delta[0]]
	Ybase = Ybaset[dfBase.y > 0.2*dfBase.delta[0]]

#	Calculate range of Xbase for interpolation
	XbaseRange = [min(Xbase),max(Xbase)]


	# estimate delta params from clauser function
	# these are dependent on dy but not uTau
	# assume kappa same for both dataframes
	for dy in range(len(dyV)):
		y = df.z.as_matrix()*1e-3 - dyV[dy]
#		Boundary layer thickness of new data - needed to limit analysis to outer region
		delta = y[df.UxMean.as_matrix()>0.98*np.mean(df.UxMean.as_matrix()[-4:])][0]
#		print(delta, dfBase.delta[0])
		#	note - clauser function only used for integral quantities
		[t, deltaStar, theta, t] = clauserEstimation(df.UxMean.as_matrix(),y,1,1)
##		original
		Xt = y/delta
		Yt = (np.mean(df.UxMean.as_matrix()[-4:]) - df.UxMean.as_matrix())/(np.mean(df.UxMean.as_matrix()[-4:])*deltaStar/delta)
##		modified fit
#		Xt = y/delta
#		Yt = (np.mean(df.UxMean.as_matrix()[-4:]) - df.UxMean.as_matrix())/np.mean(df.UxMean.as_matrix()[-4:])
		X = Xt[y > 0.2*delta]
		Y = Yt[y > 0.2*delta]

		###	first method: interpolate denticle data to flat plate data and find rms error:
#		newY = np.interp(Xbase,X,Y)
#		nrmsY = np.sqrt(np.mean((Ybase-newY)**2))/np.mean(Ybase)

		### 	Second method: fit polynomial to flat plate data and use this to calculate 'ref' values for rms
		coeffs = np.polyfit(Xbase, Ybase, 4)
		Yref = coeffs[0]*X**4 + coeffs[1]*X**3 + coeffs[2]*X**2 + coeffs[3]*X + coeffs[4]
		nrmsY = np.sqrt(np.mean((Yref-Y)**2))/np.mean(Y)

		e[dy] = nrmsY
#		plt.plot(X,Y,'x')
#		plt.plot(X,Yref,'-')
#		plt.plot(Xbase,Ybase,'o')
#		plt.show()
#	plt.plot(dyV,e,'o')
#	plt.show()
	return(dyV[e.argmin()])

###########################################################################################################
###	loglaw fit - we fit denticle surface to the log-region of flat plate data
def logLawFitting(dfBase,df,uTauV,dUPlusV,nu,dyV,saveName):
	errorGrid = np.zeros([len(uTauV),len(dyV),len(dUPlusV)])
	### estimate kappa and C from base profile
	# initialise fields
	tUPlusBase = dfBase.UxMean.as_matrix()/dfBase.uTau.as_matrix()[0]
	tyPlusBase = dfBase.y.as_matrix()*dfBase.uTau.as_matrix()[0]/dfBase.nu.as_matrix()[0]
	# limit to log-region
	UPlusBase = tUPlusBase[(tyPlusBase > 50) & (dfBase.y < 0.8*dfBase.delta.as_matrix()[0])]
	yPlusBase = tyPlusBase[(tyPlusBase > 50) & (dfBase.y < 0.8*dfBase.delta.as_matrix()[0])]
	# find alpha and beta from linearRegression
	[a,b,R] = linearReg(UPlusBase,np.log(yPlusBase))
	# use these to estimate denticle params uTau and deltaU
	#
	#	i.e denticle case: u+ = 1/kappa * ln(y+) + C + deltaU
	#
	tdelta = df.z.as_matrix()[df.UxMean.as_matrix()>0.98*np.mean(df.UxMean.as_matrix()[-4:])][0]*1e-3
	tU = df.UxMean.as_matrix()
	for k in range(len(dyV)):
		ty = df.z.as_matrix()*1e-3 - dyV[k]
		for i in range(len(uTauV)):
			# construct fields	
			tUPlus = tU/uTauV[i]
			tyPlus = ty*uTauV[i]/nu
			UPlus = tUPlus[(tyPlus > 50) & (df.z.as_matrix()*1e-3 < 0.8*tdelta)]
			yPlus = tyPlus[(tyPlus > 50) & (df.z.as_matrix()*1e-3 < 0.8*tdelta)]
			for j in range(len(dUPlusV)):
				Ufit = b*np.log(yPlus) + a
				nrms = np.sqrt(np.mean(Ufit - UPlus - dUPlusV[j])**2)/np.mean(Ufit)	
				errorGrid[i][k][j] = nrms
	
	errorGrid[np.isnan(errorGrid)] = 1
	errorGrid[errorGrid>1] = 1
	locmax = np.unravel_index(errorGrid.argmax(), errorGrid.shape)
	locmin = np.unravel_index(errorGrid.argmin(), errorGrid.shape)

	uTaunorm = (uTauV-np.min(uTauV))/(np.max(uTauV)-np.min(uTauV))
#	dynorm = (dYV-np.min(dYV))/(np.max(dYV)-np.min(dYV))
	dUPlusnorm = (dUPlusV - np.min(dUPlusV))/(np.max(dUPlusV) - np.min(dUPlusV))
	enorm = (errorGrid - errorGrid[locmin])/(errorGrid[locmax]- errorGrid[locmin])
	
#		print(enorm)
	gridToVTK(str('./' + saveName), uTaunorm,np.array([0,1]),dUPlusnorm, pointData={'enorm':enorm})
	#	gridToVTK(str('./' + saveName), kappa,deltaY,a, pointData={'error':errorGrid})

	print(str('rms error = ' + str(errorGrid[locmin])))
	return([uTauV[locmin[0]],dyV[locmin[1]],dUPlusV[locmin[2]]])
#	plt.semilogx(yPlusBase,UPlusBase,'o')
#	plt.plot(yPlusBase,b*np.log(yPlusBase) + a)
#	plt.show()
#	print(a,b,R)
		
	
##
##
##		NOTE: Does it matter if we can't predict uTau properly for flat surface?
##			Do we still get a correct relative error? 
#
################################################################################
#####	function for estimating error between composite profile and data
def genericPolyFit(X0,yTilde,U,lim,flag=False):
	y = yTilde[yTilde>lim] + X0[0]
	coeffs = np.polyfit(y,U[yTilde>lim]-X0[1]*np.log(y),3)
	Ufit = X0[1]*np.log(y) + coeffs[0]*y**3 + coeffs[1]*y**2 + coeffs[2]*y + coeffs[3]
#	print(Ufit)
	nrms = np.sqrt(
		np.mean((U[yTilde>lim]-Ufit)**2)
	)/np.mean(Ufit)
	if flag==False:
		return nrms
	else:
		return Ufit
################################################################################								

################################################################################
####
####		Script that compares probe configurations on freestream stats
####		and prints results to a .csv in latex format
####
################################################################################
import pandas as pd
import csv
from decimal import Decimal, ROUND_UP
import numpy as np
################################################################################
##
##	dfNames
raw =[
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation0/angle2p7/8hz_rotation0_dip_averaged_raw.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation0/angle0/8hz_rotation0_nodip_averaged_raw.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation45/angle2p7/8hz_rotation45_dip_averaged_raw.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation45/angle0/8hz_rotation45_nodip_averaged_raw.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation90/angle2p7/8hz_rotation90_dip_averaged_raw.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation90/angle0/8hz_rotation90_nodip_averaged_raw.pkl')
	]
fil =[
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation0/angle2p7/8hz_rotation0_dip_averaged_w50_MA_high.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation0/angle0/8hz_rotation0_nodip_averaged_w50_MA_high.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation45/angle2p7/8hz_rotation45_dip_averaged_w50_MA_high.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation45/angle0/8hz_rotation45_nodip_averaged_w50_MA_high.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation90/angle2p7/8hz_rotation90_dip_averaged_w50_MA_high.pkl'),
	pd.read_pickle('../data/processedData/dataQualityTests/8Hz/rotation90/angle0/8hz_rotation90_nodip_averaged_w50_MA_high.pkl')
	]

rotationAngles = ['0','0','45','45','90','90']
dippingAngles = ['2.7','0','2.7','0','2.7','0']
tols = ['0.001','0.001','0.01','0.01','0.001']

def rounder(numberString,toleranceString):
	return str(Decimal(numberString).quantize(Decimal(toleranceString), rounding=ROUND_UP))

with open('test.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile, delimiter='&')
	writer.writerow(['\\begin{tabular}{l l r r r r r r r r r r}'])
	writer.writerow(	['Rot.\t']
			 + 	['\tDip.\t']
			 + 	['\t$\mu_u$\t']
			 + 	['\t$\overline{\mu_u}$\t']
			 + 	['\t$\mu_v$\t']
			 + 	['\t$\overline{\mu_v}$\t']
			 + 	['\t$\sigma_u$\t']
			 + 	['\t$\overline{\sigma_u}$\t']
			 + 	['\t$\sigma_v$\t']
			 + 	['\t$\overline{\sigma_v}$\t']
			 + 	['\t$\gamma_{u,v}$\t']
			 + 	['\t$\overline{\gamma_{u,v}}$\t \\\\'])
	writer.writerow(	['angle\t']
			 + 	['\tangle\t']
			 + 	['\t$ $\t']
			 + 	['\t$ $\t']
			 + 	['\t$\\times 10^3$\t']
			 + 	['\t$\\times 10^3$\t']
			 + 	['\t$\\times 10^3$\t']
			 + 	['\t$\\times 10^3$\t']
			 + 	['\t$\\times 10^3$\t']
			 + 	['\t$\\times 10^3$\t']
			 + 	['\t$\\times 10^3$\t']
			 + 	['\t$\\times 10^3$\t \\\\'])
	writer.writerow(	['(Deg)\t']
			 + 	['\t(Deg)\t']
			 + 	['\t m/s \t']
			 + 	['\t m/s \t']
			 + 	['\t m/s \t']
			 + 	['\t m/s \t']
			 + 	['\t m/s \t']
			 + 	['\t m/s \t']
			 + 	['\t m/s \t']
			 + 	['\t m/s \t']
			 + 	['\t m\\textsuperscript{2}/s\\textsuperscript{2} \t']
			 + 	['\t m\\textsuperscript{2}/s\\textsuperscript{2} \t \\\\'])
	writer.writerow(['\hline'])
#
##	Now write data - loop through names, corresponding to rows
	for i in range(len(raw)):
		writer.writerow(
			[str('\t' + rotationAngles[i] + '\t')]
		+	[str('\t' + dippingAngles[i] + '\t')]
		+	[str('\t' + rounder(str(raw[i].loc[raw[i]['z']==5.0]['UxMean'].as_matrix()[0]),'0.001') + '\t')]
		+	[str('\t' + rounder(str(fil[i].loc[fil[i]['z']==5.0]['UxMean'].as_matrix()[0]),'0.001') + '\t')]
		+	[str('\t' + rounder(str(raw[i].loc[raw[i]['z']==5.0]['UyMean'].as_matrix()[0]*1000),'0.001') + '\t')]
		+	[str('\t' + rounder(str(fil[i].loc[fil[i]['z']==5.0]['UyMean'].as_matrix()[0]*1000),'0.001') + '\t')]
		+	[str('\t' + rounder(str(raw[i].loc[raw[i]['z']==5.0]['uxRMS'].as_matrix()[0]*1000),'0.01') + '\t')]
		+	[str('\t' + rounder(str(fil[i].loc[fil[i]['z']==5.0]['uxRMS'].as_matrix()[0]*1000),'0.01') + '\t')]
		+	[str('\t' + rounder(str(raw[i].loc[raw[i]['z']==5.0]['uyRMS'].as_matrix()[0]*1000),'0.01') + '\t')]
		+	[str('\t' + rounder(str(fil[i].loc[fil[i]['z']==5.0]['uyRMS'].as_matrix()[0]*1000),'0.01') + '\t')]
		+	[str('\t' + rounder(str(raw[i].loc[raw[i]['z']==5.0]['uv'].as_matrix()[0]*1000),'0.001') + '\t')]
		+	[str('\t' + rounder(str(fil[i].loc[fil[i]['z']==5.0]['uv'].as_matrix()[0]*1000),'0.001') + '\t \\\\')]
				) 
	writer.writerow(['\end{tabular}'])	

tmeanU, tmeanV, tstdU, tstdV, tcov = [],[],[],[],[]
for i in range(len(fil)):
	tmeanU.append(fil[i].loc[fil[i]['z']==5.0]['UxMean'].as_matrix())
	tmeanV.append(fil[i].loc[fil[i]['z']==5.0]['UyMean'].as_matrix())
	tstdU.append(fil[i].loc[fil[i]['z']==5.0]['uxRMS'].as_matrix())
	tstdV.append(fil[i].loc[fil[i]['z']==5.0]['uyRMS'].as_matrix())
	tcov.append(fil[i].loc[fil[i]['z']==5.0]['uv'].as_matrix())

meanU, meanV, stdU, stdV, cov = np.mean(tmeanU), np.mean(tmeanV), np.mean(tstdU), np.mean(tstdV), np.mean(tcov)
print(len(fil),meanU, meanV, stdU, stdV, cov)


###########################################################################################
##############		DATA QUALITY ANALYSIS SCRIPT		###############################
####
####		 
####
####
####
###########################################################################################
####
##		Initialise python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from movingAverageFilterReynoldsStresses import movingWeightedAverage as mwa
#
###########################################################################################
SA0_raw = pd.read_pickle('../data/processedData/dataQualityTests/8Hz/angle0/seeding/8hz_seeding_angle0_averaged_raw.pkl')
NSA0_raw = pd.read_pickle('../data/processedData/dataQualityTests/8Hz/angle0/noSeeding/8hz_noSeeding_angle0_averaged_raw.pkl')
SA2_raw = pd.read_pickle('../data/processedData/dataQualityTests/8Hz/angle2p7/seeding/8hz_seeding_angle2p7_averaged_raw.pkl')
NSA2_raw = pd.read_pickle('../data/processedData/dataQualityTests/8Hz/angle2p7/noSeeding/8hz_noSeeding_angle2p7_averaged_raw.pkl')

SA0_fil = pd.read_pickle('../data/processedData/dataQualityTests/8Hz/angle0/seeding/8hz_seeding_angle0_averaged_w50_MA_high.pkl')
NSA0_fil = pd.read_pickle('../data/processedData/dataQualityTests/8Hz/angle0/noSeeding/8hz_noSeeding_angle0_averaged_w50_MA_high.pkl')
SA2_fil = pd.read_pickle('../data/processedData/dataQualityTests/8Hz/angle2p7/seeding/8hz_seeding_angle2p7_averaged_w50_MA_high.pkl')
NSA2_fil = pd.read_pickle('../data/processedData/dataQualityTests/8Hz/angle2p7/noSeeding/8hz_noSeeding_angle2p7_averaged_w50_MA_high.pkl')

######	
##	Plot Ux
#plt.plot()


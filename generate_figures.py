# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-31 18:46:04
# @Last Modified by:   joaopn
# @Last Modified time: 2019-03-31 22:52:06

import analysis
import matplotlib.pyplot as plt

#Defines state dict
state_dict = {
	'poisson': {'m': 0.000, 'h': 2e-3},
	'subcritical': {'m': 0.900, 'h': 2e-4},
	'reverberant': {'m': 0.980, 'h': 4e-5},
	'critical': {'m': 0.999, 'h': 2e-6},
}


#Generates figure 1

#Sets path
str_savepath = 'figs/fig1/'
if not os.path.exists(str_savepath):
	os.makedirs(str_savepath)

#Generates timeseries figures
plt.figure()




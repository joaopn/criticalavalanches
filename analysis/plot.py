# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-26 13:40:21
# @Last Modified by:   joaopn
# @Last Modified time: 2019-03-27 11:02:55

import matplotlib.pyplot as plt
import numpy as np
import analysis.avalanche

def pS(S,label='data'):

	#Calculates p(S) (no log-binning)
	S_max = int(S.max())
	pS = np.zeros(S_max)

	for i in range(1,S_max):
		pS[i] = np.sum(S==i)

	pS = pS/np.sum(pS)

	#Plots it
	plt.loglog(pS,label=label)
	plt.xlim(1,1e3)
	plt.legend()


def timeseries_threshold(data,th):
	data_th = analysis.avalanche.threshold_data(data,th)

	X = np.arange(data.size)

	plt.plot(data)
	plt.plot(X,th*np.std(data)*np.ones(data.size))
	plt.scatter(X[data_th==1],data[data_th==1])



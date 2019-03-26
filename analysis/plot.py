# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-26 13:40:21
# @Last Modified by:   joaopn
# @Last Modified time: 2019-03-26 13:49:36

import matplotlib.pyplot as plt
import numpy as np

def pS(S,label='data'):

	#Calculates p(S) (no log-binning)
	S_max = int(S.max())
	pS = np.zeros(S_max)

	for i in range(1,S_max):
		pS[i] = np.sum(S==i)

	pS = pS/np.sum(pS)

	#Plots it
	plt.figure()
	plt.loglog(pS,label=label)
	plt.xlim(1,S_max)
	plt.legend()
	plt.xlabel('S')
	plt.ylabel('p(S)')
	plt.show()


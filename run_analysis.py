# -*- coding: utf-8 -*-
# @Author: Joao PN
# @Date:   2019-03-25 16:45:25
# @Last Modified by:   joaopn
# @Last Modified time: 2019-03-26 13:57:11

from analysis import *
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

	#Runs analysis for one timeseries as plots it
	data_binned = avalanche.run_analysis(filepath='00050.hdf5',binsize=4)

	S = avalanche.get_S(data_binned)

	plot.pS(S)

# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-24 19:03:08
# @Last Modified by:   joaopn
# @Last Modified time: 2019-03-25 10:08:12

import analysis.utils as utils
import analysis.avalanche as avalanche
import analysis.filtering as filtering
import analysis.fitting as fitting

import powerlaw
import numpy as np

def run_analysis(
	filepath=None,
	threshold=3,
	channels=100,
	timesteps = 500000,
	binsize = 4
	):


	#Analyses channel by channel
	data_th = np.zeros(timesteps)
	for ch in range(channels):
		
		#Loads data
		data_ch = utils.load_hdf5(filepath,'data/coarse',ch)

		#Thresholds data		
		data_th = data_th + avalanche.threshold_data(data_ch,threshold)


	#Plots distribution
	bin_array = np.arange(0,timesteps/binsize)
	data_pdf = np.histogram(data,bin=bin_array)[0]
	plt.plot(data_pdf)

	#data_analyzed = utils.bin_data(data_th,binsize)
	#print(data_analyzed.shape)

	#Fits and plots power-laws
	#results = powerlaw.Fit(data_analyzed)
	#print(results.power_law.alpha)
	#print(results.power_law.xmin)
	#R, p = results.distribution_compare('power_law', 'lognormal')
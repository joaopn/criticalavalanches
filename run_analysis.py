# -*- coding: utf-8 -*-
# @Author: Joao PN
# @Date:   2019-03-25 16:45:25
# @Last Modified by:   Joao PN
# @Last Modified time: 2019-03-25 20:59:57

from analysis import *
import numpy as np
import matplotlib.pyplot as plt

def run_analysis(
	filepath='',
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

	return data_th



	#data_analyzed = utils.bin_data(data_th,binsize)
	#print(data_analyzed.shape)

	#Fits and plots power-laws
	#results = powerlaw.Fit(data_analyzed)
	#print(results.power_law.alpha)
	#print(results.power_law.xmin)
	#R, p = results.distribution_compare('power_law', 'lognormal')

if __name__ == "__main__":
	data_th = run_analysis(filepath='00050.hdf5')

	#Plots distribution
	bin_array = np.arange(0,timesteps/binsize)
	data_pdf = np.histogram(data_th,bins=bin_array)[0]
	plt.loglog(data_pdf)
	plt.xlim(1,100)
	plt.show()
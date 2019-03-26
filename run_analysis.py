# -*- coding: utf-8 -*-
# @Author: Joao PN
# @Date:   2019-03-25 16:45:25
# @Last Modified by:   Joao PN
# @Last Modified time: 2019-03-26 19:54:32

from analysis import *
import numpy as np
import matplotlib.pyplot as plt
import os


def run_save_plot(filename,binsize, datatype):

	#Parameters
	filepath = 'dat/'
	fig_dir = 'dat/plot/'

	#Plots results
	plt.figure()
	plt.title(filename)
	plt.xlabel('S')
	plt.ylabel('p(S)')

	#Runs analysis for all binsizes
	for bs in binsize:
		print('Running for b = {}'.format(bs))
		data_binned = avalanche.run_analysis(
			filepath=filepath+filename+'.hdf5',
			binsize=bs,
			threshold=threshold
			)
		S = avalanche.get_S(data_binned)
		plot.pS(S,label='bs = {}'.format(bs))

	#Creates save dir
	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	#Saves fig
	str_fig = fig_dir + datatype + '_pS_' + filename + '_th{}'.format(threshold) + '.png'
	plt.savefig(str_fig)

if __name__ == "__main__":

	#Parameters
	binsize=[4]
	threshold = 3.5
	reps = 3
	datatype = 'coarse'

	#Test data
	state_dict = {
		'subcrit' : {'m' : 0.92,'h' : 1.8e-4},
		'rev' : {'m' : 1.006,'h' : 3e-5},
		'crit' : {'m' : 1.03,'h' : 2.2e-6},
	}


	#Analyzes and saves plot for all datasets
	for state in state_dict:
		for rep in range(reps):

			#Defines variables
			m = state_dict[state]['m']
			h = state_dict[state]['h']
			filename = 'm{:0.5f}_h{:0.3e}_r{:d}'.format(m,h,rep)
		
			run_save_plot(filename,binsize,datatype)
	
	

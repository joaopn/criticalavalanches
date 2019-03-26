# -*- coding: utf-8 -*-
# @Author: Joao PN
# @Date:   2019-03-25 16:45:25
# @Last Modified by:   joaopn
# @Last Modified time: 2019-03-26 15:01:35

from analysis import *
import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":

	#Parameters
	filepath = ''
	filename='00050'
	fig_dir = 'dat/plot/'
	binsize=[4,8,16]
	threshold = 3

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
	str_fig = fig_dir + 'pS_' + filename + '_th{}'.format(threshold) + '.png'
	plt.savefig(str_fig)
	
	
	

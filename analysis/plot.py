# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-26 13:40:21
# @Last Modified by:   joaopn
# @Last Modified time: 2019-03-31 23:09:15

import os
import matplotlib
if os.environ.get('DISPLAY', '') == '':
    matplotlib.use('Agg')
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

def pS_mean(S_list,label='data'):

	#Gets largest avalanche from the list (+1 for zero_index)
	S_max = int(max([Si.max() for Si in S_list]) + 1)

	#Obtains pS
	pS = np.zeros((len(S_list),S_max))

	for i in range(len(S_list)):
		for j in range(S_max):
			pS[i,j] = np.sum(S_list[i]==j)
		pS[i,:] = pS[i,:]/np.sum(pS[i,:])

	#Obtains mean and STD
	pS_mean = np.mean(pS,axis=0)
	pS_std = np.std(pS,axis=0)
	pS_up = pS_mean + pS_std/2
	pS_dw = pS_mean - pS_std/2

	#Plots confidence interval (1 std) and mean
	plt.fill_between(range(S_max),pS_up,pS_dw,alpha=0.75)
	plt.plot(range(S_max),pS_mean,label=label)
	plt.legend()

def timeseries_threshold(data,th):
	data_th = analysis.avalanche.threshold_data(data,th)

	X = np.arange(data.size)

	plt.plot(data)
	plt.plot(X,th*np.std(data)*np.ones(data.size))
	plt.scatter(X[data_th==1],data[data_th==1])

def sim_pS(m,h,d,b,datatype,reps,label_plot=None,bw_filter = False,data_dir ='dat/',threshold = 3):

	#Definitions
	if bw_filter:
		saveplot_dir = 'analyzed_filtered/'
	else:
		saveplot_dir = 'analyzed_unfiltered/'


	#Parse input
	if type(d) != list:
		d = [d]
	if type(b) != list:
		b = [b]

	for i in range(len(d)):

		#Sets up filepath and loads data
		dataset = 'm{:0.5f}_h{:0.3e}_d{:02d}_th{:0.1f}_rep{:02d}/'.format(m,h,d[i],threshold,reps)
		datafile = 'pS_' + datatype + '_b{:02d}.tsv'.format(b[i])
		str_load = data_dir + saveplot_dir + dataset + datafile

		S,pS_mean, pS_std = np.loadtxt(str_load,delimiter='\t')

		#Plots data
		if label_plot is None:
			label_i = 'm = {:0.3f}, b = {:d}, d = {:d}'.format(m, b[i], d[i])
		else:
			label_i = label_plot[i]

		pS_up = pS_mean + pS_std/2
		pS_dw = pS_mean - pS_std/2
		plt.fill_between(S,pS_up,pS_dw,alpha=0.75)
		plt.loglog(S,pS_mean,label=label_i)
		
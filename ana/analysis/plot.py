# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-26 13:40:21
# @Last Modified by:   Joao PN
# @Last Modified time: 2019-04-03 17:43:49

import os
import matplotlib
#if os.environ.get('DISPLAY', '') == '':
#    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import analysis.avalanche
import analysis.fitting

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
	plt.fill_between(range(S_max),pS_up,pS_dw,alpha=0.25)
	plt.plot(range(S_max),pS_mean,label=label)
	plt.legend()

def timeseries_threshold(data,th):
	data_th = analysis.avalanche.threshold_data(data,th)

	X = np.arange(data.size)

	plt.plot(data)
	plt.plot(X,th*np.std(data)*np.ones(data.size))
	plt.scatter(X[data_th==1],data[data_th==1])

def sim_pS(m,h,d,b,datatype,reps,label_plot=None,bw_filter = False,data_dir ='dat/',threshold = 3, plt_color='black'):

	#Parameters
	S_lim = 1e-8 #Fixes vectorial plotting bug

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

		#Fixes vector bug
		S_rem = np.argwhere(pS_mean < S_lim)
		pS_mean = np.delete(pS_mean,S_rem)
		pS_std = np.delete(pS_std,S_rem)
		S = np.delete(S,S_rem)

		#Plots data
		pS_up = pS_mean + pS_std/2
		pS_dw = pS_mean - pS_std/2
		plt.fill_between(S,pS_up,pS_dw,alpha=0.5,color=plt_color)

		if label_plot is None:
			tau_rep = analysis.fitting.tau_sim_dataset(m,h,d[i],threshold,data_dir,bw_filter)
			plt.loglog(
				S,
				pS_mean,
				label= r'$\tau$ = {:0.1f} $\pm$ {:0.1f} ms'.format(np.mean(tau_rep), np.std(tau_rep)),
				color=plt_color)
		else:
			label_i = label_plot[i]
			plt.loglog(
				S,
				pS_mean,
				label=label_plot[i],
				color=plt_color)

def analyze_pS(data, b, threshold=3):
	"""Does the avalanche analysis of a raw data matrix and plots p(S)
	
	[description]
	:param data: [n_ch,timeteps] raw data array
	:type data: float
	:param b: binsize to plot (in units of timesteps)
	:type b: float
	:param threshold: threshold used to discretize the signal, defaults to 3
	:type threshold: number, optional
	"""

	#Parameters
	n_e, timesteps = data.shape
	
	#Thresholds data
	data_th = np.zeros(timesteps)

	for ch in range(n_e):
		data_th_ch = analysis.avalanche.threshold_ch(data[ch,:], threshold)
		data_th = data_th + data_th_ch

	#Bins data
	data_binned = analysis.avalanche.bin_data(data_th,b)

	#Calculates S
	S = analysis.avalanche.get_S(data_binned)

	#Plots pS
	analysis.plot.pS(S,label="b = {:d}".format(b))


	return data_th, data_binned,S
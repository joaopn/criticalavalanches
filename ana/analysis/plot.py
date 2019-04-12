# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-26 13:40:21
# @Last Modified by:   joaopn
# @Last Modified time: 2019-04-12 03:03:39

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

def sim_pS(m,h,d,b,datatype,reps,label_plot=None,bw_filter = False,data_dir ='dat/',threshold = 3, plt_color='black', plt_std=True, offset=1):

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
	if type(label_plot) == str:
		label_plot = [label_plot]

	assert len(d) == len(b)

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
		if plt_std:
			pS_up = pS_mean + pS_std/2
			pS_dw = pS_mean - pS_std/2
			plt.fill_between(S,offset*pS_up,offset*pS_dw,alpha=0.5,color=plt_color)

		if label_plot is None:
			tau_rep = analysis.fitting.tau_sim_dataset(m,h,d[i],threshold,data_dir,bw_filter)
			str_label = r'$\tau$ = {:0.1f} $\pm$ {:0.1f} ms'.format(np.mean(tau_rep), np.std(tau_rep))
			plt.loglog(S,offset*pS_mean,label= str_label,color=plt_color)
		else:			
			plt.loglog(S,offset*pS_mean,label=label_plot[i],
				color=plt_color)

def sim_mav(m,h,b_list,data_dir,label_plot=None,bw_filter=False,threshold=3,plt_color='black', linestyle='-'):

	#Definitions
	if bw_filter:
		saveplot_dir = 'analyzed_filtered/branching_mav/'
	else:
		saveplot_dir = 'analyzed_unfiltered/branching_mav/'

	#Parse input
	if type(b_list) is not list:
		b_list = [b_list]

	#Plots for every b
	for b in b_list:
		#Sets filepath
		dataset = 'm{:0.5f}_h{:0.3e}_b{:02d}_th{:0.1f}.tsv'.format(m,h,b,threshold)
		str_load = data_dir + saveplot_dir + dataset

		#Loads data
		IED,mav_mean,mav_std = np.loadtxt(str_load,delimiter='\t')

		#Plots data
		mav_up = mav_mean + mav_std/2
		mav_dw = mav_mean - mav_std/2
		plt.fill_between(IED,mav_up,mav_dw,alpha=0.5,color=plt_color)

		if label_plot is None:

			print('Calculating tau for b = {:02d}'.format(b))

			tau_all = np.zeros(0)
			for d in IED:
				tau_d = analysis.fitting.tau_sim_dataset(m,h,int(d),threshold,data_dir,bw_filter)
				tau_all = np.concatenate((tau_all,tau_d))

			tau_mean = np.mean(tau_all)
			tau_std = np.std(tau_all)

			print('tau = {:0.1f} +- {:0.1f}'.format(tau_mean, tau_std))
			str_label = r'$\Delta$t = {:d} ms, $\tau$ = {:0.1f}$\pm${:0.1f} ms'.format(2*b,tau_mean, tau_std)

			plt.plot(IED,mav_mean,label= str_label,color=plt_color,linestyle=linestyle)
		else:
			plt.plot(IED,mav_mean,label=label_plot,color=plt_color,linestyle=linestyle)

def plot_alpha_bs(m,h,d,datatype,reps,bw_filter = False,data_dir ='dat/',threshold = 3, color_rgba=None):

	#Definitions
	if bw_filter:
		saveplot_dir = 'analyzed_filtered/'
	else:
		saveplot_dir = 'analyzed_unfiltered/'

	#Sets up filepath and loads data
	dataset = 'm{:0.5f}_h{:0.3e}_d{:02d}_th{:0.1f}_rep{:02d}/'.format(m,h,d,threshold,reps)
	datafile = 'alpha_' + datatype + '.tsv'
	str_load = data_dir + saveplot_dir + dataset + datafile	
	b, alpha_mean, alpha_std = np.loadtxt(str_load,delimiter='\t')

	#Fits data to alpha~b^-fit_exp
	fit_exp, fit_err, lin_coef = analysis.fitting.powerlaw(b,alpha_mean,alpha_std)

	#Plots fit
	X = np.arange(1,100)
	plt.plot(X,lin_coef*np.power(X,-fit_exp),linestyle='--', color='k')

	#Plots data
	if color_rgba is not None:
		plt.errorbar(b, alpha_mean, yerr=alpha_std/2, ecolor='k',fmt='none', elinewidth=1)
		plt.scatter(b,alpha_mean,marker='s',color=color_rgba)
		
	else:
		plt.errorbar(b, alpha_mean, yerr=alpha_std/2, ecolor='k', fmt='s')

	return fit_exp, fit_err


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

	#Parses input
	b = np.array(b)

	#Parameters
	n_e, timesteps = data.shape
	
	#Thresholds data
	data_th = np.zeros(timesteps)

	for ch in range(n_e):
		data_th_ch = analysis.avalanche.threshold_ch(data[ch,:], threshold)
		data_th = data_th + data_th_ch

	#Bins and plots data
	for i in range(len(b)):

		#Bins data
		data_binned = analysis.avalanche.bin_data(data_th,b[i])

		#Calculates S
		S = analysis.avalanche.get_S(data_binned)

		#Plots pS
		analysis.plot.pS(S,label="b = {:d}".format(b[i]))
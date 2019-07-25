# -*- coding: utf-8 -*-
# @Author: Joao
# @Date:   2019-07-05 17:56:44
# @Last Modified by:   joaopn
# @Last Modified time: 2019-07-25 15:22:01

"""
Module for directly handling datasets.
"""

from analysis import avalanche, fitting, plot, parser
import powerlaw
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, os, pickle

def sim_plot_pS(filepath,deltaT,datatype,str_leg='Data', threshold=3,bw_filter=True,timesteps=None,channels=None, save_fig=None,show_error=True,color='k', lineType='-', zorder=2):
	"""Plots the avalanche size distribution for a single hdf5 dataset, and a single deltaT. If [filepath] is a list, averages over datasets.
	
	Args:
	    filepath (str): Path to the .hdf5 dataset
	    deltaT (int): Binsize to use, *in timesteps*
		datatype (str): 'coarse', 'sub'
	    str_leg (str, optional): plot legend
	    threshold (int, optional): Threshold in standard deviations of the signal (for coarse) 
	    bw_filter (bool, optional): Toggles butterworth filtering (for coarse)
	    timesteps (None, optional): Length of the dataset in timesteps
	    channels (None, optional): Number of channels in the dataset
	    color (str, optional): plot color
	"""

	if type(filepath) is not list:
		filepath = [filepath]
	nreps = len(filepath)


	S_list = []
	for filepath_file in filepath:
		#Loads and thresholds data
		data_th = avalanche.analyze_sim_raw(filepath_file, threshold, datatype, bw_filter, timesteps, channels)

		#Bins data
		data_binned = avalanche.bin_data(data=data_th,binsize=deltaT)

		#Gets S and plots it
		S_list.append(avalanche.get_S(data_binned))
	plot.pS_mean(S_list,label=str_leg,lineType=lineType,color=color,show_error=show_error, zorder=zorder)    

	#Sets up figure
	fig = plt.gcf()
	ax = plt.gca()
	ax.set_xlim([1,500])
	ax.set_ylim([1e-6,1])	

	if save_fig is not None:
		#Saves figure png and pickled data
		fig.savefig(save_fig+'.pdf',bbox_inches="tight")
		with open(save_fig + '.pkl','wb') as file:
			pickle.dump(fig, file)

def sim_plot_deltaT(filepath,deltaT,datatype,threshold=3,S_fit_max=50,bw_filter=True,timesteps=None,channels=None, save_fig=None):
	"""Plots p(S), m_av and fits p(S)~S^-alpha, for a list of binsizes. If [filepath] is a list, averages over datasets.
	
	Args:
	    filepath (str): Path to the .hdf5 dataset
	    deltaT (int): Vector of binsizes (in timesteps)
	    datatype (str): 'coarse' or 'sub'
	    threshold (int, optional): Threshold in standard deviations of the signal (for coarse) 
	    S_fit_max (int, optional): Limit on the power law fit range (default 50)
	    bw_filter (bool, optional): Whether to use a Butterworth filter to bandpass the signal (for coarse)
	    timesteps (None, optional): Number of timesteps to use (default extracts from dataset)
	    channels (None, optional): Number of electrode channels to use (default extracts from dataset)   
	    save_fig (str, optional): Saves the figure under fig/[save_fig].png
	"""

	if type(filepath) is not list:
		filepath = [filepath]

	#Parameters
	timescale_ms = 2 #length of a timestep, in ms
	fig_dir = 'fig/'

	#Lets avoid list shenenigans
	deltaT = np.array(deltaT)
	nbins = deltaT.size
	nreps = len(filepath)

	S_list = []
	deltaT_ms = np.array(timescale_ms*deltaT)
	alpha_exp = np.zeros((nreps, nbins))
	m_av = np.zeros((nreps, nbins))

	#Runs analysis for each dataset
	for j in range(nreps):

		S_list_rep = []

		#Loads and thresholds data
		data_th = avalanche.analyze_sim_raw(filepath[j],threshold, datatype, bw_filter, timesteps, channels)

		#Bins data for each deltaT and calculates observables
		for i in range(nbins):

			data_binned = avalanche.bin_data(data=data_th,binsize=deltaT[i])

			#Obtains avalanche list S
			S = avalanche.get_S(data_binned)

			#Saves S to list
			S_list_rep.append(S)

			#Calculates alpha_exp
			fit = powerlaw.Fit(S, discrete=True, estimate_discrete=False, xmin=1, xmax=S_fit_max)
			alpha_exp[j,i] = fit.alpha

			#Calculates m_av
			m_av[j,i] = fitting.m_avalanche(data_binned)

		#Appends S list from each repetition to pS
		S_list.append(S_list_rep)

	#Sets up subplots
	fig = plt.figure(constrained_layout=True)
	gs = fig.add_gridspec(2,3)
	ax_alpha = fig.add_subplot(gs[0,2])
	ax_mav = fig.add_subplot(gs[1,2])
	ax_ps = fig.add_subplot(gs[0:2,0:2])

	#Plots pS_mean for each binsize
	for k in range(nbins):

		#Rebuild S_list for each binsize
		S_bin = []
		for j in range(nreps):
			S_bin.append(S_list[j][k])

		#Gets largest avalanche from the list (+1 for zero_index)
		S_max = int(max([Si.max() for Si in S_bin]) + 1)

		#Obtains pS
		pS = np.zeros((nreps,S_max))

		for i in range(nreps):
			for j in range(S_max):
				pS[i,j] = np.sum(S_bin[i]==j)
			pS[i,:] = pS[i,:]/np.sum(pS[i,:])

		#Obtains mean and STD
		pS_mean = np.mean(pS,axis=0)
		pS_std = np.std(pS,axis=0)
		pS_up = pS_mean + pS_std/2
		pS_dw = pS_mean - pS_std/2	

		#Plots pS_mean
		str_leg = r'$\Delta$t = {:d} ms'.format(timescale_ms*deltaT[k])
		plt.fill_between(range(S_max),pS_up,pS_dw,alpha=0.25,lw=0)
		plt.plot(range(S_max),pS_mean,label=str_leg)

	#Plots alpha_exp
	alpha_exp_mean = np.mean(alpha_exp,axis=0)
	alpha_exp_std = np.std(alpha_exp,axis=0)
	ax_alpha.errorbar(deltaT_ms,alpha_exp_mean,yerr=alpha_exp_std/2,fmt='o-',color='k', fillstyle='full')

	#Plots m_av
	m_av_mean = np.mean(m_av,axis=0)
	m_av_std = np.std(m_av,axis=0)
	ax_mav.errorbar(deltaT_ms,m_av_mean,yerr=m_av_std/2,fmt='o-',color='k', fillstyle='full')
	
	#Beatifies plots
	ax_ps.set_xlabel('S')
	ax_ps.set_ylabel('p(S)')
	ax_ps.set_xlim(1,1e3)	
	ax_ps.set_yscale('log')
	ax_ps.set_xscale('log')
	ax_ps.legend()
	ax_ps.set_title(datatype+'-sampled')

	ax_alpha.set_xlabel(r'$\Delta$t (ms)')
	ax_alpha.set_ylabel(r'$\alpha$')
	ax_alpha.set_xscale('log')
	ax_alpha.set_xticks([1,10])
	ax_alpha.set_xlim([1,64])
	ax_alpha.plot([1,100],[1.5,1.5],'--')	

	ax_mav.set_xlabel(r'$\Delta$t (ms)')
	ax_mav.set_ylabel(r'$m_{av}$')
	ax_mav.set_xscale('log')
	ax_mav.set_xticks([1,10])
	ax_alpha.set_xlim([1,64])
	ax_mav.plot([1,100],[1,1],'--')	

	# #Fix the horrible ticks
	# for ax in [ax_mav,ax_alpha]:
	# 	for axis in [ax.xaxis, ax.yaxis]:
	# 		formatter = matplotlib.ticker.ScalarFormatter()
	# 		formatter.set_scientific(False)
	# 		axis.set_major_formatter(formatter)

	if save_fig is not None:
		fig_path = fig_dir+'/'.join(save_fig.split('/')[:-1])+'/'
		if not os.path.exists(fig_path):
			os.makedirs(fig_path)
		ax_ps.set_title(save_fig)
		ax_ps.set_ylabel(datatype+'-sampled p(S)')
		str_save = fig_dir+save_fig+'_'+datatype+'.png'
		fig.savefig(str_save,bbox_inches="tight")		


def spk_plot_pS(filepath,deltaT,datapath,sampling_rate = 25000,str_leg='Data'):
	"""Plots the avalanche size distribution for a single hdf5 dataset, and a single deltaT
	
	Args:
	    filepath (str): Path to the .hdf5 dataset
	    deltaT (int): Binsize to use, *in timesteps*
	    datatype (str): 'coarse', 'sub' or 'both' (compares both sampling types) 
	    str_leg (str, optional): plot legend
	    threshold (int, optional): Threshold in standard deviations of the signal (for coarse) 
	    bw_filter (bool, optional): Toggles butterworth filtering (for coarse)
	    timesteps (None, optional): Length of the dataset in timesteps
	    channels (None, optional): Number of channels in the dataset
	"""

	file = h5py.File(filepath,'r')
	data_spk = file[datapath][:]

	#Loads and thresholds data
	data_th = avalanche.analyze_sim_raw(filepath, threshold, datatype, bw_filter, timesteps, channels)

	#Bins data
	data_binned = avalanche.bin_data(data=data_th,binsize=deltaT)

	#Gets S and plots it
	S = avalanche.get_S(data_binned)
	plot.pS(S,label=str_leg)

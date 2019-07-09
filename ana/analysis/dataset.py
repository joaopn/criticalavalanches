# -*- coding: utf-8 -*-
# @Author: Joao
# @Date:   2019-07-05 17:56:44
# @Last Modified by:   Joao
# @Last Modified time: 2019-07-06 13:59:53

from analysis import avalanche, fitting, plot
import powerlaw
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Qt5Agg')

def sim_plot_pS(filepath,deltaT,datatype,str_leg='Data', threshold=3,bw_filter=True,timesteps=None,channels=None):
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
	if datatype is 'both':
		for datatype_i in ['sub','coarse']:
			#Loads and thresholds data
			data_th = avalanche.analyze_sim_raw(filepath, threshold, datatype_i, bw_filter, timesteps, channels)

			#Bins data
			data_binned = avalanche.bin_data(data=data_th,binsize=deltaT)

			#Gets S and plots it
			S = avalanche.get_S(data_binned)
			plot.pS(S,label=str_leg + ' ,' + datatype_i)			

	else:
		#Loads and thresholds data
		data_th = avalanche.analyze_sim_raw(filepath, threshold, datatype, bw_filter, timesteps, channels)

		#Bins data
		data_binned = avalanche.bin_data(data=data_th,binsize=deltaT)

		#Gets S and plots it
		S = avalanche.get_S(data_binned)
		plot.pS(S,label=str_leg)

def spk_plot_pS(filepath,deltaT,datapath,str_leg='Data', sampling_rate = 25000):
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

def sim_plot_deltaT(filepath,deltaT,datatype,threshold=3,S_fit_max=50,bw_filter=True,timesteps=None,channels=None):
	"""Plots p(S), m_av and fits p(S)~S^-alpha, for a list of binsizes.
	
	Args:
	    filepath (str): Path to the .hdf5 dataset
	    deltaT (int): Vector of binsizes (in timesteps)
	    datatype (str): 'coarse' or 'sub'
	    threshold (int, optional): Threshold in standard deviations of the signal (for coarse) 
	    S_fit_max (int, optional): Limit on the power law fit range (default 50)
	    bw_filter (bool, optional): Whether to use a Butterworth filter to bandpass the signal (for coarse)
	    timesteps (None, optional): Number of timesteps to use (default extracts from dataset)
	    channels (None, optional): Number of electrode channels to use (default extracts from dataset)
	"""
	#Parameters
	timescale_ms = 2 #length of a timestep, in ms

	#Lets avoid list shenenigans
	deltaT = np.array(deltaT)
	nbins = deltaT.size

	#Loads and thresholds data
	data_th = avalanche.analyze_sim_raw(filepath, threshold, datatype, bw_filter, timesteps, channels)

	#Observables
	pS_list = []
	alpha_exp = np.zeros(nbins)
	m_av = np.zeros(nbins)

	#Bins data for each deltaT and calculates observables
	for i in range(nbins):
		data_binned = avalanche.bin_data(data=data_th,binsize=deltaT[i])

		#Calculates p(S)
		S = avalanche.get_S(data_binned)
		S_max = int(S.max())
		pS = np.zeros(S_max)
		for j in range(S_max):
			pS[j] = np.sum(S==j+1)
		pS = pS/np.sum(pS)
		pS_list.append(pS)

		#Calculates alpha_exp
		fit = powerlaw.Fit(S, discrete=True, estimate_discrete=False, xmin=1, xmax=S_fit_max)
		alpha_exp[i] = fit.alpha

		#Calculates m_av
		m_av[i] = fitting.m_avalanche(data_binned)

	#Sets up subplots
	fig = plt.figure(constrained_layout=True)
	gs = fig.add_gridspec(2,3)
	ax_alpha = fig.add_subplot(gs[0,2])
	ax_mav = fig.add_subplot(gs[1,2])
	ax_ps = fig.add_subplot(gs[0:2,0:2])

	# #Fix the horrible ticks
	# for ax in [ax_mav,ax_alpha]:
	# 	for axis in [ax.xaxis, ax.yaxis]:
	# 		formatter = matplotlib.ticker.ScalarFormatter()
	# 		formatter.set_scientific(False)
	# 		axis.set_major_formatter(formatter)

	#Plots pS
	ax_ps.set_xlabel('S')
	ax_ps.set_ylabel('p(S)')
	plt.xlim(1,1e3)	
	for i in range(nbins):
		str_leg = r'$\Delta$t = {:d} ms'.format(timescale_ms*deltaT[i])
		ax_ps.loglog(pS_list[i],label=str_leg)
	ax_ps.legend()
	ax_ps.set_title(datatype+'-sampled')

	#Plots alpha_exp
	ax_alpha.set_xlabel(r'$\Delta$t (ms)')
	ax_alpha.set_ylabel(r'$\alpha$')
	ax_alpha.set_xscale('log')
	ax_alpha.set_xticks([1,10])
	ax_alpha.set_xlim([1,64])
	ax_alpha.plot([1,100],[1.5,1.5],'--')
	ax_alpha.plot(timescale_ms*deltaT,alpha_exp,'o-',color='k', fillstyle='full')

	#Plots m_av
	ax_mav.set_xlabel(r'$\Delta$t (ms)')
	ax_mav.set_ylabel(r'$m_{av}$')
	ax_mav.set_xscale('log')
	ax_mav.set_xticks([1,10])
	ax_alpha.set_xlim([1,64])
	ax_mav.plot([1,100],[1,1],'--')
	ax_mav.plot(timescale_ms*deltaT,m_av,'o-',color='k', fillstyle='full')
	
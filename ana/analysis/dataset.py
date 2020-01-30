# -*- coding: utf-8 -*-
# @Author: Joao
# @Date:   2019-07-05 17:56:44
# @Last Modified by:   joaopn
# @Last Modified time: 2019-12-21 22:57:09

"""
Module for directly handling datasets.
"""

from analysis import avalanche, fitting, plot, parser
from scipy.optimize import curve_fit
import powerlaw
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, os, pickle
import h5py
import pandas as pd

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


def sim_plot_thresholded(filepath, datatype, deltaT, str_leg=None, shape_d = [5,7], reps = None, save_fig=False):
	"""Analyzes and plots thresholded data in the format (reps,timesteps)
	
	Args:
		filepath (str): Path to the .h5 file
		datatype (str): Path (inside the .h5 file) to the data
		deltaT (float): binsize (in timesteps)
	"""

	#Plotting fixed parameters
	exp_avg = 2.00
	zorder = 2
	lineType = '-'
	color1 = 'r'
	color2='b'
	color3=plt.cm.get_cmap('Set1', len(shape_d)).colors
	show_error = True
	title = filepath.split('/')[-1]
	fig_path = 'results/shape/'

	if datatype not in ['coarse', 'sub']:
		ValueError('Invalid datatype. Valid options: "coarse", "sub"')

	#Loads data and bins it
	data_thresholded = h5py.File(filepath,'r')
	if reps is None:
		reps = data_thresholded[datatype].shape[0]
	fig1 = plt.figure()

	if str_leg is None:
		str_leg =  datatype

	#Plots p(S)
	print('Plotting avalanche distributions')
	S_list = []
	for i in range(reps):
		data_rep = data_thresholded[datatype][i,:]
		data_bin = avalanche.bin_data(data=data_rep, binsize=deltaT)
		S_list.append(avalanche.get_S(data_bin))
	plot.pS_mean(S_list,label='p(S), '+ str_leg,lineType=lineType,color=color1,show_error=show_error, zorder=zorder)
	
	ax = plt.gca()
	ax.set_xlim([1,500])
	ax.set_ylim([1e-6,1])
	ax.set_title(title)

	#Plots p(D)
	D_list = []
	for i in range(reps):
		data_rep = data_thresholded[datatype][i,:]
		data_bin = avalanche.bin_data(data=data_rep, binsize=deltaT)
		D_list.append(avalanche.get_D(data_bin))
	plot.pS_mean(D_list,label='p(D), '+ str_leg,lineType=lineType,color=color2,show_error=show_error, zorder=zorder)
	
	ax.set_ylabel('Probability')	

	#Plots avalanche shapes
	print('Plotting avalanche shapes')
	
	fig2 = plt.figure('shapes')
	ax2 = plt.gca()
	ax2.set_xlim([1,max(shape_d)])

	fig3 = plt.figure('shapes_collapsed')
	ax3 = plt.gca()	
	ax3.set_xlim([0,1])

	for j in range(len(shape_d)):
		shape_list = []

		for i in range(reps):
			data_rep = data_thresholded[datatype][i,:]
			data_bin = avalanche.bin_data(data=data_rep, binsize=deltaT)
			shape_list.append(avalanche.get_shape(data_bin, shape_d[j]))
		label_shape = 'D = {:d}, '.format(shape_d[j]) + str_leg

		plt.figure('shapes')
		plot.shape_mean(shape_list,size_d=shape_d[j],label=label_shape,lineType=lineType,color=color3[j,:],show_error=show_error, zorder=zorder)

		plt.figure('shapes_collapsed')
		plot.shape_mean(shape_list,size_d=shape_d[j],label=label_shape,lineType=lineType,color=color3[j,:],show_error=show_error, zorder=zorder, collapse_exp=exp_avg)

	if save_fig is not None:
		#Saves figure png and pickled data
		if not os.path.exists(fig_path):
				os.makedirs(fig_path)		
		str_save1 = fig_path + title + '_' + datatype + '_distribution_bs{:d}.png'.format(deltaT)
		str_save2 = fig_path + title + '_' + datatype + '_shape_bs{:d}.png'.format(deltaT)
		str_save3 = fig_path + title + '_' + datatype + '_shape_collapsed_bs{:d}.png'.format(deltaT)
		fig1.savefig(str_save1,bbox_inches="tight")
		fig2.savefig(str_save2,bbox_inches="tight")
		fig3.savefig(str_save3,bbox_inches="tight")
		plt.close(fig1)
		plt.close(fig2)
		plt.close(fig3)

	
def sim_plot_scaling(filepath, deltaT, reps = None, xmax_S = None, xmax_D = None, xmax_avg = None, save_fig=True):

	#Plotting fixed parameters
	zorder = 2
	lineType = '-'
	color_dist = {'sub':'b', 'coarse':'r'}
	show_error = True
	title = filepath.split('/')[-1]
	fig_path = 'results/shape/'
	collapse_min_d = 4
	collapse_min_rep = 20

	#xmax
	if xmax_S is None:
		xmax_S = {'sub': 300, 'coarse':64}
	if xmax_D is None:
		xmax_D = {'sub': 100, 'coarse':8}
	if xmax_avg is None:
		xmax_avg = {'sub': 100, 'coarse':8}

	#deltaT
	if len(deltaT) == 1:
		deltaT = deltaT * 2
	elif len(deltaT) > 2:
		raise ValueError('deltaT needs one or two values')
	deltaT_dict = {'coarse':deltaT[0], 'sub':deltaT[1]}

	#Results dict
	observables = ['alpha','beta','gamma','gamma_shape']
	results = {}

	#Loads data and bins it
	data_thresholded = h5py.File(filepath,'r')
	if reps is None:
		reps = data_thresholded['coarse'].shape[0]
	#fig1 = plt.figure()

	#Creates figures
	# fig1 = plt.figure('pS')
	# fig2 = plt.figure('pD')
	# fig3 = plt.figure('average_S')
	# fig4 = plt.figure('shape_coarse')
	# fig5 = plt.figure('shape_sub')
	# fig6 = plt.figure('bar_comparison')
	fig = plt.figure(constrained_layout=True, figsize=(8,6))
	gs = fig.add_gridspec(2,3)
	ax_pS = fig.add_subplot(gs[0,0])
	ax_pD = fig.add_subplot(gs[0,1])
	ax_avgS = fig.add_subplot(gs[0,2])
	ax_collapse_coarse = fig.add_subplot(gs[1,0])
	ax_collapse_sub = fig.add_subplot(gs[1,1])
	ax_comparison = fig.add_subplot(gs[1,2])

	for datatype in ['coarse', 'sub']:

		results[datatype] = dict.fromkeys(observables)
		all_avalanches = {}

		#Calculates avalanche statistics for all reps
		for i in range(reps):
			data_rep = data_thresholded[datatype][i,:]
			data_bin = avalanche.bin_data(data=data_rep, binsize=deltaT_dict[datatype])		
			all_avalanches[i] = avalanche.get_avalanches(data_bin)

		#Obtains lists of S, D, and shape
		S_list = [[all_avalanches[i][j]['S'] for j in range(len(all_avalanches[i]))] for i in range(reps)] 
		D_list = [[all_avalanches[i][j]['D'] for j in range(len(all_avalanches[i]))] for i in range(reps)] 
		shape_list = [[all_avalanches[i][j]['shape'] for j in range(len(all_avalanches[i]))] for i in range(reps)] 

		#return shape_list

		#avgS_list = [[[all_avalanches[i][j]['S'] for j in range(len(all_avalanches[i]))] ] for i in range(reps)]

		#Calculates average avalanche size
		avgS_list = []		#List, averaged by rep
		avgS = {}

		D_max = int(max([max(Di) for Di in D_list]))
		for r in range(reps):
			S_rep = [[all_avalanches[r][V]['S'] for V in all_avalanches[r].keys() if all_avalanches[r][V]['D'] == i] for i in range(1,D_max+1)]
			avgS_list.append([np.mean(V) for V in S_rep])

			if r == 0:
				avgS['D'] = [all_avalanches[r][V]['D'] for V in sorted(all_avalanches[r].keys())]
				avgS['S'] = [all_avalanches[r][V]['S'] for V in sorted(all_avalanches[r].keys())]
			else:
				avgS['D'] = avgS['D'] + [all_avalanches[r][V]['D'] for V in sorted(all_avalanches[r].keys())]
				avgS['S'] = avgS['S'] + [all_avalanches[r][V]['S'] for V in sorted(all_avalanches[r].keys())]

		avgS['X'] = [*range(1,max(avgS['D'])+1)]
		maxD = max(avgS['D'])

		#Converts to a non-bullshit format
		avgS['D'] = np.array(avgS['D'])
		avgS['S'] = np.array(avgS['S'])
		avgS['Savg'] = np.zeros(maxD)
		avgS['Sstd'] = np.zeros(maxD)

		for i in range(maxD):
			S_select = avgS['S'][avgS['D']==i+1]
			#S_select = avgS['S'][np.where(np.equal(avgS['D'],i+1))[0]]
			if S_select.size > 0:
				avgS['Savg'][i] = np.mean(S_select)
				avgS['Sstd'][i] = np.std(S_select)
		
		avgS_X = np.array(avgS['X'][:])
		avgS_Y = avgS['Savg'][:]
		avgS_Yerr = avgS['Sstd'][:]

		#Flattens dict into lists
		S = [item for sublist in S_list for item in sublist]
		D = [item for sublist in D_list for item in sublist]
		shape_list_flat = [item for sublist in shape_list for item in sublist]

		#print(avgS['Savg'])

		#S = [all_avalanches[i]['S'] for i in range(len(all_avalanches))]
		#D = [all_avalanches[i]['D'] for i in range(len(all_avalanches))]
		#shape_av = [all_avalanches[i]['shape'] for i in range(len(all_avalanches))]

		#Gets p(S)
		# S_list = []
		# for i in range(reps):
		# 	data_rep = data_thresholded[datatype][i,:]
		# 	data_bin = avalanche.bin_data(data=data_rep, binsize=deltaT_dict[datatype])
		# 	S_list.append(avalanche.get_S(data_bin))

		#Fits p(S)
		fit = powerlaw.Fit(S, discrete=True, estimate_discrete=False, xmin=1, xmax=xmax_S[datatype])
		results[datatype]['alpha'] = fit.alpha	

		#Plots p(S)
		plt.sca(ax_pS)
		str_leg = datatype + r': $\alpha$ = {:0.2f}'.format(results[datatype]['alpha'])
		plot.pS_mean(S_list,label=str_leg,lineType=lineType,color=color_dist[datatype],show_error=show_error, zorder=zorder)
		X_fit = np.array(range(1,xmax_S[datatype]+1))
		Y_fit = np.power(X_fit, -results[datatype]['alpha'])
		plt.plot(X_fit, Y_fit/np.sum(Y_fit), '--', color=color_dist[datatype])

		# plot.pS_mean(S_list,label=str_leg,lineType=lineType,color=color_dist[datatype],show_error=show_error, zorder=zorder)

		#Gets p(D)
		# D_list = []
		# for i in range(reps):
		# 	data_rep = data_thresholded[datatype][i,:]
		# 	data_bin = avalanche.bin_data(data=data_rep, binsize=deltaT_dict[datatype])
		# 	D_list.append(avalanche.get_D(data_bin))

		#Fits p(D)
		# D = [item for sublist in D_list for item in sublist]
		fit = powerlaw.Fit(D, discrete=True, estimate_discrete=False, xmin=1, xmax=xmax_D[datatype])
		results[datatype]['beta'] = fit.alpha	

		#Plots p(D)
		plt.sca(ax_pD)
		str_leg = datatype + r': $\beta$ = {:0.2f}'.format(results[datatype]['beta'])
		plot.pS_mean(D_list,label=str_leg,lineType=lineType,color=color_dist[datatype],show_error=show_error, zorder=zorder)
		X_fit = np.array(range(1,xmax_D[datatype]+1))
		Y_fit = np.power(X_fit, -results[datatype]['beta'])
		plt.plot(X_fit, Y_fit/np.sum(Y_fit), '--', color=color_dist[datatype])

		# plot.pS_mean(D_list,label=str_leg,lineType=lineType,color=color_dist[datatype],show_error=show_error, zorder=zorder)

		#Gets estimation of gamma from pS and pD
		results[datatype]['gamma_p'] = (results[datatype]['beta'] - 1)/(results[datatype]['alpha'] - 1)

		#Censors data to fit
		bool_keep = np.less(avgS_X,xmax_avg[datatype])
		avgS_X = avgS_X[bool_keep]
		avgS_Y = avgS_Y[bool_keep]
		avgS_Yerr = avgS_Yerr[bool_keep]

		#Fits avgS
		#kwargs = {'maxfev': 10000}
		# def __pl(x,a,b):
		# 	return a*np.power(x,b)
		# results_fit, pcov = curve_fit(__pl,avgS_X,avgS_Y, method='lm', **kwargs)
		#results[datatype]['gamma'] = results_fit[1]
		results[datatype]['gamma'], fit_err, lin_coef = fitting.powerlaw(avgS_X, avgS_Y, avgS_Yerr, loglog=True)

		#Plots avgS
		plt.sca(ax_avgS)
		str_leg = datatype + r': $\gamma$ = {:0.2f}'.format(results[datatype]['gamma'])
		plot.avgS_mean(avgS_list,maxD,label=str_leg,lineType=lineType,color=color_dist[datatype],show_error=show_error, zorder=zorder)

		#Plots avgS fit
		X_fit = np.array(range(1,xmax_avg[datatype]+1))
		#plt.plot(X_fit, results_fit[0]*np.power(X_fit, results_fit[1]), '--', color=color_dist[datatype])
		plt.plot(X_fit, lin_coef*np.power(X_fit, results[datatype]['gamma']), '--', color=color_dist[datatype])

		#Fits shape scaling
		results[datatype]['gamma_shape'] = fitting.shape_collapse(shape_list, collapse_min_d, collapse_min_rep)

		#Plots shape scaling
		if datatype == 'coarse':
			ax_scaled = ax_collapse_coarse
		elif datatype == 'sub':
			ax_scaled = ax_collapse_sub
		plot.shape_collapse(shape_list,results[datatype]['gamma_shape'], collapse_min_d, collapse_min_rep, ax=ax_scaled)

	#Beautifies plots
	plt.sca(ax_pS)
	ax = plt.gca()
	ax.set_xlim([1,500])
	ax.set_ylim([1e-6,1])
	#ax.set_title(title)
	ax.set_xlabel('S')
	ax.set_ylabel('p(S)')

	plt.sca(ax_pD)
	ax = plt.gca()
	ax.set_xlim([1,100])
	ax.set_ylim([1e-6,1])
	#ax.set_title(title)
	ax.set_xlabel('D')
	ax.set_ylabel('p(D)')

	print(results)

	#Plots bar plot
	plt.sca(ax_comparison)
	ax = plt.gca()
	#ax.set_title(title)
	ax.set_ylim([0,3.5])
	str_gamma_P = r'($\beta-1$)/($\alpha-1$)'
	str_gamma_avg = r'$\gamma$'
	str_gamma_shape = r'$\gamma_{s}$'
	results_gamma = {'coarse-sampled':{str_gamma_P:results['coarse']['gamma_p'],str_gamma_avg:results['coarse']['gamma'], str_gamma_shape:results['coarse']['gamma_shape']},
	 'sub-sampled':{str_gamma_P:results['sub']['gamma_p'],str_gamma_avg:results['sub']['gamma'], str_gamma_shape:results['sub']['gamma_shape']}}
	pd.DataFrame(results_gamma).T.plot(kind='bar', ax=ax)
	plt.xticks(rotation='horizontal')



	if save_fig is True:
		#Saves figure png and pickled data
		if not os.path.exists(fig_path):
				os.makedirs(fig_path)	
		str_save = 	fig_path + title + '_' + 'bs{:d}-{:d}.png'.format(deltaT[0],deltaT[1])
		# str_save1 = fig_path + title + '_' + 'bs{:d}-{:d}_pS.png'.format(deltaT[0],deltaT[1])
		# str_save2 = fig_path + title + '_' + 'bs{:d}-{:d}_pD.png'.format(deltaT[0],deltaT[1])
		# str_save3 = fig_path + title + '_' + 'bs{:d}-{:d}_avgS.png'.format(deltaT[0],deltaT[1])
		# str_save6 = fig_path + title + '_' + 'bs{:d}-{:d}_comparison.png'.format(deltaT[0],deltaT[1])
		# fig1.savefig(str_save1,bbox_inches="tight")
		# fig2.savefig(str_save2,bbox_inches="tight")
		# fig3.savefig(str_save3,bbox_inches="tight")
		# fig6.savefig(str_save6,bbox_inches="tight")
		fig.suptitle(title)
		fig.savefig(str_save,bbox_inches="tight")
		plt.close(fig)
		# plt.close(fig1)
		# plt.close(fig2)
		# plt.close(fig3)
		# plt.close(fig4)
		# plt.close(fig5)
		# plt.close(fig6)
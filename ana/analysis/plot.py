# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-26 13:40:21
# @Last Modified by:   Joao
# @Last Modified time: 2019-12-12 02:36:35

"""

Module for plotting distributions. 

"""
import analysis.avalanche
import analysis.fitting
import matplotlib, os
import matplotlib.pyplot as plt
import numpy as np


def pS(S,label='data'):

	#Calculates p(S) (no log-binning)
	S_max = int(S.max()+1)
	pS = np.zeros(S_max)

	for i in range(1,S_max):
		pS[i] = np.sum(S==i)
	pS = pS/np.sum(pS)

	#Plots it
	plt.loglog(pS,label=label)
	plt.xlim(1,1e3)
	plt.legend()

def pS_mean(S_list,label='data',lineType='-', color='k',show_error=True, zorder=2, add_zero = False):

	if add_zero is True:
		extra_av = 2
	else:
		extra_av = 1
	#Gets largest avalanche from the list (+extra_av for zero_index)
	S_max = int(max([max(Si) for Si in S_list]) + extra_av)
	#S_max = np.amax(S_list)

	#Obtains pS
	pS = np.zeros((len(S_list),S_max))

	for i in range(len(S_list)):
		for j in range(S_max):
			pS[i,j] = np.sum(np.equal(S_list[i],j))
		pS[i,:] = pS[i,:]/np.sum(pS[i,:])
		
	#Obtains mean and STD
	pS_mean = np.mean(pS,axis=0)
	pS_std = np.std(pS,axis=0)
	pS_up = pS_mean + pS_std/2
	pS_dw = pS_mean - pS_std/2

	#Plots confidence interval (1 std) and mean
	if show_error:
		plt.fill_between(range(S_max),pS_up,pS_dw,alpha=0.25,lw=0,color=color,zorder=zorder)
	plt.plot(range(S_max),pS_mean,lineType,label=label,color=color,zorder=zorder)
	plt.yscale('log')
	plt.xscale('log')
	plt.legend(loc='upper right', frameon=False)

def avgS_mean(avgS_list,Dmax,label='data',lineType='-', color='k',show_error=True, zorder=2):
		
	#Obtains mean and STD
	avgS_mean = np.mean(avgS_list,axis=0)
	avgS_std = np.std(avgS_list,axis=0)
	avgS_up = avgS_mean + avgS_std/2
	avgS_down = avgS_mean - avgS_std/2

	#Plots confidence interval (1 std) and mean
	if show_error:
		plt.fill_between(range(1,Dmax+1),avgS_up,avgS_down,alpha=0.25,lw=0,color=color,zorder=zorder)
	plt.plot(range(1,Dmax+1),avgS_mean,lineType,label=label,color=color,zorder=zorder)
	plt.yscale('log')
	plt.xscale('log')
	plt.legend(frameon=False)	

def shape_mean(shape_list, size_d, label='data',lineType='-', color='k',show_error=True, zorder=2, collapse_exp=None):

	#Obtains mean and STD
	shape_mean = np.mean(shape_list,axis=0)
	shape_std = np.std(shape_list,axis=0)
	shape_up = shape_mean + shape_std/2
	shape_down = shape_mean - shape_std/2

	if collapse_exp is None:
		x_plot = np.array(range(1,size_d+1))
	else:
		shape_mean = shape_mean/size_d**(collapse_exp-1)
		shape_up = shape_up/size_d**(collapse_exp-1)
		shape_down = shape_down/size_d**(collapse_exp-1)
		x_plot = np.array(range(1,size_d+1))/size_d

	#Plots confidence interval (1 std) and mean
	if show_error:
		plt.fill_between(x_plot,shape_up,shape_down,alpha=0.25,lw=0,color=color,zorder=zorder)
	plt.plot(x_plot,shape_mean,lineType,label=label,color=color,zorder=zorder)
	plt.legend()

def timeseries_threshold(data,th):
	data_th = analysis.avalanche.threshold_ch(data,th)

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

		#Defines error band
		pS_up = pS_mean + pS_std/2
		pS_dw = pS_mean - pS_std/2

		#Fixes vector bug for S_mean
		S_rem = np.argwhere(pS_mean < S_lim)
		pS_mean = np.delete(pS_mean,S_rem)
		S_mean = np.delete(S,S_rem)
		S_rem = np.argwhere(pS_dw < S_lim)
		pS_dw = np.delete(pS_dw,S_rem)
		pS_up = np.delete(pS_up,S_rem)
		S_std = np.delete(S,S_rem)

		#Plots data
		if plt_std:
			plt.fill_between(S_std,offset*pS_up,offset*pS_dw,alpha=0.5,color=plt_color, lw=0)

		if label_plot is None:
			tau_rep = analysis.fitting.tau_sim_dataset(m,h,d[i],threshold,data_dir,bw_filter)
			str_label = r'$\tau$ = {:0.1f} $\pm$ {:0.1f} ms'.format(np.mean(tau_rep), np.std(tau_rep))
			plt.plot(S_mean,offset*pS_mean,label= str_label,color=plt_color)
		else:			
			plt.plot(S_mean,offset*pS_mean,label=label_plot[i],
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
		plt.fill_between(IED,mav_up,mav_dw,alpha=0.5,color=plt_color,lw=0)

		if label_plot is None:

			print('Calculating tau for b = {:02d}'.format(b))
			tau_mean, tau_std = analysis.fitting.tau_sim_dataset(m,h,IED,threshold,data_dir,bw_filter)

			print('tau = {:0.1f} +- {:0.1f}'.format(tau_mean, tau_std))
			str_label = r'$\Delta$t = {:d} ms, $\bar\tau$ = {:0.1f} ms'.format(2*b,tau_mean)

			plt.plot(IED,mav_mean,label= str_label,color=plt_color,linestyle=linestyle)
		else:
			plt.plot(IED,mav_mean,label=label_plot,color=plt_color,linestyle=linestyle)

def plot_alpha_bs(m,h,b,d,datatype,reps,bw_filter = False,data_dir ='dat/',threshold = 3, color_rgb=None):

	#Definitions
	deltaT = 2
	if bw_filter:
		saveplot_dir = 'analyzed_filtered/'
	else:
		saveplot_dir = 'analyzed_unfiltered/'

	#Sets up filepath and loads data
	dataset = 'm{:0.5f}_h{:0.3e}_d{:02d}_th{:0.1f}_rep{:02d}/'.format(m,h,d,threshold,reps)
	datafile = 'alpha_' + datatype + '.tsv'
	str_load = data_dir + saveplot_dir + dataset + datafile	
	b_data, alpha_mean_data, alpha_std_data = np.loadtxt(str_load,delimiter='\t')

	#Censors data, only using data for the given b
	b_deltat = deltaT*np.array(b,dtype=float)
	b_sort = np.argsort(b_data)
	b_index = b_sort[np.searchsorted(b_data, b_deltat, sorter=b_sort)]
	b = b_data[b_index]
	alpha_mean = alpha_mean_data[b_index]
	alpha_std = alpha_std_data[b_index]

	#Fits data to alpha~b^-fit_exp
	fit_exp, fit_err, lin_coef = analysis.fitting.powerlaw(b,alpha_mean,alpha_std)

	#Plots fit
	X = np.arange(0.1,100)
	plt.plot(X,lin_coef*np.power(X,-fit_exp),linestyle='--', color='k', zorder=1)

	#Plots data
	if color_rgb is not None:
		plt.errorbar(b, alpha_mean, yerr=alpha_std/2, ecolor='k',fmt='none', elinewidth=1)
		plt.scatter(b,alpha_mean,marker='s',color=color_rgb, s=16, zorder=2)
		
	else:
		plt.errorbar(b, alpha_mean, yerr=alpha_std/2, ecolor='k', fmt='s', s=16)

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

def sim_corr(m,h,b,d,threshold,reps,data_dir,type='corr', loc=1, color_state=[0,0,0]):

	#Parameters
	dir_corr = 'correlations/'
	bar_width = 0.125
	color_coarse = '#800080' #purple
	color_sub = '#5F9EA0' #cadetblue

	#Sets up filepath and loads data
	datafile = 'm{:0.5f}_h{:0.3e}__d{:02d}_b{:02d}_th{:0.1f}_rep{:02d}.tsv'.format(m,h,b,d,threshold,reps)
	str_load = data_dir + dir_corr + datafile
	corr_coarse,corr_sub, rate_coarse, rate_sub = np.loadtxt(str_load,delimiter='\t')

	#Calculates statistics
	if type=='corr':
		coarse_mean = np.mean(corr_coarse)
		coarse_std = np.std(corr_coarse)
		sub_mean = np.mean(corr_sub)
		sub_std = np.std(corr_sub)

	elif type == 'rate':		
		coarse_mean = np.mean(rate_coarse)
		coarse_std = np.std(rate_coarse)
		sub_mean = np.mean(rate_sub)
		sub_std = np.std(rate_sub)

	#Plots barplot at loc
	Y = [coarse_mean, sub_mean]
	Y_err = [coarse_std, sub_std]
	bar_loc = [loc-bar_width,loc+bar_width]

	plt.bar(bar_loc[0],Y[0],yerr=Y_err[0], capsize=3, width=2*bar_width, color=[1,1,1], hatch='////', edgecolor=color_state, linewidth=1)
	plt.bar(bar_loc[1],Y[1],yerr=Y_err[1], capsize=3, width=2*bar_width, color=color_state, linewidth=0.5)

def _color_picker(color_str):

	if color_str in ['red', 'r']:
		color_rgb = [255,51,51]
	elif color_str in ['blue', 'b']:
		color_rgb = [51,102,153]
	elif color_str in ['green', 'g']:
		color_rgb = [0,153,51]
	elif color_str in ['gray', 'grey']:
		color_rgb = [128,128,128]

	return np.array(color_rgb)/255

def _color_gradient_rgba(color,num_colors, brightness_lim=0.8, primary_index = 0.545):
	"""Creates gradient colors RGBA matrices for basic colors
	
	Args:
		color: color string, accepts both e.g. 'r' and 'red'
		num_colors: number of colors to generate
		brightness_lim: maximum brightness (default: {0.8})
		primary_index: index of the base color. Setting it to 1 generates a pure tone (default: {0.545})
	
	Returns:
		color RGBA matrix
		float: [num_colors,4]
	"""

	color_rgba = np.zeros((num_colors,4))
	color_rgba[:,3] = np.linspace((1-brightness_lim),1,num_colors)

	if color in ['r', 'red']:
		color_rgba[:,0] = primary_index

	if color in ['g', 'green']:
		color_rgba[:,1] = primary_index

	if color in ['b', 'blue']:
		color_rgba[:,2] = primary_index

	if color in ['gray', 'grey']:
		color_rgba[:,0] = 0.2
		color_rgba[:,1] = 0.2
		color_rgba[:,2] = 0.2 

	return color_rgba

def _convert_rgba_rgb(colors):

	#Deals with python dimentionality silliness
	colors = np.array(colors)
	if colors.ndim == 1:
		colors = np.expand_dims(colors,0)

	rgb_background = np.array([1,1,1], dtype = float)

	n_colors = colors.shape[0]
	colors_rgb = np.zeros((n_colors,3))

	for i in range(n_colors):
		alpha = colors[i,3]
		colors_rgb[i,:] = (1-alpha)*rgb_background + alpha*colors[i,:3]

	return colors_rgb

def shape_collapse(shape_list, gamma_shape, min_d, min_rep, ax=None, str_leg = None, extrapolate=False, color='r'):

	from scipy.interpolate import InterpolatedUnivariateSpline	

	#Definitions
	interp_points = 1000

	if ax is None:
		plt.figure()
	else:
		plt.sca(ax)

	#Flattens list of reps
	flat_list = np.array([item for sublist in shape_list for item in sublist])

	#List of avalanche sizes
	shape_size = np.zeros(len(flat_list))
	for i in range(len(flat_list)):
		shape_size[i] = flat_list[i].size

	max_size = shape_size.max()

	#Avalanche size count
	shape_count,_ = np.histogram(shape_size,bins=np.arange(0,max_size+2))

	#Censors data by size
	censor_d_keep = np.arange(0,max_size+1) >= min_d
	censor_rep_keep = shape_count >= min_rep
	censor_index =  np.where([a and b for a, b in zip(censor_d_keep, censor_rep_keep)])[0]

	#Defines average size matrix
	average_shape = np.zeros((censor_index.size, interp_points))

	#Defines bottom interpolation range from data, to prevent extrapolation bias
	#x_min = 1/censor_index[0]
	if extrapolate is True:
		x_min = 0
	elif extrapolate is False:
		x_min = 1/censor_index[0]
	else:
		error('extrapolate is not binary.')
	x_range = np.linspace(x_min,1,num=interp_points)

	#Averages shape for each duration and interpolates results
	y_min = 100
	for i in range(len(censor_index)):

		#Calculates average shape
		size_i = censor_index[i]
		avg_shape_i_y = np.mean(flat_list[shape_size==size_i])/np.power(size_i,gamma_shape-1)
		avg_shape_i_x = np.arange(1,size_i+1)/size_i

		if np.min(avg_shape_i_y) < y_min:
			y_min = np.min(avg_shape_i_y)

		#Interpolates results
		fx = InterpolatedUnivariateSpline(avg_shape_i_x,avg_shape_i_y)
		average_shape[i,:] = fx(x_range)

		#Plots transparent subplots
		plt.plot(avg_shape_i_x, avg_shape_i_y, alpha=0.2, color=color)

	#Plots interpolated average curve
	plt.plot(x_range, np.mean(average_shape, axis=0), color='k', linewidth=2)

	#Beautifies plot
	ax.set_xlabel('Scaled time')
	ax.set_ylabel('Scaled activity')
	plt.xlim([0,1])
	if extrapolate is True:
		_,ylim_1 = ax.get_ylim()
		ax.set_ylim([y_min, ylim_1])
	if str_leg is not None:
		#plt.text(0.95, 0.95, str_leg, horizontalalignment='right',verticalalignment='top', transform=ax.transAxes, bbox=dict(facecolor='none', alpha=0.5))
		plt.text(0.95, 0.95, str_leg, horizontalalignment='right',verticalalignment='top', transform=ax.transAxes)
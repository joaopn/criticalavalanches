# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-31 18:46:04
# @Last Modified by:   Joao
# @Last Modified time: 2019-07-27 14:26:00

import analysis
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py
import os
import argparse

from matplotlib import rcParams, rc
matplotlib.rcParams['pdf.fonttype'] = 42
rcParams['font.sans-serif'] = ['Arial']
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10
rc('legend', fontsize=8)

def states_tau():

	#Gets tau for legends, from reps=50 datasets encompassing d \in [2,10]
	state_tau = {
		'poisson': {'tau': 0.1342898650609103, 'std': 0.13768478203889073},
		'subcritical': {'tau': 18.964152286594654, 'std': 0.08660762654290584},
		'reverberant': {'tau': 98.26704463714672, 'std': 1.0030388400445385},
		'critical': {'tau': 1581.1362196543237, 'std': 114.9900207541611}
	}	

	return state_tau

def states_parameters():

	#Defines state dict
	state_dict = {
		'poisson': {'m': 0.000, 'h': 2e-3},
		'subcritical': {'m': 0.900, 'h': 2e-4},
		'reverberant': {'m': 0.980, 'h': 4e-5},
		'critical': {'m': 0.999, 'h': 2e-6},
	}

	return state_dict

def color_picker(color_str):

	if color_str in ['red', 'r']:
		color_rgb = [255,51,51]
	elif color_str in ['blue', 'b']:
		color_rgb = [51,102,153]
	elif color_str in ['green', 'g']:
		color_rgb = [0,153,51]
	elif color_str in ['gray', 'grey']:
		color_rgb = [128,128,128]

	return np.array(color_rgb)/255

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)	

def pS_filepath(state,datatype,d,b,reps,threshold,data_dir, bw_filter):

	#Definitions
	state_dict = states_parameters()
	if bw_filter:
		saveplot_dir = 'analyzed_filtered/'
	else:
		saveplot_dir = 'analyzed_unfiltered/'

	#Sets up filepath
	m = state_dict[state]['m']
	h = state_dict[state]['h']
	dataset = 'm{:0.5f}_h{:0.3e}_d{:02d}_th{:0.1f}_rep{:02d}/'.format(m,h,d,threshold,reps)
	datafile = 'pS_' + datatype + '_b{:02d}.tsv'.format(b)
	str_load = data_dir + saveplot_dir + dataset + datafile

	return str_load

def thresholded_filepath(state,d,threshold,data_dir,bw_filter):

	#Parameters
	state_dict = states_parameters()
	if bw_filter:
		dir_threshold = 'thresholded_filtered/'
	else:
		dir_threshold = 'thresholded_unfiltered/'

	#Sets up filepath
	m = state_dict[state]['m']
	h = state_dict[state]['h']
	filename = 'm{:0.5f}_h{:0.3e}_d{:02d}_th{:0.1f}.hdf5'.format(m,h,d,threshold)
	str_load = data_dir + dir_threshold + filename

	return str_load

def plot_compare_states(datatype,b,d,reps,threshold,data_dir,bw_filter):

	#Plots guiding line
	X = np.arange(5,64)
	Y = np.power(X,-1.5)
	Y = 5*Y/np.sum(Y)
	#plt.plot(X,Y,linestyle='--',color='black',label=r'$p(S) \sim S^{-1.5}$')
	plt.plot(X,Y,linestyle='--',color='black',label=None)

	#Plots states
	states = ['poisson','subcritical','reverberant', 'critical']
	colors = ['gray','blue', 'green', 'red']
	state_dict = states_parameters()

	for i in range(len(states)):

		#Get parameters		
		m = state_dict[states[i]]['m']
		h = state_dict[states[i]]['h']

		#Plots distributions
		color_data = color_picker(colors[i])
		analysis.plot.sim_pS(m,h,d,b,datatype,reps,states[i],bw_filter,data_dir,threshold,color_data)

def parametersDefault():

	#default Parameters
	figDefault = '1'
	binsizeDefault= 4
	repsDefault = 50
	dDefault = 8
	datafolderDefault = 'dat/'
	bw_filterDefault = True
	datatypeDefault = 'coarse'
	thresholdDefault = 3

	#Parse input
	parser = argparse.ArgumentParser()

	#Add single parameters
	parser.add_argument("--datafolder",type=str,nargs='?',const=1,default=datafolderDefault)
	parser.add_argument("--reps",type=int,nargs='?',const=1,default=repsDefault)
	parser.add_argument("-b",type=int,nargs='?',const=1,default=binsizeDefault)
	parser.add_argument("-d",type=int,nargs='?',const=1,default=dDefault)
	parser.add_argument("--threshold",type=float,nargs='?',const=1,default=thresholdDefault)
	parser.add_argument("--bw_filter",type=bool,nargs='?',const=1,default=bw_filterDefault)
	parser.add_argument("--fig",type=str,nargs='?',const=1,default=figDefault)
	parser.add_argument("--datatype",type=str,nargs='?',const=1,default=datatypeDefault)
	
	args = parser.parse_args()

	return args

def color_gradient(color,num_colors, brightness_lim=0.8):
	"""Creates gradient color RGB matrices for the basic colors.

	Supported colors: 'red', 'green', 'blue', 'grey'.
	
	Args:
		color: color string, e.g. "red" or "r"
		num_colors: number of colors to generate
		brightness_lim: maximum brightness. Setting it to 1 means the last color is white (default: {0.8})
	"""
	color_rgb = np.zeros((num_colors,3))
	v_num = np.linspace(brightness_lim,0,num=num_colors)

	if color in ['r', 'red']:
		color_rgb[:,0] = darkness
		color_rgb[:,1] = v_num
		color_rgb[:,2] = v_num

	if color in ['g', 'green']:
		color_rgb[:,0] = v_num
		color_rgb[:,1] = 1
		color_rgb[:,2] = v_num

	if color in ['b', 'blue']:
		color_rgb[:,0] = v_num
		color_rgb[:,1] = v_num
		color_rgb[:,2] = 1

	if color in ['gray', 'grey']:
		color_rgb[:,0] = 0.499 + v_num/2
		color_rgb[:,1] = 0.499 + v_num/2
		color_rgb[:,2] = 0.499 + v_num/2

	return color_rgb

def color_gradient_rgba(color,num_colors, brightness_lim=0.8, primary_index = 0.545):
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

def convert_rgba_rgb(colors):

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

def figure_mav(data_dir,b,bw_filter):

	#Parameters
	fig_size = [7,4]
	Xticks = [2, 4, 6, 8, 10]
	Xticklabels = np.array(Xticks)*50
	linestyles = ['-','--', ':', '-.'] #Handles up to len(b) = 4

	#Parse input
	if type(b) is not list:
		b = [b]

	#Sets path
	if bw_filter:
		str_savepath = 'figs/fig2_filtered/'
	else:
		str_savepath = 'figs/fig2_unfiltered/'
	if not os.path.exists(str_savepath):
		os.makedirs(str_savepath)

	#Sets up figure
	plt.figure(figsize=(fig_size[0]/2.54,fig_size[1]/2.54))
	plt.xlabel(r'$d_E$ ($\mu$m)')
	plt.ylabel(r'$m_{av}$')
	plt.xlim(2,10)
	plt.ylim(0.6,1.8)
	plt.yticks([0.6,1.,1.4, 1.8])
	plt.xticks(Xticks)	
	plt.grid(False)
	ax = plt.gca()
	ax.set_axisbelow(True)
	ax.set_xticklabels(Xticklabels)

	#Plots guiding line
	#plt.plot([1,10],[1,1],'--', color='k')
	#plt.grid(True)

	#Plots states
	states = ['subcritical','reverberant','critical']
	colors = ['blue','green','red']


	for j in range(len(b)):
		for i in range(len(states)):

			print('Plotting Fig. 2 for the ' + states[i] + ' state.')
			#Get parameters
			state_dict = states_parameters()
			m = state_dict[states[i]]['m']
			h = state_dict[states[i]]['h']

			color_state = color_picker(colors[i])

			str_label = r'$\Delta$t = {:d} ms, '.format(2*b[j]) + states[i]
			analysis.plot.sim_mav(m,h,b[j],data_dir,
				bw_filter=bw_filter, plt_color=color_state,
				linestyle=linestyles[j], label_plot=str_label)

		# if j ==1:
		# 	plt.legend(frameon=True)

	#Saves figure
	set_size(fig_size[0]/2.54,fig_size[1]/2.54)
	if len(b) == 1:
		str_save = str_savepath + 'mav_b{:02d}.pdf'.format(b)
	else:
		str_save = str_savepath + 'mav_combined.pdf'
	plt.savefig(str_save,bbox_inches="tight")
	plt.close()

def figure_1(data_dir,b,d,reps,bw_filter):

	#Parameters
	N = 160000
	threshold = 3
	rep_activity = 0 #timeseries to plot the activity from
	fig_activity_size = [5,3] #in cm
	fig_pS_size = [6,6]
	fig_pS_ylim = 1e-5

	#Sets path
	if bw_filter:
		str_savepath = 'figs/fig1_filtered/'
	else:
		str_savepath = 'figs/fig1_unfiltered/'
	if not os.path.exists(str_savepath):
		os.makedirs(str_savepath)

	"""Generates timeseries figures"""
	print('Fig1: generating timeseries figures')
	states = ['poisson','subcritical', 'reverberant', 'critical']
	colors = ['gray','blue', 'green','red']
	states_tau_dict = states_tau()

	for i in range(len(states)):

		#Sets up figure
		plt.figure(figsize=(fig_activity_size[0]/2.54,fig_activity_size[1]/2.54))
		plt.xlim((0,10000))
		plt.ylim((0,4))
		plt.ylabel('Hz')

		#Loads data
		str_thresholded = thresholded_filepath(states[i],d,threshold,data_dir,bw_filter)
		file = h5py.File(str_thresholded,'r')
		activity_data = file['activity'][rep_activity,:]
		file.close()

		#Plots data
		color_data = color_picker(colors[i])
		str_leg = r'$\bar\tau$ = {:0.1f} ms'.format(states_tau_dict[states[i]]['tau'])
		plt.plot(activity_data*1000/N/2,color=color_data, linewidth=1,label=str_leg)
		plt.legend(loc=1, frameon=False)

		#Removes x labels
		fig_gca = plt.gca()
		fig_gca.axes.xaxis.set_ticklabels([])
		fig_gca.set_yticks([0,2,4])

		#Saves plot
		str_save = str_savepath + 'activity_d{:02d}_rep{:02d}_'.format(d,rep_activity) + states[i] +'.pdf'
		plt.savefig(str_save,bbox_inches="tight")
		plt.close()


	"""Generates p(S) figures"""

	#Generates subsampled figure, removes legend
	print('Fig1: generating subsampled figure')
	plt.figure(figsize=(fig_pS_size[0]/2.54,fig_pS_size[1]/2.54))
	plot_compare_states('sub',b,d,reps,threshold,data_dir,bw_filter)
	#ax = plt.gca()
	#ax.get_legend().remove()
	#plt.xlabel('Avalanche size S')
	#plt.ylabel('p(S)')
	plt.xlim(1,300)
	plt.ylim(fig_pS_ylim,1)
	plt.yticks([1e-4,1e-2,1])

	str_save = str_savepath + 'sub_b{:02d}_d{:02d}.pdf'.format(b,d)
	plt.savefig(str_save,bbox_inches="tight")
	plt.close()


	#Generates coarse sampled figure
	print('Fig1: generating coarse sampled figure')
	plt.figure(figsize=(fig_pS_size[0]/2.54,fig_pS_size[1]/2.54))
	plot_compare_states('coarse',b,d,reps,threshold,data_dir,bw_filter)
	#plt.legend(loc=3, frameon=False)
	# ax = plt.gca()
	# ax.get_legend().remove()
	#plt.xlabel('Avalanche size S')
	#plt.ylabel('p(S)')
	plt.yscale('log')
	plt.xscale('log')
	plt.xlim(1,300)
	plt.ylim(fig_pS_ylim,1)
	plt.yticks([1e-4,1e-3,1e-2,1e-1,1])
	plt.minorticks_on()

	str_save = str_savepath + 'coarse_b{:02d}_d{:02d}.pdf'.format(b,d)
	plt.savefig(str_save,bbox_inches="tight")
	plt.close()

def figure_2(data_dir,reps,bw_filter):

	#Sets path
	if bw_filter:
		str_savepath = 'figs/fig2_filtered/'
	else:
		str_savepath = 'figs/fig2_unfiltered/'
	if not os.path.exists(str_savepath):
		os.makedirs(str_savepath)

	#Parameters
	fig_size = [3,2.7]
	b_mav = [2,4]
	b_ps = [2,2,2]
	d_ps = [2,5,10]
	threshold = 3
	
	#Saves (m_av vs d) figure
	figure_mav(data_dir,b_mav,bw_filter)

	#Plots comparison of different (b,d) configurations


	#Plots states
	states = ['subcritical', 'reverberant','critical']
	colors = ['blue', 'green','red']
	#offset_list = [1,5e-3,2.5e-5]
	state_dict = states_parameters()

	for j in range(len(d_ps)):

		#Sets up figure
		plt.figure(figsize=(fig_size[0]/2.54,fig_size[1]/2.54))
		for i in range(len(states)):

			#Get parameters		
			m = state_dict[states[i]]['m']
			h = state_dict[states[i]]['h']

			#Plots distributions
			color_plot = color_picker(colors[i])
			str_leg = r'$Delta$t = {:02d}, '.format(2*b_ps[j]) + states[i]
			analysis.plot.sim_pS(m,h,d_ps[j],b_ps[j],'coarse',reps,str_leg,bw_filter,data_dir,threshold,color_plot,True)	

		#Plots comparison line
		X = np.arange(2,10)
		Y = np.power(X,-1.5)
		Y = 0.05*Y/np.sum(Y)
		plt.plot(X,Y,linestyle='--',color='black',label=None)

		#Adjusts figure
		#plt.xlabel('Avalanche size S')
		#plt.ylabel('p(S)')
		plt.yscale('log')
		plt.xscale('log')
		plt.xlim(1,100)
		plt.ylim(1e-5,1)
		plt.yticks([1e-5,1e-4,1e-3,1e-2,1e-1,1],['',r'$10^{-4}$','',r'$10^{-2}$','',r'$10^{0}$'])

		#plt.yticks([1e-4,1e-2,1],[])

		#Sets minor and major ticks
		ax = plt.gca()
		subs_loc = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8,0.9)
		locmin = matplotlib.ticker.LogLocator(base=10.0,subs=subs_loc,numticks=6)
		ax.yaxis.set_minor_locator(locmin)
		ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

		#Set size correctly
		set_size(fig_size[0]/2.54,fig_size[1]/2.54)

		#Saves figures
		str_save = str_savepath + 'ps_comparison_d{:02d}_b{:02d}.pdf'.format(d_ps[j], b_ps[j])
		plt.savefig(str_save,bbox_inches="tight")
		plt.close()

def figure_3(data_dir,d,reps,bw_filter, datatype):

	#Parameters
	b = [1,2,4,8]
	threshold = 3
	fig_pS_size = [4,4]
	fig_alpha_size = [2.5,2]
	states = ['poisson', 'subcritical', 'reverberant', 'critical']
	state_colors = {
		'poisson':'gray', 
		'subcritical':'blue', 
		'reverberant': 'green', 
		'critical':'red'}

	#Sets path
	if bw_filter:
		str_savepath = 'figs/fig3_filtered/'
	else:
		str_savepath = 'figs/fig3_unfiltered/'
	if not os.path.exists(str_savepath):
		os.makedirs(str_savepath)

	#Gets state variables
	state_dict = states_parameters()
	states_tau_dict = states_tau()
	
	#Plots every state
	for state in states:
		m = state_dict[state]['m']
		h = state_dict[state]['h']
		print('Fig2: generating  plots for state:' + state)

		#Gets line colors
		plt_color = color_gradient_rgba(state_colors[state],len(b))
		plt_color = convert_rgba_rgb(plt_color)

		#Plots and fits alpha = alpha(b) for non-poisson dynamics
		if state is not 'poisson':

			#Sets up figure
			plt.figure(figsize=(fig_alpha_size[0]/2.54,fig_alpha_size[1]/2.54))
			plt.xlabel(r'$\Delta$t (ms)')
			plt.ylabel(r'$\alpha$', rotation='horizontal')
			plt.xscale('log')
			plt.yscale('log')
			plt.minorticks_off() #solves tick bu... "intended changes" in pyplot
			plt.xlim(1,100)
			plt.ylim(1.1,2.1)			
			plt.xticks([1,10,100])
			plt.yticks([1.1,1.5,2.1],['1.1', '1.5', '2.1'])
			ax1 = plt.gca()
			ax1.spines['right'].set_visible(False)
			ax1.spines['top'].set_visible(False)

			#Plots and fits alpha(b)
			fit_exp, fit_err = analysis.plot.plot_alpha_bs(m,h,b,d,datatype,reps,bw_filter,data_dir,threshold, plt_color)

			#Saves alpha vs b fit to file
			str_save = str_savepath+'alpha_fitresults_' + datatype + '_d{:02d}.txt'.format(d)
			with open(str_save,'a') as file:
				str_fitdata = state + ': {:0.3f} +- {:0.3f}\n'.format(fit_exp,fit_err)
				file.write(str_fitdata)

			#Saves plot of alpha vs b
			set_size(fig_alpha_size[0]/2.54,fig_alpha_size[1]/2.54)
			str_save = str_savepath + 'fit_' + datatype + '_' + state + '_d{:02d}.pdf'.format(d)
			plt.savefig(str_save,bbox_inches="tight", transparent=True)
			plt.close()

		#Sets up pS figure
		plt.figure(figsize=(fig_pS_size[0]/2.54,fig_pS_size[1]/2.54))

		#Calculates tau
		#tau_rep = analysis.fitting.tau_sim_dataset(m,h,int(d),threshold,data_dir,bw_filter)
		#str_title = r'$\bar\tau$ = {:0.1f} ms'.format(np.mean(tau_rep))
		#plt.title(str_title)
		#str_title = state + r' ($\bar\tau$ = {:0.1f} ms)'.format(states_tau_dict[state]['tau'])
		str_title = state
		#plt.title(str_title)

		#Plots it
		for i in range(len(b)):

			bs_leg = r"$\Delta$t = {:d} ms".format(2*b[i])

			analysis.plot.sim_pS(m,h,d,b[i],datatype,reps,
				label_plot=bs_leg,
				bw_filter=bw_filter,
				data_dir=data_dir,
				plt_color=plt_color[i,:],
				plt_std=False)	

		if state == 'poisson' and datatype == 'coarse':
			plt.legend(frameon=False)

		#Sets up figure
		#plt.xlabel('Avalanche size S')
		#plt.ylabel('p(S)')
		plt.yscale('log')
		plt.xscale('log')
		if datatype == 'sub':
			plt.xlim(1,1000)
		else:
			plt.xlim(1,100)
		plt.ylim(1e-7,1)
		ytick_str = ['',r'$10^{-6}$','',r'$10^{-4}$','',r'$10^{-2}$','',r'$10^{0}$']
		plt.yticks([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1],ytick_str)

		#Sets minor and major ticks
		ax = plt.gca()
		subs_loc = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8,0.9)
		locmin = matplotlib.ticker.LogLocator(base=10.0,subs=subs_loc,numticks=8)
		ax.yaxis.set_minor_locator(locmin)
		ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())	

		#Saves plot
		set_size(fig_pS_size[0]/2.54,fig_pS_size[1]/2.54)
		str_save = str_savepath + datatype + '_' + state + '_d{:02d}.pdf'.format(d)
		plt.savefig(str_save,bbox_inches="tight")
		plt.close()

def figure_corr(data_dir,b,d,reps, threshold):

	#Parameters
	states = ['critical','reverberant','subcritical','poisson']
	state_labels = ['Critical', 'Reveberating', 'Subcritical', 'Poisson']
	colors = ['red', 'green','blue', 'gray']

	thresholds = [2,3]
	figSize = [6,2.5]
	plot_types = ['corr', 'rate']

	#Save path
	str_savepath = 'figs/correlation/'
	if not os.path.exists(str_savepath):
		os.makedirs(str_savepath)

	#Gets m and h for each state
	state_dict = states_parameters()

	#Plots for all thresholds
	for threshold in thresholds:
		#Plots data for both correlation and rate
		for plot_type in plot_types:

			#Sets up figure
			plt.figure(figsize=(figSize[0]/2.54,figSize[1]/2.54))		

			if plot_type == 'corr':
				plt.ylabel('Correlation')
				plt.ylim((0,0.7))
			else:
				plt.ylabel('Event rate (Hz)')

			#Plots all barplots for the corr plot
			for i in range(len(states)):
				color_plot = color_picker(colors[i])
				m = state_dict[states[i]]['m']
				h = state_dict[states[i]]['h']
				analysis.plot.sim_corr(m,h,d,b,threshold,reps,data_dir,type=plot_type, loc=i+1,color_state=color_plot)

			#Sets up ticks	
			ax = plt.gca()				
			plt.xticks(ticks=np.arange(1,len(states)+1),labels=state_labels)	
			ax.set_xticklabels(state_labels, fontsize=8)

			#Saves plot
			set_size(figSize[0]/2.54,figSize[1]/2.54)
			str_save = str_savepath + plot_type + '_d{:02d}_b{:2d}_th{:0.1f}_rep{:02d}.pdf'.format(d, b, threshold, reps)
			plt.savefig(str_save,bbox_inches="tight")
			plt.close()		

def figure_corr_d(data_dir,b,reps, threshold):
	#Parameters
	d_list = np.arange(1,11)
	pass

if __name__ == "__main__":

	parameters = parametersDefault()

	#Plots selected figures
	if parameters.fig == '1':
		figure_1(
			data_dir=parameters.datafolder,
			b=parameters.b,
			d=parameters.d,
			bw_filter=parameters.bw_filter, 
			reps=parameters.reps)

	elif parameters.fig == '2':
		figure_2(
			data_dir=parameters.datafolder,
			reps=parameters.reps,
			bw_filter = parameters.bw_filter)

	elif parameters.fig == '3':
		figure_3(
			data_dir=parameters.datafolder,
			d=parameters.d,
			bw_filter=parameters.bw_filter, 
			reps=parameters.reps,
			datatype = parameters.datatype)

	elif parameters.fig == 'corr':
		figure_corr(data_dir=parameters.datafolder,
			d=parameters.d,
			b=parameters.b,
			reps=parameters.reps,
			threshold=parameters.threshold)
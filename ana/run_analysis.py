# -*- coding: utf-8 -*-
# @Author: Joao PN
# @Date:   2019-03-25 16:45:25
# @Last Modified by:   joaopn
# @Last Modified time: 2019-07-23 02:34:51

from analysis import avalanche, plot, fitting, parser
import numpy as np
import matplotlib.pyplot as plt
import os, argparse
import glob
import h5py
import powerlaw

def parametersDefault():

	#default Parameters
	binsizeDefault="1,2,4"
	thresholdDefault = 3
	repsDefault = 1
	datatypeDefault = 'coarse'
	datafolderDefault = 'dat/'
	modeDefault = 'save_plot'
	bw_filterDefault = True

	#Parse input
	parser = argparse.ArgumentParser()

	#Add single parameters
	parser.add_argument("--mode",
		type=str,   nargs='?', const=1, default=modeDefault)
	parser.add_argument("-t", "--threshold",
		type=float, nargs='?', const=1, default=thresholdDefault)
	parser.add_argument("--reps",
		type=int,   nargs='?', const=1, default=repsDefault)
	parser.add_argument("--datatype",
		type=str,   nargs='?', const=1, default=datatypeDefault)
	parser.add_argument("--datafolder",
		type=str,   nargs='?', const=1, default=datafolderDefault)
	parser.add_argument("--datamask",
		type=str,   nargs='?', const=1, default=None)
	parser.add_argument("--bw_filter",
		type=bool,  nargs='?', const=1, default=bw_filterDefault)


	#Adds string of binsizes
	parser.add_argument("-b","--binsize",type=str,nargs='?',const=1,default=binsizeDefault)
	args = parser.parse_args()
	args.binsize = [int(item) for item in args.binsize.split(',')]

	return args

def save_plot(data_dir,filename,threshold,datatype,reps,binsize,bw_filter):

	#Save location
	if bw_filter:
		fig_dir = data_dir + 'plot_filtered/'
	else:
		fig_dir = data_dir + 'plot_unfiltered/'

	#Plots results
	plt.figure()
	plt.title(filename)
	plt.xlabel('S')
	plt.ylabel('p(S)')
	plt.yscale('log')
	plt.xscale('log')
	plt.xlim(1,1e3)
	plt.ylim(1e-5,1)
	X = np.arange(5,64)
	Y = np.power(X,-1.5)
	Y = 5*Y/np.sum(Y)
	plt.plot(X,Y,linestyle='--',color='black',label=r'$p(S)~S^{-1.5}$')

	#Runs analysis for all binsizes
	for bs in binsize:
		print(filename + ' ,running for b = {}'.format(bs))

		S_list = []
		for rep in range(reps):

			#Creates filepath
			filepath = data_dir + filename + '_r{:02d}.hdf5'.format(rep)

			#Analyzes rep
			data_thresholded = avalanche.analyze_sim_raw(
				filepath=filepath,
				threshold=threshold,
				datatype=datatype,
				bw_filter=bw_filter
				)
			data_binned = avalanche.bin_data(data=data_thresholded,binsize=bs)

			S_list.append(avalanche.get_S(data_binned))

		plot.pS_mean(S_list,label='bs = {}'.format(bs))

	#Creates save dir
	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	#Saves fig
	str_fig = fig_dir + datatype + '_pS_' + filename + '_th{}'.format(threshold) + '.png'
	plt.savefig(str_fig)
	plt.close()

def save_threshold(data_dir,filename,threshold,reps,bw_filter,timesteps=None):
	"""Thresholds datasets in [data_dir], numbered in files ending with '_rep%2d.hdf5'. Results are saved to [data_dir]/thresholded_[filtered/unfiltered]/filename.hdf5, where each columnn corresponds to a repetition.
	
	Args:
	    data_dir (str): Location to search for hdf5 datasets
	    filename (str): dataset name format
	    threshold (float): Threshold to use
	    reps (int): Number of repetitions in the dataset
	    bw_filter (bool): Whether to bandpass the coarse signal with a butterworth 4th order filter
	    timesteps (None, optional): Timeseries length to use, in timesteps.
	"""

	#Definitions
	if bw_filter:
		dir_threshold = 'thresholded_filtered/'
	else:
		dir_threshold = 'thresholded_unfiltered/'
	data_save = data_dir + dir_threshold
	str_savefile = data_save + filename + '_th{:0.1f}.hdf5'.format(threshold)
	datatypes = ['coarse', 'sub']

	#Gets timesteps
	if timesteps is None:
		filepath_0 = data_dir + filename + '_r00.hdf5'
		file_raw = h5py.File(filepath_0,'r')
		timesteps = file_raw['data/activity'].shape[0]
		file_raw.close()

	#Creates path to .h5 save file
	if not os.path.exists(data_save):
		os.makedirs(data_save)
	file = h5py.File(str_savefile,'w') #Overwrites old file


	#Saves thresholded data ('coarse' and 'sub')
	for datatype in datatypes:

		#Creates .h5 dataset
		file.create_dataset(datatype,shape=(reps,timesteps),dtype=int,chunks=(1,timesteps),compression=1,maxshape=(None,timesteps))
		file[datatype].attrs['threshold'] = threshold

		for rep in range(reps):

			#Creates filepath
			filepath = data_dir + filename + '_r{:02d}.hdf5'.format(rep)

			#Analyzes rep
			data_thresholded = avalanche.analyze_sim_raw(
				filepath=filepath,
				threshold=threshold,
				datatype=datatype,
				bw_filter=bw_filter
				)
			file[datatype][rep,:] = np.int32(data_thresholded)

	#Copies population activity data
	file.create_dataset('activity',shape=(reps,timesteps),dtype=float,chunks=(1,timesteps),compression=1,maxshape=(None,timesteps))
	for rep in range(reps):
		filepath = data_dir + filename + '_r{:02d}.hdf5'.format(rep)
		file_raw = h5py.File(filepath,'r')
		file['activity'][rep,:] = file_raw['data/activity'][:]
		file_raw.close()

	#Flushes and closes file
	file.close()

def save_ps(data_dir, filename, binsize, bw_filter, reps=None):
	"""Saves avalanche-size distribution data, using the thresholded .hdf5 files generated by save_thresholded. 
	
	Args:
	    data_dir (str): Location to search for hdf5 datasets
	    filename (str): dataset name format
	    binsize (list): Bin-sizes (in timesteps) to use for the avalanches
	    bw_filter (bool): Whether to bandpass the coarse signal with a butterworth 4th order filter
	    reps (None, optional): Number of repetitions to use
	"""
	#Definitions
	if bw_filter:
		dir_threshold = 'thresholded_filtered/'
		saveplot_dir = 'analyzed_filtered/'
	else:
		dir_threshold = 'thresholded_unfiltered/'
		saveplot_dir = 'analyzed_unfiltered/'

	datatypes = ['coarse', 'sub']

	#Loads file
	file_path = data_dir + dir_threshold + filename + '.hdf5'
	file = h5py.File(file_path,'r')

	#Gets reps from file
	if reps is None:
		reps = file['coarse'].shape[0]

	#Runs the analysis for each b
	for k in range(len(binsize)):

		#Analyzes both 'coarse' and 'sub'
		for datatype in datatypes:

			#Obtains S for each repetition
			S_list = []
			for rep in range(reps):

				#Loads data, analyzes and bins it
				data_thresholded = file[datatype][rep,:]
				data_binned = avalanche.bin_data(data=data_thresholded,binsize=binsize[k])
				S_list.append(avalanche.get_S(data_binned))

			#Gets largest avalanche from the list (+1 for zero_index)
			S_max = int(max([Si.max() for Si in S_list]) + 1)

			#Obtains pS
			pS = np.zeros((len(S_list),S_max))
			for i in range(len(S_list)):
				for j in range(S_max):
					pS[i,j] = np.sum(S_list[i]==j+1)
				pS[i,:] = pS[i,:]/np.sum(pS[i,:])

			#Obtains mean and STD
			pS_mean = np.mean(pS,axis=0)
			pS_std = np.std(pS,axis=0)
			X = np.arange(S_max)+1

			#Saves plot data
			str_savefolder = data_dir + saveplot_dir + filename + '_rep{:02d}/'.format(reps)
			if not os.path.exists(str_savefolder):
				os.makedirs(str_savefolder)
			str_savefile = 'pS_' + datatype + '_b{:02d}.tsv'.format(binsize[k])
			str_save = str_savefolder + str_savefile
			np.savetxt(str_save,(X,pS_mean,pS_std),delimiter='\t',header='S\tpS_mean\tpS_std')

def save_mav(data_dir, filename, d_list, binsize, threshold, bw_filter, reps=None):

	#Definitions
	if bw_filter:
		dir_threshold = 'thresholded_filtered/'
		save_dir = 'analyzed_filtered/branching_mav/'
	else:
		dir_threshold = 'thresholded_unfiltered/'
		save_dir = 'analyzed_unfiltered/branching_mav/'

	#Runs it for each binsize
	for k in range(len(binsize)):

		#Defines save variables
		mav_mean = np.zeros(len(d_list))
		mav_std = np.zeros(len(d_list))
		IED = np.array(d_list)

		#Runs analysis for each d
		for d_id in range(len(d_list)):

			#Loads file
			file_path = data_dir + dir_threshold + filename + "d{:02d}_th{:0.1f}.hdf5".format(d_list[d_id],threshold)
			file = h5py.File(file_path,'r')

			#Gets reps from file
			if reps is None:
				reps = int(file['coarse'].shape[0])

			#Obtains m_av for each repetition
			mav_reps = np.zeros(reps)
			for rep in range(reps):

				#Loads data, analyzes and bins it
				data_thresholded = file['coarse'][rep,:]
				data_binned = avalanche.bin_data(data=data_thresholded,binsize=binsize[k])
				mav_reps[rep] = fitting.m_avalanche(data_binned)


			#Obtains mean and STD
			mav_mean[d_id] = np.mean(mav_reps)
			mav_std[d_id] = np.std(mav_reps)

			print('b = {:d}, d = {:d}: m_av = {:0.3f} +- {:0.3f}'.format(binsize[k], d_list[d_id], mav_mean[d_id],mav_std[d_id]))

		#Saves data
		if not os.path.exists(data_dir + save_dir):
			os.makedirs(data_dir + save_dir)
		str_savefile = filename + 'b{:02d}_th{:0.1f}.tsv'.format(binsize[k], threshold)
		str_save = data_dir + save_dir + str_savefile
		np.savetxt(str_save,(IED,mav_mean,mav_std),delimiter='\t',header='d\tmav_mean\tmav_std')

def save_ps_alpha(data_dir, filename, binsize, bw_filter, reps=None, xmax=50):

	#Parameters
	timescale = 2

	#Definitions
	if bw_filter:
		dir_threshold = 'thresholded_filtered/'
		savedata_dir = 'analyzed_filtered/'
	else:
		dir_threshold = 'thresholded_unfiltered/'
		savedata_dir = 'analyzed_unfiltered/'

	datatypes = ['coarse', 'sub']
	

	#Loads file
	file_path = data_dir + dir_threshold + filename + '.hdf5'
	file = h5py.File(file_path,'r')

	#Gets reps from file
	if reps is None:
		reps = file['coarse'].shape[0]

	#Analyzes both 'coarse' and 'sub'
	for datatype in datatypes:

		#Defines alpha matrix
		alpha_fit = np.zeros((len(binsize), reps))

		#Runs the analysis for each b
		for k in range(len(binsize)):
			#Obtains S for each repetition
			for rep in range(reps):

				#Loads data, bins it and fits the avalanches
				data_thresholded = file[datatype][rep,:]
				data_binned = avalanche.bin_data(data=data_thresholded,binsize=binsize[k])
				S = avalanche.get_S(data_binned)
				fit = powerlaw.Fit(S, discrete=True, estimate_discrete=False, xmin=1, xmax=xmax)
				alpha_fit[k,rep] = fit.alpha				

		#Obtains mean and STD
		alpha_mean = np.mean(alpha_fit,axis=1)
		alpha_std = np.std(alpha_fit,axis=1)
		X = timescale*np.array(binsize)

		#Saves plot data
		str_savefolder = data_dir + savedata_dir + filename + '_rep{:02d}/'.format(reps)
		if not os.path.exists(str_savefolder):
			os.makedirs(str_savefolder)
		str_savefile = 'alpha_' + datatype + '.tsv'
		str_save = str_savefolder + str_savefile
		np.savetxt(str_save,(X,alpha_mean,alpha_std),delimiter='\t',header='b\talpha_mean\talpha_std')	

def save_corr(data_dir, filename, d_list, binsize, threshold, bw_filter,reps=None):

	#Parameters
	elec_base = 0
	bw_freqs = [0.1,200]
	fs = 500
	str_coarse = 'data/coarse'
	str_sub = 'data/sub'
	str_activity = 'data/activity'

	#Definitions
	save_dir = 'correlations/'

	#Parse input
	if type(binsize) is not list:
		binsize = [binsize]

	#Runs for every binsize
	for bs in binsize:

		#Runs it for every d
		for d in d_list:
			corr_coarse = np.zeros(reps)
			corr_sub = np.zeros(reps)
			rate_coarse = np.zeros(reps)
			rate_sub = np.zeros(reps)

			#Runs it for every rep
			for i in range(reps):

				#Loads coarse and sub data
				filepath = data_dir + filename + 'd{:02d}_r{:02d}.hdf5'.format(d,i)
				file = h5py.File(filepath,'r')
				timesteps = file[str_activity].shape[0]

				data_coarse = file[str_coarse][elec_base:elec_base+2,:]
				data_sub = file[str_sub][elec_base:elec_base+2,:]

				#Filters the signal
				if bw_filter:
					data_coarse[0,:] = avalanche.filter_bw_ch(data_coarse[0,:],bw_freqs,fs)
					data_coarse[1,:] = avalanche.filter_bw_ch(data_coarse[1,:],bw_freqs,fs)

				#Thresholds coarse data and bins data
				data_coarse_bin_0 = avalanche.bin_data(avalanche.threshold_ch(data_coarse[0,:], threshold), bs)
				data_coarse_bin_1 = avalanche.bin_data(avalanche.threshold_ch(data_coarse[1,:], threshold), bs)
				data_sub_bin_0 = avalanche.bin_data(avalanche.convert_timestamps(data_sub[0,:],timesteps), bs)
				data_sub_bin_1 = avalanche.bin_data(avalanche.convert_timestamps(data_sub[1,:],timesteps), bs)
														
				#Calculates correlations
				corr_coarse[i] = np.corrcoef(data_coarse_bin_0, data_coarse_bin_1)[0,1]
				corr_sub[i] = np.corrcoef(data_sub_bin_0, data_sub_bin_1)[0,1]

				#Calculates rates
				timescale_s = 2*1e-3*bs*data_coarse_bin_0.size
				rate_coarse[i] = (np.sum(data_coarse_bin_0) + np.sum(data_coarse_bin_1))/2/timescale_s
				rate_sub[i] = (np.sum(data_sub_bin_0) + np.sum(data_sub_bin_1))/2/timescale_s

				file.close()	

			print('d = {:02d}: coarse = {:0.2f} ({:0.2f} Hz), sub = {:0.2f} ({:0.2f} Hz)'.format(d,np.mean(corr_coarse), np.mean(rate_coarse),np.mean(corr_sub), np.mean(rate_sub[i])))

			#Saves results
			str_savefolder = data_dir + save_dir + '/'
			if not os.path.exists(str_savefolder):
				os.makedirs(str_savefolder)
			str_savefile = str_savefolder + filename + '_d{:02d}_b{:02d}_th{:0.1f}_rep{:02d}.tsv'.format(d,bs,threshold,reps)
			with open(str_savefile, 'w') as f:
				np.savetxt(f,(corr_coarse,corr_sub, rate_coarse, rate_sub),delimiter='\t',header='coarse\tsub\trate_coarse\trate_sub')	

if __name__ == "__main__":

	#Parse input
	args       = parametersDefault()
	mode       = args.mode
	binsize    = args.binsize
	threshold  = args.threshold
	reps       = args.reps
	datatype   = args.datatype
	datafolder = args.datafolder
	datamask   = args.datamask
	bw_filter  = args.bw_filter

	#Does the requested operation
	if mode == 'save_plot':
		dataset_list = parser.sim_find_unique(datafolder, datamask)
		print('Analyzing and plotting')
		for i in range(len(dataset_list)):
			save_plot(data_dir=datafolder,
				filename=dataset_list[i],
				binsize=binsize,
				threshold=threshold,
				datatype=datatype,
				reps=reps,
				bw_filter=bw_filter)

	elif mode == 'threshold':
		dataset_list = parser.sim_find_unique(datafolder, datamask)
		for i in range(len(dataset_list)):
			print('Thresholding and saving: ' + dataset_list[i])
			save_threshold(data_dir=datafolder,
				filename=dataset_list[i],
				threshold=threshold,
				reps=reps,
				bw_filter=bw_filter)

	elif mode == 'save_ps':
		dataset_list = parser.sim_find_thresholded(datafolder, bw_filter, datamask)
		for i in range(len(dataset_list)):
			print('Analysing thresholded data and saving pS: ' + dataset_list[i])
			save_ps(datafolder,
				filename=dataset_list[i],
				binsize=binsize,
				bw_filter=bw_filter)

	elif mode == 'save_ps_alpha':
		dataset_list = parser.sim_find_thresholded(datafolder, bw_filter, datamask)
		for i in range(len(dataset_list)):
			print('Analysing thresholded data and fitting pS: ' + dataset_list[i])
			save_ps_alpha(datafolder,
				filename=dataset_list[i],
				binsize=binsize,
				bw_filter=bw_filter)

	elif mode == 'save_mav':
		dataset_list, d_list = parser.sim_find_thresholded_no_d(datafolder, bw_filter, datamask)
		for i in range(len(dataset_list)):
			print('Analysing thresholded data: ' + dataset_list[i])
			save_mav(datafolder,
				filename=dataset_list[i],
				d_list=d_list,
				binsize=binsize,
				threshold=threshold,
				bw_filter=bw_filter)		

	elif mode == 'save_corr':
		dataset_list,d_list = parser.sim_find_unique_no_d(datafolder)
		for i in range(len(dataset_list)):
			print('Calculating timeseries correlations: ' + dataset_list[i])
			save_corr(
				data_dir=datafolder, 
				filename=dataset_list[i], 
				d_list=d_list, 
				binsize=binsize,
				threshold=threshold,
				bw_filter=bw_filter,
				reps=reps)
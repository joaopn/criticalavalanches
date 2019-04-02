"""
Module for the avalanche analysis of MEA datasets.

"""
# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-22 12:54:07
# @Last Modified by:   joaopn
# @Last Modified time: 2019-04-02 00:43:53

import numpy as np
import h5py
import os
from scipy.signal import butter, lfilter

def threshold_ch(data, threshold):

	#Demeans data
	data = data - np.mean(data)

	#Defines threshold
	th = threshold*np.std(data)

	#Thresholds signal
	data[data<th] = 0

	#Finds crossings of the signal to positive
	id_cross=np.where(np.sign(data[:-1]) != np.sign(data[1:]))[0] + 1
	id_cross_plus = id_cross[data[id_cross]>0]

	#Defines output array
	data_thresholded = np.zeros(data.shape)
	
	#For every positive excursion, finds max
	for id in range(len(id_cross_plus)-1):
		id_max = id_cross_plus[id]+np.argmax(data[id_cross_plus[id]:id_cross_plus[id+1]])
		data_thresholded[id_max] = 1

	return data_thresholded

def convert_timestamps(data,timesteps):

	data_converted = np.zeros(timesteps)
	data_converted[data] = 1
	return data_converted

def bin_data(data,binsize):

	timesteps = data.size
	bin_array = np.arange(0,timesteps,binsize)
	data_binned = np.zeros(bin_array.size)

	for i in range(data_binned.size-1):
		data_binned[i] = np.sum(data[bin_array[i]:bin_array[i+1]])

	return data_binned

def filter_bw_ch(data,freqs=[0.1,200],fs=500):

	#Parameters
	order = 4

	#Filters signal (simple butterworth)
	b, a = butter(order, [2*freqs[0]/fs, 2*freqs[1]/fs], btype='band')
	data_filt = lfilter(b, a, data)

	return data_filt
	
def get_S(data):

	#Finds crossings of the signal to positive
	id_cross=np.where(np.sign(data[:-1]) != np.sign(data[1:]))[0] + 1
	id_cross_plus = id_cross[data[id_cross]>0]


	#Calculates size S of avalanches
	n_avalanches = id_cross_plus.size

	S = np.zeros(n_avalanches-1)

	for i in range(n_avalanches-1):
		S[i] = np.sum(data[id_cross_plus[i]:id_cross_plus[i+1]])

	return S

def analyze_sim_raw(
	filepath,
	threshold,
	datatype,
	bw_filter,
	timesteps=None,
	channels=None
	):

	#Parameters
	data_dir = 'data/'
	fs = 500
	bw_freqs = [0.1, 200]

	#Gets timesteps and number of channels
	file = h5py.File(filepath,'r')

	if timesteps is None:
		timesteps = file[data_dir + 'activity'].shape[0]
	if channels is None:
		channels = int(file['meta/num_elec'][:])

	#Analyses channel by channel
	data_th = np.zeros(timesteps)
	for ch in range(channels):
		
		#Loads data
		data_ch = file[data_dir + datatype][ch,:]

		#Analyzes sub and coarse channel data accordingly
		if datatype == 'coarse':
			if bw_filter:
				data_ch = filter_bw_ch(data_ch,bw_freqs,fs)		
			data_th_ch = threshold_ch(data_ch,threshold)			
		elif datatype == 'sub':
			data_th_ch = convert_timestamps(data_ch,timesteps)

		data_th = data_th + data_th_ch

	return data_th

def save_sim_pS(data_dir,dataset,binsize,reps=None, bw_filter=False):

	#Definitions
	if bw_filter:
		dir_threshold = 'thresholded_filtered/'
		saveplot_dir = 'analyzed_filtered/'
	else:
		dir_threshold = 'thresholded_unfiltered/'
		saveplot_dir = 'analyzed_unfiltered/'

	datatypes = ['coarse', 'sub']

	#Loads file
	file_path = data_dir + dir_threshold + dataset + '.hdf5'
	file = h5py.File(file_path,'r')

	#Gets reps from file
	if reps is None:
		reps = file['coarse'].shape[0]

	#Analyzes both 'coarse' and 'sub'
	for datatype in datatypes:

		#Obtains S for each repetition
		S_list = []
		for rep in range(reps):

			#Loads data and bins it
			data_thresholded = file[datatype][rep,:]
			data_binned = bin_data(data=data_thresholded,binsize=binsize)
			S_list.append(get_S(data_binned))

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
		str_savefolder = data_dir + saveplot_dir + dataset + '_rep{:02d}/'.format(reps)
		if not os.path.exists(str_savefolder):
			os.makedirs(str_savefolder)
		str_savefile = 'pS_' + datatype + '_b{:02d}.tsv'.format(binsize)
		str_save = str_savefolder + str_savefile
		np.savetxt(str_save,(X,pS_mean,pS_std),delimiter='\t',header='S\tpS_mean\tpS_std')
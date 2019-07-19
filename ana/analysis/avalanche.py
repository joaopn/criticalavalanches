# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-22 12:54:07
# @Last Modified by:   Joao
# @Last Modified time: 2019-07-06 13:23:07

"""

Module for the avalanche analysis of MEA datasets.


"""

import numpy as np
import h5py, os
from scipy.signal import butter, lfilter

def threshold_ch(data, threshold):

	#Demeans data
	data = data - np.mean(data)

	#Defines threshold
	th = threshold*np.std(data)

	#Finds crossings of the signal to positive
	id_cross=np.where(np.sign(data[:-1]) != np.sign(data[1:]))[0] + 1
	id_cross_plus = id_cross[data[id_cross]>0]

	#Defines output array
	data_thresholded = np.zeros(data.shape)
	
	#For every positive excursion, finds max and add and event if larger than th
	for id in range(len(id_cross_plus)-1):
		id_max = id_cross_plus[id]+np.argmax(data[id_cross_plus[id]:id_cross_plus[id+1]])
		if data[id_max] > th:
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
	id_cross=	np.where(np.sign(data[:-1]) != np.sign(data[1:]))[0] + 1
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
		channels = file[data_dir+datatype].shape[0]

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
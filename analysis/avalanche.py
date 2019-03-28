"""
Module for the avalanche analysis of MEA datasets.

"""
# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-22 12:54:07
# @Last Modified by:   Joao PN
# @Last Modified time: 2019-03-28 16:06:57

import numpy as np
import h5py

def threshold_data(data, threshold):
	import numpy as np


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

def run_analysis(
	filepath,
	threshold,
	channels,
	binsize,
	datatype
	):

	#Parameters
	data_dir = 'data/'

	#Gets timesteps
	file = h5py.File(filepath,'r')
	timesteps = file[data_dir + 'activity'].shape[0]

	#Analyses channel by channel
	data_sum = np.zeros(timesteps)
	for ch in range(channels):
		
		#Loads data
		data_ch = file[data_dir + datatype][ch,:]
		print(np.sum(data_ch))

		#Analyzes sub and coarse channel data accordingly
		if datatype == 'coarse':		
			data_th = threshold_data(data_ch,threshold)			
		elif datatype == 'sub':
			data_th = convert_timestamps(data_ch,timesteps)

		data_sum = data_sum + data_th
	
	data_binned = bin_data(data_sum,binsize)
	print('Total events: {:d}'.format(int(np.sum(data_binned))))

	return data_binned
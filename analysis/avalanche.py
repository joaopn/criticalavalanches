"""
Module for the avalanche analysis of MEA datasets.

"""
# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-22 12:54:07
# @Last Modified by:   Joao PN
# @Last Modified time: 2019-03-22 19:42:24

def threshold_data(data, threshold):
	import numpy as np


	#Demeans data
	data = data - np.mean(data)

	#Defines threshold
	th = threshold*np.std(data)

	#Thresholds signal
	data[data<th] = 0

	#Finds crossings of the signal to positive
	id_cross=np.where(np.sign(a[:-1]) != np.sign(a[1:]))[0] + 1
	id_cross_plus = id_cross[id_cross > 0]

	#Defines output array
	data_thresholded = np.zeros(len(id_cross_plus),1)
	
	#For every positive excursion, finds max
	for id in range(len(crossing_plus)):
		id_max = id_cross_plus[id]+np.argmax(data[id_cross_plus[id]:id_cross_plus[id+1]-1])
		data_thresholded[id_max] = True

	return data_thresholded


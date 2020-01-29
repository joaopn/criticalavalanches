# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-04-01 01:44:18
# @Last Modified by:   joaopn
# @Last Modified time: 2019-12-21 23:18:20

import numpy as np
import h5py
import powerlaw
from scipy.optimize import curve_fit, leastsq

def tau_linear(data,deltaT = 2):

	#Defines X and Y
	X = np.array(data[:-1],dtype=float) #A(t)
	Y = np.array(data[1:],dtype=float) #A(t+1)

	m = np.cov(X,Y)[0,1]/np.var(X)

	if m > 0:
		tau = -deltaT/np.log(m)
	else:
		tau = 0

	return tau

def tau_sim_dataset(m,h,d,threshold,data_dir,bw_filter):
	
	#Sets up filepaths
	if bw_filter:
		dir_threshold = 'thresholded_filtered/'
	else:
		dir_threshold = 'thresholded_unfiltered/'


	#Parse input
	if type(d) not in [list, np.ndarray]:
		d = [d]

	# Calculates tau for every d and rep
	tau_all = np.zeros(0)
	for d_i in d:
		#Sets up filepath
		filename = 'm{:0.5f}_h{:0.3e}_d{:02d}_th{:0.1f}.hdf5'.format(m,h,d_i,threshold)		
		str_load = data_dir + dir_threshold + filename

		#Loads file
		file = h5py.File(str_load,'r')
		reps = int(file['activity'].shape[0])

		tau_rep = np.zeros(reps)

		for i in range(reps):
			tau_rep[i] = tau_linear(file['activity'][i,:])

		tau_all = np.concatenate((tau_all,tau_rep))
	
	tau_mean = np.mean(tau_all)
	tau_std = np.std(tau_all)

	return tau_mean,tau_std

def m_avalanche(data):
	"""Estimates m_{av}
	
	Calculates the branching parameter m_{av} estimated from avalanches of a thresholded and binned timeseries.
	:param data: thresholded 1xtimesteps timeseries
	:type data: float
	"""

	#Finds timesteps with activity
	act_id = np.nonzero(data)[0]

	#Removes last timestep if active
	if act_id[-1] == len(data)-1:
		act_id = np.delete(act_id,-1)

	#Calculates A(t+1)/A(t) for every point of activity and averages it
	R = data[act_id+1]/data[act_id]
	m_av = np.mean(R)

	return m_av

def powerlaw(X,Y,Yerr, loglog=False):
	#Parameters
	kwargs = {'maxfev': 100000}

	if loglog:
		def lin(x,a,b):
			return a*x+b

		bool_keep = np.less(0,Y)
		X = X[bool_keep]
		Y = Y[bool_keep]
		Yerr = Yerr[bool_keep]

		X = np.log10(X)
		Y = np.log10(Y)
		Yerr = np.divide(Yerr,Y)

		results, pcov = curve_fit(lin,X,Y, **kwargs)
		#results, pcov = curve_fit(lin,X,Y, method='lm', p0=[0,2], **kwargs)
		lin_coef = 10**results[1]
		fit_exp = results[0]

	else:
		def pl(x,a,b):
			return a*np.power(x,b)

		#Fits courve with a LM algorithm
		#results, pcov = curve_fit(pl,X,Y,sigma=Yerr, method='lm', p0=[2,0.1], **kwargs)
		results, pcov = curve_fit(pl,X,Y, **kwargs)
		lin_coef = results[0]
		fit_exp = results[1]

	#Gets error
	fit_err_all = np.sqrt(np.diag(pcov))

	#Gets relevant variables

	fit_err = fit_err_all[1]

	return fit_exp, fit_err, lin_coef

def shape_collapse(shape_list, min_d, min_rep):

	from scipy.interpolate import interp1d, InterpolatedUnivariateSpline	
	from scipy.optimize import minimize

	#Definitions
	interp_points = 1000
	gamma_x0 = 0.5
	opt_bounds = (-1,5)

	#Flattens list of reps
	flat_list = np.array([item for sublist in shape_list for item in sublist])

	#List of avalanche sizes
	shape_size = np.zeros(len(flat_list))
	for i in range(len(flat_list)):
		shape_size[i] = flat_list[i].size

	max_size = shape_size.max()

	#Avalanche size count
	shape_count,_ = np.histogram(shape_size,bins=np.arange(0,max_size+2))

	#Censors flat_list
	# censor_d_keep = list(shape_size>min_d)
	# censor_rep_keep = [True for i in range(len(flat_list))]
	# shape_count_rem = shape_count < min_rep
	# for i in range(len(flat_list)):
	# 	if shape_count_rem[int(shape_size[i])]:
	# 		censor_rep_keep[i] = False
	# censor_all =  [a and b for a, b in zip(censor_d_keep, censor_rep_keep)]
	# flat_list_censored = flat_list[censor_all]

	#Censors data by size
	censor_d_keep = np.arange(0,max_size+1) >= min_d
	censor_rep_keep = shape_count >= min_rep
	censor_index =  np.where([a and b for a, b in zip(censor_d_keep, censor_rep_keep)])[0]

	#Defines average size matrix
	average_shape = np.zeros((censor_index.size, interp_points))

	#Defines bottom interpolation range from data, to prevent extrapolation bias
	x_min = 1/censor_index[0]
	x_range = np.linspace(x_min,1,num=interp_points)

	#Averages shape for each duration and interpolates results
	for i in range(len(censor_index)):

		#Calculates average shape
		size_i = censor_index[i]
		avg_shape_i_y = np.mean(flat_list[shape_size==size_i])
		avg_shape_i_x = np.arange(1,size_i+1)/size_i

		#Interpolates results
		fx = InterpolatedUnivariateSpline(avg_shape_i_x,avg_shape_i_y)
		average_shape[i,:] = fx(x_range)

	#Error function for optimization	
	def _error(gamma_shape, *params):
		average_shape, censor_index = params
		#shape_scaled = np.divide(average_shape.T,np.power(censor_index, gamma_shape)).T
		shape_scaled = np.zeros((censor_index.size, interp_points))
		for i in range(censor_index.size):
			shape_scaled[i,:] = average_shape[i,:]/np.power(censor_index[i],gamma_shape)

		err = np.mean(np.var(shape_scaled, axis=0))/np.power((np.max(np.max(shape_scaled))-np.min(np.min(shape_scaled))),2)
		return err

	#Minimizes error
	minimize_obj = minimize(_error, x0=[gamma_x0], args=(average_shape,censor_index), bounds=[opt_bounds])

	return minimize_obj.x - 1
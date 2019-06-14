"""
@author: Jérémie Gaffarel & Rudi Smits - TU Delft - Faculty of Aerospace Engineering - BSc 2 (2018-2019)

This module is meant to filter the noisy measurements of MAVs, by the mean of an Extended Kalman FIlter (EKF)
"""

import numpy as np

def get_h(x):
	# Observation model
	Pn = -68
	h = np.zeros((8,1))
	h[0] = Pn - 2 * 10 * np.log10(np.sqrt(x[0] ** 2 + x[1] ** 2 + x[8] ** 2))
	h[1] = x[2]
	h[2] = x[3]
	h[3] = x[6]
	h[4] = np.cos(x[6] - x[7]) * x[4] - np.sin(x[6] - x[7]) * x[5]
	h[5] = np.sin(x[6] - x[7]) * x[4] + np.cos(x[6] - x[7]) * x[5]
	h[6] = x[7]
	h[7] = x[8]
	return h

def get_H(x):
	# Jacobian matrix of get_h()
	H = np.zeros((8,9))
	H[0][0] = -20 * x[0] / (np.log(10) * (x[0] ** 2 + x[1] ** 2 + x[8] ** 2))
	H[0][1] = -20 * x[1] / (np.log(10) * (x[0] ** 2 + x[1] ** 2 + x[8] ** 2))
	H[0][8] = -20 * x[8] / (np.log(10) * (x[0] ** 2 + x[1] ** 2 + x[8] ** 2))
	H[1][2] = 1
	H[2][3] = 1
	H[3][6] = 1
	H[4][4] = np.cos(x[6] - x[7])
	H[4][5] = -np.sin(x[6] - x[7])
	H[4][6] = np.sin(x[7] - x[6]) * x[4] - np.cos(x[7] - x[6]) * x[5]
	H[4][7] = np.sin(x[6] - x[7]) * x[4] + np.cos(x[6] - x[7]) * x[5]
	H[5][4] = np.sin(x[6] - x[7])
	H[5][5] = np.cos(x[6] - x[7])
	H[5][6] = np.cos(x[7] - x[6]) * x[4] + np.sin(x[7] - x[6]) * x[5]
	H[5][7] = -np.cos(x[6] - x[7]) * x[4] + np.sin(x[6] - x[7]) * x[5]
	H[6][7] = 1
	H[7][8] = 1
	return H

def get_A(dt):
	# Model matrix
	A = np.identity(9)
	A[0][2] = -dt
	A[0][4] = dt
	A[1][3] = -dt
	A[1][5] = dt
	return A

def init_P():
	# Covariance matrix
	return np.identity(9)

def init_Q(s_p=0.1, s_v=0.5, s_psi=0.5, s_z=0.5):
	# Process noise matrix
	Q = np.zeros((9,9))
	Q[0:2, 0:2] = np.identity(2) * s_p ** 2
	Q[2:6, 2:6] = np.identity(4) * s_v ** 2
	Q[6:8, 6:8] = np.identity(2) * s_psi ** 2
	Q[8:9, 8:9] = np.identity(1) * s_z ** 2
	return Q

def init_R(s_m=5, s_v=0.2, s_psi=0.2, s_z=0.2):
	# Measurement noise matrix
	R = np.zeros((8, 8))
	R[0:1, 0:1] = np.identity(1) * s_m ** 2
	R[1:5, 1:5] = np.identity(4) * s_v ** 2
	R[5:7, 5:7] = np.identity(2) * s_psi ** 2
	R[7:8, 7:8] = np.identity(1) * s_z ** 2
	return R

def kalman_filter(time, rssi, gt, gt_rel, s_q=[0.1, 0.5], s_r=[5,0.2], optimisation=False):
	""" Kalman filter
	* Inputs:
		- time: list containing the time at wich each measurement was taken
		- rssi: noisy RSSI between the two MAVs
		- gt: ground-true measurements
		- gt_rel: ground-true relative measurements
		- s_q = [s_a, s_b]: (optional) s_a and s_b are here the Standard Deviations (SD) used in the Q matrix of the filter
		- s_r = [s_a, s_b]: (optional) same for R matrix
		- optimisation: boolean: (optional) if True, no additional noise is added to the ground-true measurements
	* Outputs:
		- x_filter: dict containing x/y/z_rel, vx/vy, vx/vy_other, psi, psi_other
		- d/b/vel/rssi_f: filtered range/bearing/absolute velocity/rssi
		- d/b/vel_g: ground-true range/bearing/velocity
		- d_unf: unfiltered range
	"""
	n = len(time)
	x_filter = []
	rssi_out = []
	# x = [x, y, vx, vy, vx_other, vy_other, psi, psi_other, z]
	# Initialised vx and vy shall be non-zero
	x = np.mat([[0.], [0.], [1.], [1.], [0.], [0.], [0.], [0.], [0.]])
    	# Initialise matrices; R, Q are constant
	P = np.mat(init_P())
	Q = init_Q(s_p=s_q[0], s_v=s_q[1], s_psi=s_q[1], s_z=s_q[1])
	R = init_R(s_m=s_r[0], s_v=s_r[1], s_psi=s_r[1], s_z=s_r[1])

	x_filter = np.zeros((n, 9))
	Z = np.empty((8,n))

	noise = []
	# Create *random* noise lists, with a normal distribution (mean=0, SD=0.2)
	for i in range(8):
		if not optimisation:
			noise.append(np.random.normal(0, 0.2, n))
		else:
			noise.append(np.zeros(n))	# Do not use noise while optimising

	# Initialise the measurement array
	# Add the random noise to the ground truth measurements. Because the Kalman Filter needs noise to work properly
	Z[0] = rssi				# RSSI
	Z[1] = gt[0]["vx"] + noise[1]		# v_x MAV 1
	Z[2] = gt[0]["vy"] + noise[2]		# v_y MAV 1
	Z[3] = gt[0]["psi"] + noise[3]		# psi 
	Z[4] = gt[1]["vx"] + noise[4]		# vx_other
	Z[5] = gt[1]["vy"] + noise[5]		# vx_other
	Z[6] = gt[1]["psi"] + noise[6]		# psi_other
	Z[7] = noise[7]				# z_rel (assumed to be 0 at all times)
	Z = Z.T

	# Run the filter steps for each measurement
	for t in range(n):
		# Compute the time step
		if t > 0:
			dt = time[t] - time[t-1]
		else:
			dt = 0.3

		z = Z[t].reshape((8,1))			# get the measurement
		A = get_A(dt)				# get the model, with the correct time step

		# Prediction step
		x_p = A * x				# state prediction
		z_p = get_h(x_p)			# measurement predicition

		P = A * P * A.T + Q			# covariance matrix prediction
		
		# Update step
		H = get_H(x_p)				# compute the Jacobian of the predicted state matrix
		P1 = P * H.T
		K = P1 * np.linalg.inv(H * P1 + R)	# compute the Kalman Gain
		
		x = x_p + K * (z - z_p)			# update the state
		P = (np.identity(9) - K * H) * P	# update the covariance matrix

		x_filter[t] = x.reshape(9)		# save the filtered measurement
		
	col = ["x_rel", "y_rel", "vx", "vy", "vx_other", "vy_other", "psi", "psi_other", "z_rel"]
	x_filter = dict(zip(col, x_filter.T))
	d_g = np.sqrt(gt_rel["x"] ** 2 + gt_rel["y"] ** 2 + gt_rel["z"] ** 2)			# Ground true range
	b_g = np.arctan2(gt_rel["y"], gt_rel["x"])						# Ground true bearing
	d_f = np.sqrt(x_filter["x_rel"] ** 2 + x_filter["y_rel"] ** 2 + x_filter["z_rel"] ** 2)	# Filtered range
	b_f = np.arctan2(x_filter["y_rel"], x_filter["x_rel"])					# Filtered bearing
	d_unf = np.power([10] * len(time), (rssi + 68) / (-20))					# Unfiltered range
	vel_f = np.sqrt((x_filter["vx"]) ** 2 + (x_filter["vy"]) ** 2)				# Filtered velocity
	vel_g = np.sqrt((gt[1]["vx"] - gt[0]["vx"]) ** 2 + (gt[1]["vy"] - gt[0]["vy"]) ** 2)	# Ground true velocity
	rssi_f = -68 - 20 * np.log10(d_f)														# Filtered RSSI

	return x_filter, d_g, d_f, d_unf, b_g, b_f, vel_g, vel_f, rssi_f

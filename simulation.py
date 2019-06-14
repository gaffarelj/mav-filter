"""
@author: Jérémie Gaffarel & Rudi Smits - TU Delft - Faculty of Aerospace Engineering - BSc 2 (2018-2019)

This module is meant to simulate MAV flights
"""

import numpy as np

def simulate(time, simcase=1):
	"""Simulate
	* Inputs:
		- time list
		- simulation case: int from 1 to 10
	* Outputs:
		- rssis: simulated rssi between the two MAVs, with noise
		- gt: ground truth simulated positions of the two MAVs
		- mav: simulated positions of the two MAVs
		- mav_re: relative simulated positions of the two MAVs
	"""
	n = len(time)

	dt = time[1] - time[0]						# Timestep is calculated
	mav = []									# Empty list is initialised for dict with the x and y velocities + positions of the two mav's
	vx, vy, x0, y0 = create_vx(time, simcase)	# Velocities and initial positions are generated for the given simulation case

	# Initialistion of the dict containing the measured positions
	mav.append({"vx": vx[0], "vy": vy[0], "x":[x0[0]] * n, "y":[y0[0]] * n})
	mav.append({"vx": vx[1], "vy": vy[1], "x":[x0[1]] * n, "y":[y0[1]] * n})	
	mav_rel = {"x":np.ones(n), "y":np.ones(n), "z":np.zeros(n)}					
	rssis = np.ones(n)    # array is initialised for the measured RSSI
	dists = np.zeros(n)   # array is initialised for the distance between the two MAVs

	for i in range(1, n):
		x0n = mav[0]["x"][i - 1] + mav[0]["vx"][i] * dt		# x position at time t is the pos at t + velocity*dt
		y0n = mav[0]["y"][i - 1] + mav[0]["vy"][i] * dt		# same for y
		# same for the second MAV
		x1n = mav[1]["x"][i - 1] + mav[1]["vx"][i] * dt
		y1n = mav[1]["y"][i - 1] + mav[1]["vy"][i] * dt
		x_rel = x1n - x0n									# relative velocity between the two MAVs in the x direction
		y_rel = y1n - y0n									# same in y direction
		dist = np.sqrt(x_rel ** 2 + y_rel ** 2)				# relative distance between the two MAVs
		# If the distance is nonzero, the RSSI can be calculated with the log function. Otherwise it is set to -88 [dB]
		if dist != 0:
			rssi = -68 - 20 * np.log10(dist)
		else:
			rssi = -88
		rssi += np.random.normal(0, 1.5)					# random noise is added to the computed RSSI

		# Computed positions are saved in the dict
		mav[0]["x"][i] = x0n
		mav[0]["y"][i] = y0n
		mav[1]["x"][i] = x1n
		mav[1]["y"][i] = y1n
		mav_rel["x"][i] = x_rel
		mav_rel["y"][i] = y_rel
		rssis[i] = rssi
	gt = []
	gt.append({"vx":mav[0]["vx"], "vy":mav[0]["vy"], "psi":[0] * n})
	gt.append({"vx":mav[1]["vx"], "vy":mav[1]["vy"], "psi":[0] * n})

	return rssis, gt, mav, mav_rel

def create_vx(time, simcase=1):
	n = len(time)
	vx = [0]*2				# list for the x velocities for mav 0 and mav 1 * 2
	vy = [0]*2				# list for the y velocities for mav 0 and mav 1 * 2
	x = [1,1]				# set default initial position (will be overwriten for some simcase)
	y = [1,1]		

	# Switch between the 10 different simulation cases
	if simcase == 1:
		noiselevel = 0.5	# define the SD of the noise in the velocity
		# Define vx and vy, for both MAVs, following variations of sinuosidal
		vx[0] = time ** 0.5 * np.sin(0.15 * time) + np.random.normal(0, noiselevel, n)
		vx[1] = time ** 0.5 * np.sin(0.15 * time - 3) + np.random.normal(0, noiselevel, n)
		vy[0] = time ** 0.5 * np.cos(0.15 * time + 1.5) + np.random.normal(0, noiselevel, n)
		vy[1] = time ** 0.5 * np.cos(0.15 * time) + np.random.normal(0, noiselevel, n)
	elif simcase == 2:
		noiselevel = 0.2
		vx[0] = 0.9 * np.sin(0.5 * time) + np.random.normal(0, noiselevel, n)
		vx[1] = np.sin(0.5 * time - 3) + np.random.normal(0, noiselevel, n)
		vy[0] = 1.1 * np.cos(0.5 * time + 1.5) + np.random.normal(0, noiselevel, n)
		vy[1] = 0.8 * np.cos(0.5 * time) + np.random.normal(0, noiselevel, n)
	elif simcase == 3:
		noiselevel = 0.2
		vx[0] = np.sin(time) + np.random.normal(0, noiselevel / 2, n)
		vy[0] = np.sin(time) + np.random.normal(0, noiselevel / 2, n)
		vx[1] = np.cos(time) + np.random.normal(0, noiselevel, n)
		vy[1] = 0 * time + np.random.normal(0, noiselevel / 5, n)
		x = [0, 1]
		y = [1, 2]
	elif simcase == 4:
		noiselevel = 0.01
		vx[0] = [0.02] * n + np.random.normal(0, noiselevel/2, n)
		vy[0] = 0.1 * np.sin(time) + np.random.normal(0, noiselevel/2, n)
		vx[1] = [-0.02] * n + np.random.normal(0, noiselevel/2, n)
		vy[1] = 0.1 * np.sin(time) + np.random.normal(0, noiselevel/2, n)
		x = [1, 4.5]
		y = [1, 1]
	elif simcase == 5:
		noiselevel = 0.04
		vx[0] = [0.02] * n + np.random.normal(0, noiselevel / 2, n)
		vy[0] = 0.1 * np.sin(time * 0.2) + np.random.normal(0, noiselevel / 2, n)
		vx[1] = [-0.02] * n + np.random.normal(0, noiselevel / 2, n)
		vy[1] = 0.1 * np.sin(time * 0.2) + np.random.normal(0, noiselevel / 2, n)
		x = [1, 5]
		y = [1, 1]
	elif simcase == 6:
		noiselevel = 0.001
		vx[0] = np.array([0.02] * n)
		vy[0] = 0 * time + np.random.normal(0, noiselevel / 2, n)
		vx[1] = np.array([-0.02] * n)
		vy[1] = 0 * time + np.random.normal(0, noiselevel / 2, n)
		x = [1, 5]
		y = [1, 1]
	elif simcase == 7:
		noiselevel = 0.005
		vx[0] = np.array([0.02] * n)
		vy[0] = 0 * time + np.random.normal(0, noiselevel / 2, n)
		vx[1] = np.array([-0.02] * n)
		vy[1] = 0 * time + np.random.normal(0, noiselevel / 2, n)
		x = [1, 4]
		y = [0.75, 1]
	elif simcase == 8:
		noiselevel = 0.005
		vx[0] = np.array([-0.02] * n)
		vy[0] = 0 * time + np.random.normal(0, noiselevel / 2, n)
		vx[1] = np.array([0.02] * n)
		vy[1] = 0 * time + np.random.normal(0, noiselevel / 2, n)
		x = [3, 2]
		y = [0.75, 1]
	elif simcase == 9:
		noiselevel = 0.01
		sf = 2				# Speed factor
		vx[0] = np.sin(0.1*sf*time) + np.random.normal(0, noiselevel / 2, n)
		vy[0] = np.cos(0.1*sf*time) + np.random.normal(0, noiselevel / 2, n)
		vx[1] = np.sin(0.1*sf*time+np.pi/2) + np.random.normal(0, noiselevel / 2, n)
		vy[1] = np.cos(0.1*sf*time+np.pi/2) + np.random.normal(0, noiselevel / 2, n)
		x = [0, 11]
		y = [0, 10]
	elif simcase == 10:
		noiselevel = 0.01
		vx[0] = np.array([0] * n )
		vy[0] = np.array([0] * n )
		vx[1] = np.array([0] * n )
		vy[1] = np.array([0] * n )
		x = [1, 1]
		y = [1, 1.5]

	return vx, vy, x, y

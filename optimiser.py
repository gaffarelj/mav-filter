"""
@author: Jérémie Gaffarel - TU Delft - Faculty of Aerospace Engineering - BSc 2 (2018-2019)

This module is meant to optimise the Kalman Filter
"""
import filter
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import ticker
from scipy.interpolate import griddata

def optimise(time, rssi_unf, gt, gt_rel, fnumber=1, sim=True, mesh_ref=15, R_range=[0.5,15,0.01,1.5], Q_range=[0.01,1,0.01,1.5]):
	""" Optimise
	* Inputs:
		- time, rssi_unf, gt, gt_rel: inputs for the filter
		- fnumber: int: flight number or simulation case (id of the file that will be saved)
		- sim: bool: specify if a simulation is being optimised or a real flight
		- mesh_ref: mesh refinement for finding the optimum SD combination
		- R/Q_range: range of SDs to be tested
	"""
	# Mesh refinement factor has to be at least 2
	if mesh_ref > 1:
		m = mesh_ref
	else:
		m = 2
	# Generate lists of all combinations of SDs to be tested for a and b
	SRa = np.linspace(R_range[0], R_range[1], m)
	SRb = np.linspace(R_range[2], R_range[3], m)
	# Get best SDs for R, and the error that is obtained with those best SDs
	best_r, error = opti(SRa, SRb, [0.1, 0.5], m, True, time, rssi_unf, gt, gt_rel, fnumber, sim)
	print("Optimised SD for R matrix are:", best_r)
	print("Giving an average range error of:", error, "[m]")
	# Same process for Q
	SQa = np.linspace(Q_range[0], Q_range[1], m)
	SQb = np.linspace(Q_range[2], Q_range[3], m)
	best_q, error = opti(SQa, SQb, [5, 0.2], m, False, time, rssi_unf, gt, gt_rel, fnumber, sim)
	print("Optimised SD for Q matrix are:", best_q)
	print("Giving an average range error of:", error, "[m]")
	return best_r, best_q, error

def run_filter(s_r, s_q, time, rssi_unf, gt, gt_rel):
	""" Run filter
	* Run the filter and outputs the average error in range
	"""
	n = len(time)
	x_filter, d_g, d_f, d_unf, b_g, b_f, vel_g, vel_f, rssi_f = filter.kalman_filter(time, rssi_unf, gt, gt_rel, s_r=s_r, s_q=s_q, optimisation=True)
	error = (sum(np.fabs(d_f - d_g))) / n
	return error

def opti(Sa, Sb, other, m, opti_r, time, rssi_unf, gt, gt_rel, fnumber, sim):
	""" Optimisation
	* Inputs:
		- Sa: list with all of the SD to be tested for the first SD
		- Sb: same for the second SD
		- other: SDs values for the matrix that is not being optimised
		- m: mesh refinement factor
		- opti_r: if True, R is being optimised, otherwise Q is
		- time, rssi_unf, gt, gt_rel: intputs for the Kalman Filter
		- fnumber: flight or simulation number
		- sim: if True, simulation is being tested
	"""
	z = []
	x = []
	y = []
	opti_pr = "R optimisation:" if opti_r else "Q optimisation:"
	# Try every combination of SDs in the input ranges
	for i in range(m):
		print(opti_pr, round(100 * i / m, 2), "%")
		for j in range(m):
			if opti_r:
				error = run_filter([Sa[i], Sb[j]], other, time, rssi_unf, gt, gt_rel)
			else:
				error = run_filter(other, [Sa[i], Sb[j]], time, rssi_unf, gt, gt_rel)
			# Save every combinations of SDs and the error associated with them
			x.append(Sa[i])
			y.append(Sb[j])
			z.append(error)
		
	# Make a cubic interpolation of the computed values on a 100x100 mesh
	grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
	grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')

	# Plot the optimisation output
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.view_init(azim=32, elev=44)
	surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap=plt.cm.Spectral)
	ax.set_xlabel("SD-a")
	ax.set_ylabel("SD-b")
	ax.set_zlabel("Average range error [m]")
	situation = "sims" if sim else "flights"
	opti_m = "r" if opti_r else "q"
	cb = plt.colorbar(surf, shrink=0.5)
	tick_locator = ticker.MaxNLocator(nbins=5)
	cb.locator = tick_locator
	cb.update_ticks()
	plt.savefig("plots/optimisation/" + situation + "/" + str(fnumber) + "_" + opti_m + ".pdf")
	plt.show()
	plt.close()

	# Find the x and y coordinates of the minimum z-value in the mesh, that is the minimum error
	xmin, ymin = np.unravel_index(np.argmin(np.fabs(grid_z)), grid_z.shape)
	best = [grid_x[xmin,ymin], grid_y[xmin,ymin]]
	# Compute the error again for these best SDs
	if opti_r:
		error = run_filter(best, other, time, rssi_unf, gt, gt_rel)
	else:
		error = run_filter(other, best, time, rssi_unf, gt, gt_rel)
	return best, error

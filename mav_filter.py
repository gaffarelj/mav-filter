c = "@authors: Jérémie Gaffarel & Rudi Smits - TU Delft - Faculty of Aerospace Engineering - BSc 2 (2018-2019) \n\n\t\t*** Kalman Filter of noise on measurements of relative position between 2 MAVs ***	\n\n\n                             ::\n                        .:+yh+\n                   /ohmNNmo-\n                +hNNNNNm+`\n              :mNNdsmNNd-.:s:\n             -mms-:mNNNNNNNs\n             /m.  -ss/yNNm/\n             :o     `sNmo`\n              `    :s+-\n         mmmmmmmmmmmmmmmd .mmmmm`     .mmmmm`   -mmmmmmmmhy+.                  ymm     sdNNN+\n         mmmmmmNNNNmmmmmd .NNNNN`     .NNNNN`   :NNs:::::+dNNs                 hNN    /NNo```  +hm\n              hNNNN/      .NNNNN`     .NNNNN`   :NN+       +NNh     ./+++:     hNN  -:oNNs:: ::hNN/:-\n              hNNNN/      .NNNNN`     .NNNNN`   :NN+        yNN:  :dNhooymms   hNN  sydNNdyy symNNyys\n              hNNNN/      .NNNNN`     .NNNNN`   :NN+        oNN+ -mN+    -NNs  hNN    /NN/     yNN`\n              hNNNN/      .NNNNN`     .NNNNN`   :NN+        yNN: sNNddddddNNh  hNN    /NN/     yNN`\n              hNNNN/      `dNNNNs`   `sNNNNd    :NN+       +mNh  oNN:.....--.  hNN    /NN/     yNN`\n              hNNNN/       -dNNNNmhhhmNNNNd-    :NNs:::::odNNs`  .mNy.   oNN+  hNN    /NN/     yNN`\n              ymmmm/        `/ydmNNNNNmhs:      -mmmmmmmmhy+/     .smNddmmy/   ymm    :mm/     /dNNd+\n                                 `````                               ````                        ```\n\n"
print(c)

import data
import filter
import simulation
import optimiser
import numpy as np
from matplotlib import pyplot as plt

def plot(time, x_f, x_g, x_u, ylabel, fname, simchoice, it0=0, ite=-1):
	if simchoice == 1:
		fname = "plots/sims/" + fname + ".pdf"
	else:
		fname = "plots/flights/" + fname + ".pdf"
	plt.plot(time[it0:ite], x_f[it0:ite], linewidth=0.75, label="Filter", marker="1")
	plt.plot(time[it0:ite], x_g[it0:ite], linewidth=1, label="Ground truth", linestyle="-")
	if len(x_u) > 1:
		plt.scatter(time[it0:ite], x_u[it0:ite], label="Unfiltered", color=(0.078, 0.568, 0, 0.35), marker="2")
	plt.ylabel(ylabel)
	plt.xlabel("Time [s]")
	plt.grid()
	plt.legend(loc="upper right")
	plt.savefig(fname)
	plt.show()
	plt.close()

def save_csv(data, line):
	fname = "results.csv"
	f = open(fname, "r")
	lines = f.readlines()
	f.close()
	while len(lines) < line + 1:
		lines.append("\n")
	lines[line] = (',').join(data) + "\n" 
	with open(fname, 'w') as file:
		file.writelines(lines)

def get_choice(s, e, message):
	choice = -1
	while not (choice in np.arange(s, e + 1)):
		try:
			choice = int(input(message))
		except:
			pass
	return choice

def get_correlation(noise, para, name):
	corre = np.corrcoef(para, noise)[0][-1]
	print("\t", name, corre)
	return corre

sim_choice = get_choice(0, 1, "Enter 0 for flight data or 1 for the simulation: ")
if sim_choice == 1:
	sim_case = get_choice(1, 10, "Simulation to run [1-10]: ")
	time = np.arange(0, 120, 0.2)
	print("Running the simulation...")
	rssi_unf, gt, onb_avg, gt_rel = simulation.simulate(time, sim_case)
else:
	f_choice = get_choice(1, 4, "Flight to filter [1-4]: ")
	print("Opening Data file...")
	sim_case = f_choice
	gt, onb, onb_avg, gt_rel, time = data.get_data(f_choice)
	rssi_unf = np.array(onb_avg["rssi"])

#t0 = 50
#te = 110
#it0 = np.where(time >= t0)[0][0]
#ite = np.where(time >= te)[0][0]
it0 = 1
ite = -1

opti_choice = get_choice(0, 1, "Enter 0 for direct filtering or 1 for optimisation: ")
if opti_choice == 0:
	s_r, s_q = [11.44823232, 0.160505051], [0.0275, 0.183080808]
else:
	s_r, s_q, error = optimiser.optimise(time, rssi_unf, gt, gt_rel, fnumber=sim_case, sim=(sim_choice == 1))

# x_filter = ["x_rel", "y_rel", "vx", "vy", "vx_other", "vy_other", "psi", "psi_other", "z_rel"]
print("Running the Kalman Filter...")
x_filter, d_g, d_f, d_unf, b_g, b_f, vel_g, vel_f, rssi_f = filter.kalman_filter(time, rssi_unf, gt, gt_rel, s_r=s_r, s_q=s_q)

# Getting correlation between noise and different parameters
noise = d_unf - d_f
print("Correlations with noise:")
corre_range = get_correlation(noise, d_f, "range")
corre_vel = get_correlation(noise, vel_f, "velocity")
corre_bearing = get_correlation(noise, b_f, "bearing")
if sim_choice == 0 and opti_choice == 1:
	save_csv(["Flight", "R a", "R b", "Q a", "Q b", "average range error [m]", "correlations:", "range", "velocity", "bearing"], 0)
	save_csv([str(i) for i in [sim_case, s_r[0], s_r[1], s_q[0], s_q[1], error, "", corre_range, corre_vel, corre_bearing]], sim_case)
elif sim_choice == 1 and opti_choice == 1:
	save_csv(["Simulation", "R a", "R b", "Q a", "Q b", "average range error [m]", "correlations:", "range", "velocity", "bearing"], 10)
	save_csv([str(i) for i in [sim_case, s_r[0], s_r[1], s_q[0], s_q[1], error, "", corre_range, corre_vel, corre_bearing]], sim_case + 10)

# PLOTS of x-relative, y-relative, filter vs no filter and relative bearing
if sim_choice == 1:
	plt.plot(onb_avg[0]["x"][it0:ite], onb_avg[0]["y"][it0:ite], label="MAV 1", linewidth=1, marker="1")
	plt.plot(onb_avg[1]["x"][it0:ite], onb_avg[1]["y"][it0:ite], label="MAV 2", linewidth=1, marker="2")
	plt.ylabel("Y position [m]")
	plt.xlabel("X position [m]")
	plt.grid()
	plt.legend(loc="upper right")
	plt.show()
	plt.savefig("plots/sims/" + str(sim_case) + "_path.pdf")
	plt.close()
plot(time, vel_f, vel_g, [0], "Relative velocity [m/s]", str(sim_case) + "_rel_velocity", sim_choice, it0=it0, ite=ite)
plot(time, d_f, d_g, d_unf, "Relative distance [m]", str(sim_case) + "_rel_range", sim_choice, it0=it0, ite=ite)
plot(time, b_f, b_g, [0], "Relative bearing [rad]", str(sim_case) + "_rel_bearing", sim_choice, it0=it0, ite=ite)

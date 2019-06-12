# MAV Filter
This code is meant to open Micro Air Vehicle (MAV) flights data, simulate flights, filter the noise in the measurements, and find correlation between that noise and other parameters.

## Prerequisites
This code requires Python 3 to be installed on your machine. It can be downloaded [here](https://www.python.org/downloads).

The Numpy, Matplotlib and SciPy modules are also required to be installed. This can be done by running the following commands if Pip was installed with Python:
```
pip install numpy
pip install matplotlib
pip install scipy
```

## Modules
The main file is [mav_filter.py](mav_filter.py). Alongside it, 4 different modules were written.

### Data
[data.py](data.py) is a module meant to open the MAV flights meausurements.

As an input, it takes an integer ranging from 1 to 4, to open the data of flights 1 to 4.

It then outputs the followings: 
* onb[_i_]: list of parameters measured on-board the MAV:
    * id: id of MAV _i_
    * id_other: id of the other MAV, _j_
    * rssi: RSSI of the signal between _i_ and _j_
    * x/y/z_rel: relative position of MAV _i_ to _j_
    * vx/vy: relative velocity of MAV _i_ to _j_
    * vx/vy_other: relative velocity of the other MAV _j_ to _i_
    * psi: relative angle wrt the magnetic north of MAV _i_
    * psi_other: relative angle wrt the magnetic north of MAV _j_
* gt[_i_]: list of parameters measured by an [OptiTrack motion capture system](https://optitrack.com/):
    * x/y/z: absolute positions of MAV _i_
    * vx/vy: absolute velocity of MAV _i_
    * psi: relative angle wrt the magnetic north of MAV _i_
* onb_avg: list containing the averaged measurements of the two MAVs:
    * rssi
    * x/y/z_rel
* gt_rel: list containing the ground-true relative postisions of the two MAVs
    * x/y/z
    
Here is how this module can be called:
```
gt, onb, onb_avg, gt_rel, time = data.get_data(3)   # For data of flight 3
```

If one want to get the data files that were used in this project, a request can be made to the authors of the article from which this work originated. (cf [Acknowledgments](#acknowledgments))

### Simulation
[simulation.py](simulation.py) is a module meant to simulate MAV flights.

As inputs, it takes a time list for which to run the simulation, and an integer ranging from 1 to 10, to specify the simulation case.

It then outputs the followings:
* rssis: Simulated rssi between the two MAVs, with noise
* gt: Ground truth simulated positions of the two MAVs
* mav: Simulated positions of the two MAVs
* mav_re: Relative simulated positions of the two MAVs

Here is how this module can be called:
```
time = np.arange(0, 120, 0.2)
rssi, gt, onb_avg, gt_rel = simulation.simulate(time, 5)  # For simulation case 5
```

### Filter
[filter.py](filter.py) is a module meant to filter the noisy measurements of MAVs, by the mean of an Extended Kalman FIlter (EKF).

It takes the following inputs:
* time: list containing the time at wich each measurement was taken
* rssi: noisy RSSI between the two MAVs
* gt: ground-true measurements
* gt_rel: ground-true relative measurements
* Optional:
    * s_q = [s_a, s_b]: s_a and s_b are here the Standard Deviations (SD) used in the Q matrix of the filter
    * s_r = [s_a, s_b]: same for R matrix
    * optimisation: boolean: if True, no additional noise is added to the ground-true measurements
    
It then outputs the followings:
* x_filter: dict containing x/y/z_rel, vx/vy, vx/vy_other, psi, psi_other
* d/b/vel/rssi_f: filtered range/bearing/absolute velocity/rssi
* d/b/vel_g: ground-true range/baering/velocity
* d_unf: unfiltered range

Here is how this module can be called:
```
s_r, s_q = [5, 0.2], [0.1, 0.5]
x_filter, d_g, d_f, d_unf, b_g, b_f, vel_g, vel_f, rssi_f = filter.kalman_filter(time, rssi_unf, gt, gt_rel, s_r=s_r, s_q=s_q)
```

### Optimiser
[optimiser.py](optimiser.py) is a module meant to find the SD combination that shall be used in the EKF to get a minimum error.

It takes the following inputs:
* time
* rssi_unf: unfiltered RSSI
* gt
* gt_rel
* Optional:
    * fnumber: int: flight number or simulation case (id of the file that will be saved)
    * sim: bool: specify if a simulation is being optimised or a real flight
    * mesh_ref: mesh refinement for finding the optimum SD combination
    
It then outputs the followings:
* s_r: list containing the optimum SD for the R matrix
* s_q: same for the Q matrix
* error: average range error in [m] obtained whilst using the optimum SD combination

Here is how this module can be called:
```
r, q = [0.5,15,0.01,1.5], [0.01,1,0.01,1.5]
s_r, s_q, error = optimiser.optimise(time, rssi_unf, gt, gt_rel, mesh_ref=15, R_range=r, Q_range=q)
```

## Results
[results.csv](results.csv) contains the results of the filter and optimisation made to run for the 4 different flights and the 10 simulations.

The [plots](plots) folder contains plots of all of the flights, simulations and optimisations.

### Correlations
The correlations between the noise in the range between the MAVs and different parameters are reported in the [results](results.csv) file. 
The following table report the average correlations of the 4 flights.

| Parameter | Average correlation | Interpretation |
| --- | --- | --- |
| Range | -0.766270552 | The range between the MAVs has a strong influence on the noise in the range |
| Velocity | 0.017549867 | The relative velocity of the MAVs has no influence on the noise in the range (velocities < 2 [m/s]) |
| Bearing | -0.210140569 | The bearing of the MAVs has a slight influence on the noise in the range |

All of the correlations in the simulation are close to 0. That is because the noise is randomly added to the range in the simulation. This explains why no parameters influence the noise in the range in the simulations.

### Optimisation
The SD giving the lowest average error in range in the filter are also reported in the [results](results.csv) file. 

Averaging the results, the best SD for the R matrix are: 11.45 for the first SD, and 0.16 for the 3 last ones.

The best SD for the Q matrix are: 0.0275 for the first SD, and 0.183 for the 3 last ones.

These results were computed using a mesh refinement parameter of 15, and ranges of [0.5, 15] and [0.01, 1.5] for the SD of the R matrix, and ranges [0.01, 1] and [0.01, 1.5] for the SD of the Q matrix.

### Plots
The flights bearing, velocity and range plots can be found in the [flights](plots/flights) folder.

The same plots for the simulations, with the path of the MAVs in addition can be found in the [sims](plots/sims) folder.

Finally, plots of the SD vs the average error in range were produced during the optimisation process, and can be found in the [optimisation](plots/optimisation) folder.

## License
This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details.

## Authors
This code was produced as a student work for the Test, Analysis and Simulation (AE2223-I) course given to the second-year bachelors in the faculty of Aerospace Engineering, at the Delft University of Technology.
It was mostly written by the following two students:
* Jérémie Gaffarel
* Rudi Smits

## Acknowledgments
The [filter](filter.py) module of this project was inspired by the one from [Paparazzi](https://github.com/coppolam/paparazzi/tree/paper_AURO_ARdroneexperiments_optitrack/sw/airborne/modules/relativeavoidancefilter).

All of the work done in this project originated from the following article:

Coppola, M., McGuire, K.N., Scheper, K.Y.W. et al. Auton Robot (2018) 42: 1787. [https://doi.org/10.1007/s10514-018-9760-3](https://doi.org/10.1007/s10514-018-9760-3)

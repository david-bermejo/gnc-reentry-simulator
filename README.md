# Mars Reentry Simulator
High fidelity simulator for a lifting vehicle during Mars reentry, designed for my Final Master Thesis: ***Design and implementation of a lifting body
GNC for human Martian Entry***. It uses the HORUS-2B aerodynamic database, and provides a 6 degrees of freedom (6-DoF) description of the spacecraft motion. The
GNC subsystem tracks a pre-optimized trajectory considering vehicle and human safety constraints, with guidance provided by a Linear Quadratic Regulator (LQR)
and a bank reversal switching law. Navigation is done using an Extended Kalman filter, fusing Inertial Measurement Unit (IMU) increments and pseudorange measurements
from Martian surface beacons. The control module is operated via two longitudinal and lateral LQR controllers. Control commands are issued to an actuator allocation
function, operated by five aerosurfaces and a 14-thruster Reaction Control System. Monte Carlo simulations validated the system’s performance and adaptability to
stochastic uncertainties.

---

![Simulator design overview](/doc/sim.png)

## Initialization
The simulator environment is initialized by running the script:
```
load_simulation.m
```
The script constructs ```simdb``` structure inside MATLAB workspace, which contains all neccessary substructures given by each script located inside ```input``` folder. Once done, the script
automatically opens the simulator model in Simulink, which is located in ```models\simulator.slx```.

## Monte Carlo Simulation Campaign
 A Monte Carlo (MC) simulation campaign is performed to validate the designed Guidance, Navigation and Control subsystem. The MC method is an algorithm used to predict a
set of outcomes based on an estimated range of values for given input variables. In the context of this thesis, uncertainty is introduced to the initial state of the winged
vehicle and the atmospheric density profile, and simulation results are post-processed in bulk to verify the performance and robustness of the GNC. Atmospheric uncertainty
is modeled as an stochastic Ornstein-Uhlenbeck process in discrete time, while individual variables of spacecraft’s initial state are perturbed with zero-mean Gaussian distributions.
A MATLAB script has been programmed, which allows to execute the simulations in parallel using MATLAB's *Parallel Computing Toolbox*. The simulation campaign is simply run as:
```
simulation_campaign.m
```
A default of 6 parallel workers and 84 batches are selected. The simulation campaign has been run on a computer with a Ryzen 5 5600x processor and 32GB of RAM
memory, where the runtime was 80 min.

## Post-Processing
Finally, a post-processing script allows to obtain the neccessary figures to validate the designed GNC subsystem, which can be simply executed as:
```
postpro.m
```

## Conclusions
Finally, it is shown that Monte Carlo simulations validated the system’s performance and adaptability to stochastic uncertainties. Although wind dynamics were not considered, the GNC
architecture exhibits promising trajectory tracking and robustness against disturbances, suggesting a solid foundation for further research into more advanced GNC techniques for future Martian
missions.
# cyanobacteria-filament-fragmentation
 Code and data for the paper "Phenotypic plasticity, life cycles, and the evolutionary transition to multicellularity"

# 01 Simulation code

Contains the code for each model. The provided code is intended to run on a cluster managed by slurm (see run.sh).  The output of a single run (saved in /Results) is an optimisation history - parameters of the model and corresponding deviation from experimental data.

The calculation "Control no growth" can be run on a personal machine, as it computes only the deviation for a control model.


# 02 Simulation results
Contains the results of model optimisation. The results files are stored in /XYZ_model/Results/.

## Data format
The results files record the history of optimisation: the first line corresponds to the initial parameters values, each next line corresponds to the found values for which the deviation between the model and the experimental data was lower than in the previous line.

In each line, the first 10 elements of each line represent model parameters. Some of these are not used and labelled as empty_1, empty_2 and so on - ignore them. Some parameters are not optimised, as they are scaling factors. See the model code for details. 

The last element "Deviation" represents the difference between experimental and simulated population dynamics.

# 03 Data processing
Contains the code producing figures presented in the paper.
Note that the plotted value of deviation between experiment and simulations is scaled to the control model value. That control value is computed by /01 Simulation code/Control no growth/Control_no_growth.py
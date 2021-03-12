# How to plot

Before to run the user need to edit the fist three lines in example_plots.R to add: 
  - dir: Simulation directory
  - N: Pupularion size 
  - M: Number of simulations in the folder

# Output 

The script example_plots.R generate three figures:
  - rank_degree.png: The number of risky contacts individuals have, as a function of the order in which they are
infected for the M simulations
  - infect_curve.png: Mean of the M infected curves with the standard deviation
  - herd_immunity.png: Histogram of the necessary simulations to have M successfully simulations

The generated figures will be located in *dir* folder. 

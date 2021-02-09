# The-role-of-connectivity-on-COVID-19-preventive-approaches

This repository contains simulation scripts used in the article "The role of connectivity on COVID-19 preventive approaches" (in prep.). 
They allow to simulate the spread of an epidemic in a population where the contact structure is given by a random graph (Erdős-Rényi or Power-Law degree distribution graph). Vertices of the graph represent individuals of the populations and edges represent risky contacts during the time they are infectious. On top of this structure, the epidemics is propagated under a classical SIRD model. The details of of the parameters and their interpretation can be found in the Supplementary Materials of the paper or at [add link to preprint].


# How to run the simulations
All scripts run on an Unix-based system. They involve parallelized computations, which will not (directly) work under Windows operating system.

The script control.R allows to simulate the spread of an epidemic without any human intervention. It can be run by modifying the following command line: 

Rscript control.R -m 30 -c 30 -n 20000 -e 44 -d 0.07, 0.005 -i 0.05 -t 300 -x 50 -q 0.17   

where : 
 - m is the number of independent simulations you want to run
 - c is the number of cores to be used to run the script
 - n is the population size
 - e is the expected degree of the graph (average number of contacts per individual)
 - d are the death probabilities
 - lambda is a parameter related to the exponent of the power law degree distribution (see Qiao et al. (2018)). To run a simulation under an Erdős-Rényi graph set lambda to -1. 
 - i is the probability that an infected individual infects on of its neighbors
 - t is the maximum time for the simulations (
 - filename is a prefix for the name of the output folders. 

The script lockdown.R allows to simulate the spread of an epidemic when a lockdown takes place. It can be run by  modifying  the following command:

Rscript lockdown.R -

The script lockdown.R allows to simulate the spread of an epidemic and vaccinate part of the population at some moment. It can be run by  modifying  the following command:

Rscript vaccination.R -

All the other scripts are dependencies of these three scripts. 

# Structure of the output

The scripts will create different folders (all contained in the main folder "filename") that contain the output : 

-  G : contains m files called 1.txt, 2.txt ... Each one contains the graph used in one simulation. It is a list of all the edges in the graph, e.g. if the first line is "17 28" it means that individual 17 is connected to individual 28 in the graph. 


# Plots 


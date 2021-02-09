# The-role-of-connectivity-on-COVID-19-preventive-approaches

This repository contains simulation scripts used in the article "The role of connectivity on COVID-19 preventive approaches" (in prep.). 
They allow to simulate the spread of an epidemic in a population where the contact structure is given by a random graph (Erdős-Rényi or Power-Law degree distribution graph). Vertices of the graph represent individuals of the populations and edges represent risky contacts during the time they are infectious. On top of this structure, the epidemics is propagated under a classical SIRD model. The details of of the parameters and their interpretation can be found in the Supplementary Materials of the paper or at [add link to preprint].


# How to run the simulations
All scripts run on an Unix-based system. They involve parallelized computations, which will not (directly) work under Windows operating system.

The script control.R allows to simulate the spread of an epidemic without any human intervention. It can be run by modifying the following command line: 

Rscript control.R -rep, n, lambda, N, pD, pI, filename   

where : 
 - rep is the number of independent simulations you want to run
 - N is the population size
 - lambda is a parameter related to the exponent of the power law degree distribution (see Qiao et al. (2018)). To run a simulation under an Erdős-Rényi graph set lambda to -1. 
 - pI is the probability that an infected individual infects on of its neighbors
 - filename is a prefix for the name of the output folders. 

The script lockdown.R allows to simulate the spread of an epidemic when a lockdown takes place. It can be run by  modifying  the following command:

Rscript lockdown.R -

The script lockdown.R allows to simulate the spread of an epidemic and vaccinate part of the population at some moment. It can be run by  modifying  the following command:

Rscript vaccination.R -


# Structure of the output

# Plots 


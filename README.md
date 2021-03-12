# The role of connectivity on COVID-19 preventive approaches

This repository contains simulation scripts used in the article "The role of connectivity on COVID-19 preventive approaches" (in prep.). 
They allow to simulate the spread of an epidemic in a population where the contact structure is given by a random graph (Erdős-Rényi or Power-Law degree distribution graph). Vertices of the graph represent individuals of the populations and edges represent risky contacts during the time they are infectious. On top of this structure, the epidemics is propagated under a classical SIRD model. The details of of the parameters and their interpretation can be found in the Supplementary Materials of the paper or at [add link to preprint].


# How to run the simulations
All scripts run on an Unix-based system. They can involve parallelized computations, which will not (directly) work under Windows operating system. Windows users have to be sure that option -c 0 is selected (see below). 
To run the scripts, the user needs to clone this repository and to have installed R as well as the following libraries : Rcpp, igraph, parallel and getopt.

Before running any simulation, the user needs to change the first line in control.R, vaccination.R and lockdown.R to add the path to the directory that contains all the scripts.


The script control.R allows to simulate the spread of an epidemic without any human intervention. 
It can be run by modifying the following command line: 

Rscript control.R -m 30 -c 30 -n 20000 -e 44 -r 0.1 -d1 0.005 -d2 0.07 -q 0.17  -i 0.05 -t 300 -x 50 -l 0.5 -s foldername  

where : 
 - m is the number of independent simulations you want to run,
 - c is the number of cores to be used to run the script. If this value is set to 0 the script will not use parallelized computations (and will be suitable for Windows users),
 - n is the population size,
 - e is the expected degree of the graph (average number of contacts per individual),
 - r is the recovery rate (i.e. 1 over the average time an individual remains infected),
 - d1 d2 are the death probabilities (first value is for the highly connected individuals, second value for the lowly connected individuals),
 - q is the fraction of "vulnerable" individuals (i.e. less  connected and with a higher death rate),
 - i is the probability that an infected individual infects on of its neighbors,
 - t is the maximum time for the simulations (if -t 0 the simulations won't stop until the epidemics die out),
 - x is the minimal number of infected individuals that has to be reached to consider a simulation as successful (to avoid considering epidemics that die out very quickly without reaching an exponential phase), 
 - l is a parameter related to the exponent of the power law degree distribution (see Qiao et al. (2018)). To run a simulation under an Erdős-Rényi graph set lambda to -1,
 - s is a prefix for the name of the output folders. 



The script lockdown.R allows to simulate the spread of an epidemic when a lockdown takes place. It can be run by  modifying  the following command:

Rscript lockdown.R -m 30 -c 30 -n 20000 -e 44 -r 0.1 --d1 0.005 --d2 0.07 -q 0.17  -i 0.05 -t 300 -x 50 -l 0.5 -s foldername --oLD 1 -v 0.1 --pLD 0.5 --mLD 10 --dLD 30

where the additional parameters are : 
 - oLD is strategy of the Lockdown (1 or 2, see Supplementary Methods),
 - v is the "intervention time", i.e. cumulative fraction of the population that has to be infected  before the intervention is made,
 - pLD is the probability to remove each vertex in strategy 1, 
 - mLD is the maximum number of edges per vertex in strategy 2,
 - dLD is the duration of the lockdown in days.
 
 

The script vaccination.R allows to simulate the spread of an epidemic and vaccinate part of the population at some moment. It can be run by  modifying  the following command:

Rscript vaccination.R  -m 30 -c 30 -n 20000 -e 44 -r 0.1  --d1 0.005 --d2 0.07 -q 0.17 -i 0.05 -t 300 -x 50 -l 0.5 -s foldername -b 0.25 -v 0.3 -w 1 --pE 0.9 -o 1

where the additional parameters are: 
 - b is the fraction of the population that we can vaccinate with the available doses,
 - v is the  "vaccination time", i.e. cumulative fraction of the population that has to be infected  before the intervention is made,
 - w indicates which categories of individuals you can vaccinate (1: only susceptibles S, 2: susceptible, recovered or infected SIR),
 - pE is the efficacy of the vaccine,
 - o is the strategy,  1: Uniform, 2: Neighbor, 3: Among most connected, 4: Most connected, 5: Among least connected, 6: Least connected.
 


All the other scripts in this repository are dependencies of these three scripts. 

# Structure of the output

The scripts will create different folders (all contained in the main folder "foldername") that contain the output : 

G : contains m files called 1.txt, 2.txt ... Each one contains the graph used in one simulation. It is a list of all the edges in the graph, e.g. if the first line is "17 28" it means that individual 17 is connected to individual 28 in the graph.

DG:  contains m files called 1.txt, 2.txt ... Each file corresponds to one simulation and has one line per individual. The first value is its degree in the graph (number of contacts). The second and third values recall the graph structure and the population size.

SIR : contains m files called 1.txt, 2.txt ... Each file corresponds to one simulation. Each line corresponds to one day. Columns 1 to 4 are the number of Susceptible, Infected, Recovered and Dead at the end of this day. 

DF : contains m files called 1.txt, 2.txt ... Each file corresponds to one simulation and has one line per individual. The first column is the label of the individual. The second column is the status of the individual (S, I, R or D) at the end of the simulation. The third column is the time of infection. The fourth one is the recovery time (or death time). The fifth is the label of the individual who infected this individual. The sixth is the number of individuals infected by this individual. The seventh is the last time he infected someone else. The last column corresponds to its degree in the grap (number of contacts).

HI : contains m files called 1.txt, 2.txt ...  These are single line files, each one corresponding to one simulation. The first value is the cumulative number of infected individuals during the epidemic. The second and third values recall the graph structure and the population size.

NW : contains m files called 1.txt, 2.txt ... Each one contains a Newick format tree representing the transmission tree.






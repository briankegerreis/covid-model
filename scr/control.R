# Set COVID code directory 

source("base_sim.R")
source("base_helpers.R")
Rcpp::sourceCpp("Block.cpp")

library(getopt)
spec=matrix(c("num_sim",'m',1,"integer","number of simulations to run",
              "num_cores",'c',1,"integer","number of cores to use",
              "pop_size",'n',1,"integer","population size",
              "d_connect",'e',1,"integer","expected degree of the graph",
              "lamb",'l',1,"double","tuning parameter for power law degree distribution. Set to -1 for Erdos-Renyi graph",
              "death_high",'1',1,"double","death probability for highly connected individuals",
              "death_low",'2',1,"double","death probability for less connected individuals",
              "recovery_rate",'r',1,"double","recovery rate (1/duration of infection)",
              "p_infect",'i',1,"double","probability of infecting a neighbor",
              "t_max",'t',1,"integer","maximum length of simulations. Set to 0 for no limit (epidemic must die out)",
              "f_vulnerable",'f',1,"double","fraction of vulnerable individuals, defined as less connected with higher death rate",
              "min_cases","x",1,"integer","minimum number of infections for a successful simulation",
              "pre","s",1,'character',"output folder prefix")
            ,byrow=TRUE,ncol=5)

opt=getopt(spec)
popSize=opt$pop_size
d_connect=opt$d_connect
lamb = opt$lamb
p_death=c(opt$death_high,opt$death_low)
recovery_rate=opt$recovery_rate
p_infect=opt$p_infect
numCores=opt$num_cores
numSim=opt$num_sim 
pre=opt$pre
if(opt$t_max == 0){
  t_max = Inf
}else{
  t_max = opt$t_max
}
f_vulnerable=opt$f_vulnerable
min_cases=opt$min_cases

invisible(COVID_control(num_sim, num_cores,
                        pop_size, d_connect, lamp,
                        p_death, recovery_rate, p_infect,
                        t_max, f_vulnerable, min_cases, pre))
# invisible(COVID_control(M,mc,N,d,lamb,pD,rR,pI,tStop,q,nIs,pre))










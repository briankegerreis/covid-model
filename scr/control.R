# Set COVID code directory 

codedir <- '/N/project/Covid/scr/'
source(paste0(codedir,'COVID_control.R'))
source(paste0(codedir,'covidcontrol.R'))
Rcpp::sourceCpp(paste0(codedir,"Block.cpp"))

library(getopt)
spec=matrix(c("numSim",'m',1,"integer","number of simulations to run",
              "mc",'c',1,"integer","number of cores to use",
              "popSize",'n',1,"integer","population size",
              "dConnect",'e',1,"integer","expected degree of the graph",
              "lamb",'l',1,"double","tuning parameter for power law degree distribution. Set to -1 for Erdos-Renyi graph",
              "d1",'1',1,"double","death probability for highly connected individuals",
              "d2",'2',1,"double","death probability for less connected individuals",
              "rR",'r',1,"double","recovery rate (1/duration of infection)",
              "pI",'i',1,"double","probability of infecting a neighbor",
              "tStop",'t',1,"double","maximum length of simulations. Set to 0 for no limit (epidemic must die out)",
              "q",'q',1,"double","fraction of vulnerable individuals, defined as less connected with higher death rate",
              "x","x",1,"double","minimum number of infections for a successful simulation",
              "pre","s",1,'character',"output folder prefix")
            ,byrow=TRUE,ncol=5)

opt=getopt(spec)
N=opt$popSize
d=opt$dConnect
lamb = opt$lamb
pD=c(opt$d1,opt$d2)
rR=opt$rR
pI=opt$pI
mc=opt$mc 
M=opt$numSim 
pre=opt$pre
if(opt$tStop == 0){
  tStop = Inf
}else{
  tStop = opt$tStop
}
q=opt$q
nIs=opt$x


invisible(COVID_control(M,mc,N,d,lamb,pD,rR,pI,tStop,q,nIs,pre))

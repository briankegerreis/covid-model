codedir <- '/N/project/Covid/scr/'
source(paste0(codedir,'COVID_vaccination.R'))
source(paste0(codedir,'covidcontrol.R'))
Rcpp::sourceCpp(paste0(codedir,"Block.cpp"))

library(getopt)
spec=matrix(c("numSim",'m',1,"integer",
              "mc",'c',1,"integer",
              "popSize",'n',1,"integer",
              "dConnect",'e',1,"integer",
              "lamb",'l',1,"double",
              "d1",'1',1,"double",
              "d2",'2',1,"double",
              "rR",'r',1,"double",
              "pI",'i',1,"double",
              "tStop",'t',1,"double",
              "q",'q',1,"double",
              "x","x",1,"double",
              "pre","s",1,'character',
              "tV","v",1,"double",
              "D","b",1,"double",
              "who","w",1,"integer",
              "oV","o",1,"double",
              "pE","3",1,"double")
            ,byrow=TRUE,ncol=4)

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
tV=N*opt$tV
D=N*opt$D
oV=opt$oV
who=if(opt$who == 1){c('S')}else{c('S','I','R')}
pE = opt$pE

invisible(COVID_vacc(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,tV,D,oV,who))

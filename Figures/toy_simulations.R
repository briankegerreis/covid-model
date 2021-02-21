codedir <- '/N/project/Covid/scr/'

Rcpp::sourceCpp(paste0(codedir,"Block.cpp"))
source(paste0(codedir,'covidcontrol.R'))
source(paste0(codedir,'COVID_lockdown.R'))
source(paste0(codedir,'COVID_vaccination.R'))
source(paste0(codedir,'COVID_control.R'))

N=100
d=6
pD=c(0.1,0.1)
rR=0.1
pI=0.5
mc=0
M=30
tStop = Inf
q=0.17
nIs=10

# Erdos Control
invisible(COVID_control(M,mc,N,d,-1,pD,rR,pI,tStop,q,nIs,'ER'))
# PLD Control
invisible(COVID_control(M,mc,N,d,0.5,pD,rR,pI,tStop,q,nIs,'PLD'))

### Lockdown
library(igraph)
# Lockdwon paramets

tV = 0.1*N
dLD = 14
pLD = 0.6
mLD = 4

# Lockdown type 1 
oLD = 1
# Erdos
invisible(COVID_ld(M,mc,N,d,-1,pD,rR,pI,tStop,q,nIs,'ER',oLD,tV,pLD,mLD,dLD))
# PLD
invisible(COVID_ld(M,mc,N,d,0.5,pD,rR,pI,tStop,q,nIs,'PLD',oLD,tV,pLD,mLD,dLD))

# Lockdown type 2
oLD = 2
# Erdos
invisible(COVID_ld(M,mc,N,d,-1,pD,rR,pI,tStop,q,nIs,'ER',oLD,tV,pLD,mLD,dLD))
# PLD
invisible(COVID_ld(M,mc,N,d,0.5,pD,rR,pI,tStop,q,nIs,'PLD',oLD,tV,pLD,mLD,dLD))

# Vaccination 

tV= 0.1*N
D= 0.1*N
who = c('S','I','R')
pE = 0.9

# Uniform Vaccination groups 'S' 'I' 'R'
oV= 1
# Erdos
invisible(COVID_vacc(M,mc,N,d,-1,pD,rR,pI,t,q,nIs,'ER',tV,D,oV,who))
# PLD
invisible(COVID_vacc(M,mc,N,d,0.5,pD,rR,pI,t,q,nIs,'PLD',tV,D,oV,who))

# Most popular Vaccination groups 'S' 'I' 'R'
oV= 4
# Erdos
invisible(COVID_vacc(M,mc,N,d,-1,pD,rR,pI,t,q,nIs,'ER',tV,D,oV,who))
# PLD
invisible(COVID_vacc(M,mc,N,d,0.5,pD,rR,pI,t,q,nIs,'PLD',tV,D,oV,who))

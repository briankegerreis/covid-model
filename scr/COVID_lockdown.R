### Lockdown 
# oLD:
#     1 - Lockdown type 1 Delete edges according to bernoulli pLD
#     2 - Lockdown type 2 Bound all the vertices connections to mLD
# oB: Option to begin
#     1 - Begin lockdown if the number of infections is greater than dnI
# oF: Future parameter
# dnI: Number of individuals infected to start the lockdown or time in which start the lockdown
# dLD: Number of days in lockdown
# mLD: Bound for lookdown
# pLD: Probability of delee an edge or select a vertix

prune <- function(g,M,i){
  edgs <- igraph::E(g)[from(i)]
  if(length(edgs)>M){
    redg <- sample(edgs,length(edgs)-M)
    g <- igraph::delete_edges(g,redg)
  }
  return(g)
}

graphPrune1 <- function(g,p){
  m <- rbinom(1,length(E(g)),p)
  g <- igraph::delete.edges(g,sample(E(g),m))
  return(g)
}


graphPrune2 <- function(g,M){
  ix<-which(degree(g)>M)
  for (i in ix) {
    g <- prune(g,M,i)
  }
  return(g)
}


graphPrune <- function(g,M,p,oLD){
  if(oLD == 1){
    g <- graphPrune1(g,p)
  }
  else{
    g <- graphPrune2(g,M)
  }
}

check_begin <- function(t,nI,dnI){
  return(nI > dnI)
}

check_finish <- function(t,tlD,dLD){
  return(t > tlD+dLD)
}

covid_ld <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,oLD,dnI,pLD,mLD,dLD,ij){
  hiaux <- c()
  while (T) {
    if(lamb > 0){
      typ <- rep(0,N)
      g <- block_power_matrix(matrix(d/N),lamb,typ)
    }
    else{g <- igraph::erdos.renyi.game(N,d/N)}
    ns <- setdiff(1:N,igraph::V(g))
    if(length(ns)>0){
      g <- igraph::add_vertices(g,length(ns))
    }
    nb <- igraph::adjacent_vertices(g,1:N)
    dg <- igraph::degree(g)
    qM <- quantile(dg,c(q))
    df <- data.frame(ID = 1:N,
                     typ = rep('S',N),
                     tI = rep(NA,N),
                     tR = rep(NA,N),
                     f = rep(NA,N),
                     ch = rep(0,N),
                     lt = rep(NA,N),
                     dg = dg,
                     stringsAsFactors = F)
    k <- 1
    rI <- adjustrI(pI,rR)
    str <- ' '
    li = c(sample(which(dg >0),1))
    lrI = c(rI)
    ls = list(nb[[li]])
    lsl = c(length(nb[[li]]))
    df$f[li] <- 0
    df$typ[li] <- 'I'
    df$lt[li] <- 0
    df$tI[li] <- 0
    str <- paste0(' ',li,' ;')
    t <- 0
    ldown <- F
    ldfin <- F
    nI <-1
    while(T){
      if(!ldown && check_begin(t,nI,dnI)){
        g1 <- graphPrune(g,mLD,pLD,oLD)
        nb <- igraph::adjacent_vertices(g1,igraph::V(g1))
        ls <- lapply(li, function(i) nb[[i]][df$typ[nb[[i]]] == 'S'])
        lsl <- sapply(ls, function(ni) length(ni))
        tLD <- t
        ldown <- T
      }
      if(ldown && !ldfin && check_finish(t,tLD,dLD)){
        nb <- igraph::adjacent_vertices(g,igraph::V(g))
        ls <- lapply(li, function(i) nb[[i]][df$typ[nb[[i]]] == 'S'])
        lsl <- sapply(ls, function(ni) length(ni))
        ldfin <- T
      }
      lenI <- length(lrI)
      rateInfect <- sum(lrI*lsl)
      rateReco <- rR*lenI
      ratemin <- rateInfect + rateReco
      if(ratemin == 0){
        break
      }
      else{
        texp <- rexp(1,rate = ratemin)
        t <- t+texp
        if(tStop < t){
          break
        }
        if(runif(1) < rateInfect/ratemin){
          nI <- nI+1
          ix <- if(length(lrI)>1){sample(1:lenI,1,prob = rI*lsl/rateInfect)}else{1}
          ni <- if(length(ls[[ix]]) == 1){ls[[ix]][1]}else{sample(ls[[ix]],1)}
          lni <- unlist(lapply(li,function(i){if(ni %in% nb[[i]]){i}}))
          df$f[ni] <- li[ix]
          df$typ[ni] <- 'I'
          df$tI[ni] <- t
          df$ch[df$f[ni]]<-df$ch[df$f[ni]]+1 
          str <- newick(df$f[ni],ni,t-df$lt[df$f[ni]],str)
          df$lt[ni] <- t
          df$lt[df$f[ni]] <-t 
          inix <-which(li %in% lni) 
          ls[inix] = lapply(inix, function(i){ls[[i]][ls[[i]] != ni]})
          lsl[inix] <- lsl[inix]-1
          li <- c(li,ni)
          lrI <- c(lrI,rI)
          ls = append(ls,list(nb[[ni]][df$typ[nb[[ni]]] == 'S']))
          lsl = c(lsl,length(nb[[ni]][df$typ[nb[[ni]]] == 'S']))
        }
        else{
          ix <- if(length(lrI)>1){sample(1:lenI,1)}else{1}
          df$typ[li[ix]] <-  if(rbinom(1,1, probOfDeath(pD,df$dg[li[ix]],qM))==1){'D'}else{'R'}
          df$tR[li[ix]] <- t
          str <- newick(li[ix],F,t-df$lt[li[ix]],str)
          li <- li[-ix]
          ls <- ls[-ix]
          lsl <- lsl[-ix]
          lrI <- lrI[-ix]
        }
      }
      k <- k+1
    }
    hiaux <- c(hiaux,sum(df$typ != 'S'))
    if(nI  > nIs || t>tStop){
      break
    }
  }
  write.table(str, file = paste0(pre,N,'/LD','_',oLD,'/NW','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(tLD, file = paste0(pre,N,'/LD','_',oLD,'/tLD/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(data.frame(hi = hiaux,lvl = pre, N= N), 
              file = paste0(pre,N,'/LD','_',oLD,'/HI','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(df, 
              file = paste0(pre,N,'/LD','_',oLD,'/DF','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(which(df$typ == 'S'), 
              file = paste0(pre,N,'/LD','_',oLD,'/G','/S',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(data.frame(dg = dg,lvl = pre, N= N ),
              file = paste0(pre,N,'/LD','_',oLD,'/DG','/',ij,'.txt'),
              sep ='\t', row.names = F,col.names = F,append = F)
  write.table(igraph::as_edgelist(g),
              file = paste0(pre,N,'/LD','_',oLD,'/G','/',ij,'.txt'),sep='\t',
              row.names = F,col.names = F,append = T)
  write.table( dfSIR(dfpreSIR(df),N),
               paste0(pre,N,'/LD','_',oLD,'/SIR','/',ij,'.txt'),
               sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
}

COVID_ld <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,oLD,dnI,pLD,mLD,dLD){
  dir.create(paste0(pre,N))
  dir.create(paste0(pre,N,'/LD','_',oLD))
  dir.create(paste0(pre,N,'/LD','_',oLD,'/tLD'))
  dir.create(paste0(pre,N,'/LD','_',oLD,'/HI'))
  dir.create(paste0(pre,N,'/LD','_',oLD,'/DG'))
  dir.create(paste0(pre,N,'/LD','_',oLD,'/G'))
  dir.create(paste0(pre,N,'/LD','_',oLD,'/NW'))
  dir.create(paste0(pre,N,'/LD','_',oLD,'/DF'))
  dir.create(paste0(pre,N,'/LD','_',oLD,'/SIR'))
  if(mc == 0){
    invisible(lapply(1:M, function(i) 
      covid_ld(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,oLD,dnI,pLD,mLD,dLD,i) ))
  }
  else(
    invisible(parallel::mclapply(1:M,function(i) 
      covid_ld(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,oLD,dnI,pLD,mLD,dLD,i) ))
  )
}

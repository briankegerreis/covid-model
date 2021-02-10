covid <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,ij){
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
    nI <-1
    while(T){
      lenI <- length(lrI)
      rateInfect <- sum(lrI*lsl)
      rateReco <- rR*lenI
      ratemin <- rateInfect + rateReco
      if(ratemin == 0){
        break
      }
      else{
        texp <- rexp(1,rate = ratemin)
        if(tStop < t){
          break
        }
        t <- t+texp
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
  write.table(str, file = paste0(pre,N,'/','NW','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(data.frame(hi = hiaux,lvl = pre, N= N), 
              file = paste0(pre,N,'/','HI','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(df, 
              file = paste0(pre,N,'/','DF','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(which(df$typ == 'S'), 
              file = paste0(pre,N,'/','G','/S',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(data.frame(dg = dg,lvl = pre, N= N ),file = paste0(pre,N,'/','DG','/',ij,'.txt'),
              sep ='\t', row.names = F,col.names = F,append = F)
  write.table(igraph::as_edgelist(g),file = paste0(pre,N,'/','G','/',ij,'.txt'),sep='\t',
              row.names = F,col.names = F,append = T)
  write.table(which(df$typ == 'S'), 
              file = paste0(pre,N,'/','G','/S',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table( dfSIR(dfpreSIR(df),N),file = paste0(pre,N,'/','SIR','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
}

COVID_control <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre){
  dir.create(paste0(pre,N))
  dir.create(paste0(pre,N,'/','HI'))
  dir.create(paste0(pre,N,'/','DG'))
  dir.create(paste0(pre,N,'/','G'))
  dir.create(paste0(pre,N,'/','NW'))
  dir.create(paste0(pre,N,'/','DF'))
  dir.create(paste0(pre,N,'/','SIR'))
  if(mc == 0){
    invisible(lapply(1:M, function(i) covid(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,i)))
  }
  else(
    invisible(parallel::mclapply(1:M,function(i) covid(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,i),mc.cores = mc))
  )
}

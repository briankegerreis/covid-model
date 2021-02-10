
powB <- function(B,lamb){
  n <- dim(B)[1]
  BP <-matrix(rep(0,n**2),nrow = n)
  for (i in 1:n) {
    for (j in i:n) {
      f <- function(x) x/B[i,j]-(lamb-log(x))^2/lamb^2
      rot <- uniroot(f, c(0, 1), tol = 1e-10)
      BP[i,j] = BP[j,i] = rot$root
    }
  }
  return(BP)
}

block_power_matrix <- function(B,lamb,typ){
  lam <- rexp(n = length(typ),rate = lamb )
  BP <- powB(B,lamb)
  le <- blockListCpp(BP,lam,typ)
  return(igraph::graph_from_edgelist(t(do.call(rbind,le)),directed = F))
}

# SIR 

# Add new i daugther j. 
# texp: time since the last i reproduction
# str:  netwik tree

# Add time remainig for i recuperation if j = F

newick <- function(i,j = F,texp,str){
  str1 <- strsplit(str,paste0(' ',i,' '))
  if(j != F){
    str <- paste0(str1[[1]][1],paste0('( ',i,' , ',j,' ):',texp,' '),str1[[1]][2])
  }
  else{
    str <- paste0(str1[[1]][1],paste0(' ',i,':',texp,' '),str1[[1]][2])
  }
  return(str)
}

# Calculate rate of infecti?n necesay to have pI probability of infection given that
# rR is the rate of recovery

adjustrI <- function(pI,rR){
  return((pI/(1-pI))*rR)
}

# Calculate probability of infection in terms of infection and recovery rate.

infectP <- function(rI,rR){
  return(rI/(rI+rR))
}

probOfDeath<- function(pD,dg,qM){
  if(dg>qM){
    return(pD[1])
  }else{
    return(pD[2])
  }
}


act <- function(i,dfsir,df,ncl){
  if(df[i,1] == 'R'){
    indx <- 3
  }
  else if(df[i,1] == 'D'){
    indx <- 4}
  else if(df[i,1] == 'I'){
    dfsir$S[(df[i,2]+1):ncl] <- dfsir$S[(df[i,2]+1):ncl]-1
    dfsir$I[(df[i,2]+1):ncl] <- dfsir$I[(df[i,2]+1):ncl]+1
    return(dfsir)
  }
  else{
    dfsir$S[(df[i,2]+1):ncl] <- dfsir$S[(df[i,2]+1):ncl]-1
    return(dfsir)
  }
  dfsir$S[(df[i,2]+1):ncl] <- dfsir$S[(df[i,2]+1):ncl]-1
  dfsir$I[(df[i,2]+1):(df[i,3]+1)] <- dfsir$I[(df[i,2]+1):(df[i,3]+1)]+1
  dfsir[(df[i,3]+2):ncl,indx] <- dfsir[(df[i,3]+2):ncl,indx]+1
  return(dfsir)
}

dfSIR <- function(df,N){
  ncl <- max(c(df[,3],df[,2]),na.rm = T)+2
  dfsir <- data.frame(S = rep(N,ncl),I = rep(0,ncl),R = rep(0,ncl), 
                      D = rep(0,ncl))
  
  for(i in 1:dim(df)[1] ){dfsir <- act(i,dfsir,df,ncl)}
  return(dfsir)  
} 

dfpreSIR <- function(df){
  df1 <- df[which(df$typ != 'S'),c(2,3,4)]
  df1[,2] <- ceiling(df1[,2])
  df1[,3] <- ceiling(df1[,3])
  return(df1)
}


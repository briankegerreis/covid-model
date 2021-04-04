# M: number of sims
# mc: number of cores
# N: population
# d: expected degree
# lamb: network degree tuning parameter
# pD: [p1, p2], where p1 is death probability for highly connected individuals, p2 is for less connected indivs
# rR: recovery rate
# pI: prob of infection
# t: max time
# q: fraction of vulnerable indivs
# nIs: minimum number of infections for a successful sim
# pre: output folder prefix
# ij: I think this keeps track of which sim is running

#covid <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,ij){
covid = function(num_sim, num_cores,
                 pop_size, d_connect, lamb,
                 p_death, recovery_rate, p_infect,
                 t_max, f_vulnerable, min_cases, pre, ij){
  hiaux <- c()
  while (T) {
    # set up the graph
    if(lamb > 0){
      typ <- rep(0,pop_size)
      g <- block_power_matrix(matrix(d_connect/pop_size),lamb,typ)
    }
    else{g <- igraph::erdos.renyi.game(pop_size,d_connect/pop_size)}
    ns <- setdiff(1:pop_size,igraph::V(g))
    if(length(ns)>0){
      g <- igraph::add_vertices(g,length(ns))
    }
    nb <- igraph::adjacent_vertices(g,1:N)
    neighbors_list = adjacent_vertices(g, 1:pop_size)
    dg <- igraph::degree(g)
    vulnerable_quantile <- quantile(dg,c(f_vulnerable))
    # typ: (S)usceptible, (I)nfected, (R)ecovered, (D)ead
    # tI: time of infection
    # tR: 
    # f: who they got it from?
    # ch: how many people they infected
    # lt: 
    # dg: degree
    #df <- data.frame(ID = 1:N,
    #                 typ = rep('S',N),
    #                 tI = rep(NA,N),
    #                 tR = rep(NA,N),
    #                 f = rep(NA,N),
    #                 ch = rep(0,N),
    #                 lt = rep(NA,N),
    #                 dg = dg,
    #                 stringsAsFactors = F)
    # can wrap these in a function like initialize_base_graph(g, params)
    set_vertex_attr(g, "typ", value=rep("S",pop_size))
    set_vertex_attr(g, "t_infected", value=rep(NA,pop_size))
    set_vertex_attr(g, "t_recovered", value=rep(NA,pop_size))
    set_vertex_attr(g, "parent", value=rep(NA,pop_size))
    set_vertex_attr(g, "n_children", value=rep(0,pop_size))
    set_vertex_attr(g, "lt", value=rep(NA,pop_size))
    set_vertex_attr(g, "degree", value=dg)
    set_vertex_attr(g, "p_infect", value=rep(0,pop_size))
    set_vertex_attr(g, "recovery_rate", value=rep(0,pop_size))
    set_vertex_attr(g, "infection_rate", value=rep(0,pop_size))
    k <- 1
    patient_zero = sample(which(dg>0),1) # randomly select patient zero
    set_vertex_attr(g, "typ", patient_zero, "I")
    # these two could be attributes of new directed edges
    set_vertex_attr(g, "t_infected", patient_zero, 0)
    set_vertex_attr(g, "parent", patient_zero, 0)
    # don't know what this is
    set_vertex_attr(g, "lt", patient_zero, 0)
    # these three could be "inherited" upon infection, plus some variation due to mutations
    set_vertex_attr(g, "p_infect", patient_zero, p_infect)
    set_vertex_attr(g, "recovery_rate", patient_zero, recovery_rate)
    set_vertex_attr(g, "infection_rate", patient_zero, adjustrI(p_infect, recovery_rate))
    #susceptible_neighbors = adjacent_vertices(g, patient_zero)
    #susceptible_neighbor_counts = length(susceptible_neighbors)
    str <- paste0(' ',patient_zero,' ;')
    t <- 0
    total_infections = 1
    while(T){
      infected_patients = which(vertex_attr(g, "typ")=="I")
      num_infected = length(infected_patients)
      infection_rates = vertex_attr(g, "infection_rate", infected_patients)
      # probably not, need to loop
      s_neighbor_list = vector("list", num_infected)
      for (i in 1:length(infected_patients)) {
        all_neighbors = neighbors(g, infected_patients[i])
        s_neighbors = all_neighbors[which(vertex_attr(g, "typ", all_neighbors)=="S")]
        s_neighbor_list[[i]] = s_neighbors
      }
      s_neighbor_counts = lengths(s_neighbor_list)
      rateInfect <- sum(infection_rates*s_neighbor_counts) # expected infections today: infection rate * number of neighbors
      recovery_rates = vertex_attr(g, "recovery_rate", infected_patients)
      #rateReco <- rR*lenI # expected recoveries today: recovery rate * number of infected
      rateReco = sum(recovery_rates)
      ratemin <- rateInfect + rateReco
      if(ratemin == 0){
        break
      }
      else{
        texp <- rexp(1,rate = ratemin) # I think this is a modified time step based on the number of people
        if(t > t_max){
          break
        }
        t <- t+texp
        # 03APR2021: made it here
        if(runif(1) < rateInfect/ratemin){ # an infection occurs
          nI <- nI+1
          ix <- if(length(lrI)>1){sample(1:lenI,1,prob = rI*lsl/rateInfect)}else{1} # find the spreader, weighted by connectivity
          ni <- if(length(ls[[ix]]) == 1){ls[[ix]][1]}else{sample(ls[[ix]],1)} # find the new case among the spreader's neighbors
          lni <- unlist(lapply(li,function(i){if(ni %in% nb[[i]]){i}})) # don't know what this does
          df$f[ni] <- li[ix] # assign who spread it to new case
          df$typ[ni] <- 'I' # change status to infected
          df$tI[ni] <- t # assigned time of infection
          df$ch[df$f[ni]]<-df$ch[df$f[ni]]+1 # add 1 to spreader's infected count
          str <- newick(df$f[ni],ni,t-df$lt[df$f[ni]],str) # I think this makes a tree to keep track on infections
          df$lt[ni] <- t # don't know what this means
          df$lt[df$f[ni]] <-t # or this
          inix <-which(li %in% lni) # or this
          ls[inix] = lapply(inix, function(i){ls[[i]][ls[[i]] != ni]}) # or this
          lsl[inix] <- lsl[inix]-1 # or this
          li <- c(li,ni) # add new case to the list of infected
          lrI <- c(lrI,rI) # add new infection rate to the list of infected
          ls = append(ls,list(nb[[ni]][df$typ[nb[[ni]]] == 'S'])) # add new susceptible neighbors
          lsl = c(lsl,length(nb[[ni]][df$typ[nb[[ni]]] == 'S'])) # add connectivity of new case
        }
        else{ # a recovery or death occurs
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
    hiaux <- c(hiaux,sum(df$typ != 'S')) # keep track of how many people are not susceptible
    if(nI  > nIs || t>tStop){ # arrive here by having zero new events at this time point or by running out of time. # restart sim if not enough infections happened
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
  write.table( dfSIR(dfpreSIR(df),N),file = paste0(pre,N,'/','SIR','/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  write.table(igraph::as_edgelist(g),
              file = paste0(pre,N,'/G','/',ij,'.txt'),sep='\t',
              row.names = F,col.names = F,append = T)
}

#COVID_control <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre){
COVID_control = function(num_sim, num_cores,
                         pop_size, d_connect, lamb,
                         p_death, recovery_rate, p_infect,
                         t_max, f_vulnerable, min_cases, pre){
  dir.create(paste0(pre,pop_size))
  dir.create(paste0(pre,pop_size,'/','HI'))
  dir.create(paste0(pre,pop_size,'/','DG'))
  dir.create(paste0(pre,pop_size,'/','G'))
  dir.create(paste0(pre,pop_size,'/','NW'))
  dir.create(paste0(pre,pop_size,'/','DF'))
  dir.create(paste0(pre,pop_size,'/','SIR'))
  if(mc == 0){
    #invisible(lapply(1:M, function(i) covid(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,i)))
    invisible(lapply(1:num_sim, function(i) covid(num_sim,num_cores,pop_size,d_connect,lamb,p_death,recovery_rate,p_infect,t_max,f_vulnerable,min_cases,pre,i)))
  }
  else(
    #invisible(parallel::mclapply(1:M,function(i) covid(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,i),mc.cores = mc))
    invisible(parallel::mclapply(1:num_sim,function(i) covid(num_sim,num_cores,pop_size,d_connect,lamb,p_death,recovery_rate,p_infect,t_max,f_vulnerable,min_cases,pre,i),mc.cores = num_cores))
  )
}

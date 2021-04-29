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
                 t_max, f_vulnerable, min_cases, pre, ij) {
  hiaux <- c()
  while (T) {
    # set up the graph
    if(lamb > 0) {
      typ <- rep(0,pop_size)
      g <- block_power_matrix(matrix(d_connect/pop_size),lamb,typ)
    }
    else {g <- igraph::erdos.renyi.game(pop_size,d_connect/pop_size)}
    ns <- setdiff(1:pop_size,igraph::V(g))
    if(length(ns)>0) {
      g <- igraph::add_vertices(g,length(ns))
    }
    nb <- igraph::adjacent_vertices(g,1:N)
    neighbors_list = adjacent_vertices(g, 1:pop_size)
    dg <- igraph::degree(g)
    vulnerable_quantile <- quantile(dg,c(f_vulnerable))
    # these could be set in a yaml config file
    pop_attrs = list(typ="S", t_infected=NA, t_resolved=NA, parent=NA, n_children=0, lt=NA, p_infect=0, recovery_rate=0, infection_rate=0,
                     vax=FALSE, vax_efficacy=0)
    #personal_attrs = list(degree=dg, p_death=ifelse(dg<vulnerable_quantile, p_death[2], p_death[1]))
    personal_attrs = list(degree=dg, risk=ifelse(dg<vulnerable_quantile, "hi", "lo"))
    g = initialize_population_attributes(g, pop_attrs) %>%
      initialize_person_attributes(g, personal_attrs)
    k <- 1
    patient_zero = Sample(which(dg>0),1) # randomly select patient zero
    # these could be set in a yaml config file
    covid_attrs = list(p_infect=p_infect, recovery_rate=recovery_rate, p_death_lo=p_death[1], p_death_hi=p_death[2], vax_resistance=0)
    g = initialize_patient_zero(g, patient_zero, covid_attrs)
    str <- paste0(' ',patient_zero,' ;')
    t <- 0
    total_infections = 1
    while(T) {
      #### g = vaccinate(g, vacc_rate, strategy, test_for_eligibility) # vaccinate 1% per day, either most connected or random, OK vaccinating R patients or S only?
      infected_patients = which(vertex_attr(g, "typ")=="I")
      num_infected = length(infected_patients)
      rates = calculate_infection_recovery_rates(g, infected_patients)
      rateInfect = rates$rateInfect
      rateReco = rates$rateReco
      ratemin <- rateInfect + rateReco
      if(ratemin == 0) {
        break # the outbreak ended
      }
      else {
        texp <- rexp(1,rate = ratemin) # I think this is a modified time step based on the number of people
        if(t > t_max) {
          break
        }
        t <- t+texp
        if(runif(1) < rateInfect/ratemin) { # an infection occurs
          total_infections = total_infections + 1
          # because all of the infection and spread coefficients can't be vertex attributes, need to first find which patient spreads the disease (spreader_ix), then figure out which member of the population that is (spreader)
          spreader_ix = Sample(1:length(infected_patients), 1, prob=rates$infection_coefficients/rateInfect) # find the spreader, weighted by infection rate, connection strength to susceptible neighbors, and susceptibility of neighbors
          spreader = infected_patients[spreader_ix]
          catch_probabilities = rates$catch_coefficients[[spreader_ix]]/sum(rates$catch_coefficients[[spreader_ix]]
          new_case = Sample(susceptible_neighbors(g, spreader), 1, prob=catch_probabilities)) # find the new case among the spreader's susceptible neighbors
          #### g = infect_and_mutate(g, spreader, new_case, mutation_prob, mutation_rate)
          g = infect_new_case(g, spreader, new_case, t, names(covid_attr))
          # str <- newick(df$f[ni],ni,t-df$lt[df$f[ni]],str) # I think this makes a tree to keep track of infections
        }
        else { # a recovery or death occurs
          # is this right?
          patient = Sample(infected_patients, 1, prob=get_vertex_attr(g,"recovery_rate",infected_patients)/rateReco)
          outcome = if(rbinom(1,1, vertex_attr(g, "p_death", patient)==1){"D"}else{"R"}
          g = set_vertex_attr(g, "typ", patient, outcome)
          g = set_vertex_attr(g, "t_resolved", patient, t)
          # str <- newick(li[ix],F,t-df$lt[li[ix]],str)
        }
      }
      k <- k+1
    }
    hiaux <- c(hiaux,length(which(vertex_attr(g, "typ")!="S")) # keep track of how many people are not susceptible
    if(nI  > nIs || t>tStop) { # arrive here by having zero new events at this time point or by running out of time. # restart sim if not enough infections happened
      break
    }
  }
  # write.table(str, file = paste0(pre,N,'/','NW','/',ij,'.txt'),
  #             sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  # write.table(data.frame(hi = hiaux,lvl = pre, N= N), 
  #             file = paste0(pre,N,'/','HI','/',ij,'.txt'),
  #             sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  # write.table(df, 
  #             file = paste0(pre,N,'/','DF','/',ij,'.txt'),
  #             sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  # write.table(which(df$typ == 'S'), 
  #             file = paste0(pre,N,'/','G','/S',ij,'.txt'),
  #             sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  # write.table(data.frame(dg = dg,lvl = pre, N= N ),file = paste0(pre,N,'/','DG','/',ij,'.txt'),
  #             sep ='\t', row.names = F,col.names = F,append = F)
  # write.table(igraph::as_edgelist(g),file = paste0(pre,N,'/','G','/',ij,'.txt'),sep='\t',
  #             row.names = F,col.names = F,append = T)
  # write.table( dfSIR(dfpreSIR(df),N),file = paste0(pre,N,'/','SIR','/',ij,'.txt'),
  #             sep = "\t",row.names = FALSE,col.names = FALSE,append=F)
  # write.table(igraph::as_edgelist(g),
  #             file = paste0(pre,N,'/G','/',ij,'.txt'),sep='\t',
  #             row.names = F,col.names = F,append = T)
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

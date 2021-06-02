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
                 t_max, f_vulnerable, min_cases, pre, 
                 mutation, vaccination, ij) {
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
    nb <- igraph::adjacent_vertices(g,1:pop_size)
    neighbors_list = adjacent_vertices(g, 1:pop_size)
    dg <- igraph::degree(g)
    vulnerable_quantile <- quantile(dg,c(f_vulnerable))
    # these could be set in a yaml config file
    pop_attrs = list(typ="S", t_infected=NA, t_resolved=NA, parent=NA, n_children=0, lt=NA,
                     p_infect=0, recovery_rate=0, infection_rate=0,
                     susceptibility=1, vax=FALSE, vax_efficacy=0)
    #personal_attrs = list(degree=dg, p_death=ifelse(dg<vulnerable_quantile, p_death[2], p_death[1]))
    personal_attrs = list(degree=dg, risk=ifelse(dg<vulnerable_quantile, "hi", "lo"))
    edge_attrs = list(weight=1)
    g = initialize_population_attrs(g, pop_attrs) %>%
      initialize_personal_attrs(personal_attrs) %>%
      initialize_edge_attrs(edge_attrs)
    k <- 1
    patient_zero = Sample(which(dg>0),1) # randomly select patient zero
    covid_attrs = list(p_infect=p_infect, recovery_rate=recovery_rate,
                       p_death_lo=p_death[1], p_death_hi=p_death[2], vax_resistance=0,
                       strain_id=1)
    g = initialize_patient_zero(g, patient_zero, covid_attrs)
    newick_tree <- paste0(' ',patient_zero,' ;')
    newick_list = vector("list", pop_size) # make a list that is longer than it needs to be to hold all possible strains
    newick_list[[1]] = newick_tree
    strain_id = 2
    next_strain_id = 2
    t <- 0
    total_infections = 1
    vax_done = FALSE
    while(T) {
      #### g = vaccinate(g, vacc_rate, strategy, test_for_eligibility) # vaccinate 1% per day, either most connected or random, OK vaccinating R patients or S only?
      infection_status = vertex_attr(g, "typ")
      infected_patients = which(infection_status=="I")
      num_infected = length(infected_patients)
      if (num_infected==0) {
        break
      }
      if (vaccination$enabled & !vax_done) {
        if (total_infections/pop_size >= vaccination$threshold_fraction) {
          vax_result = vaccination_control(g, vaccination$fraction_to_vaccinate*pop_size,
                                           vaccination$fraction_to_vaccinate*pop_size,
                                           pop_size, vaccination$strategy,
                                           vaccination$eligible_groups,
                                           vaccination$efficacy)
          g = vax_result$g
          vax_done = TRUE
          print(paste0("vaccinated at t: ",t))
        }
      }
      rates = calculate_infection_recovery_rates(g, infected_patients, neighbors_list)
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
          spreader_ix = Sample(1:length(infected_patients), 1, prob=rates$spread_coefficients/rateInfect) # find the spreader, weighted by infection rate, connection strength to susceptible neighbors, and susceptibility of neighbors
          spreader = infected_patients[spreader_ix]
          catch_probabilities = rates$catch_coefficients[[spreader_ix]]/sum(rates$catch_coefficients[[spreader_ix]])
          new_case = Sample(susceptible_neighbors(g, spreader, neighbors_list, infection_status), 1, prob=catch_probabilities) # find the new case among the spreader's susceptible neighbors
          #### g = infect_and_mutate(g, spreader, new_case, mutation_prob, mutation_rate)
          # move newick here before we set spreader.lt to t
          # tree_list
          # with original strain and each new strain, add a strain_id
          # use strain_id to modify tree_list[[strain_id]]
          # then write a tree file for each strain at the end
          # and maintain one overall tree
          newick_tree = newick(spreader, new_case, t-vertex_attr(g,"lt",spreader), newick_tree)
          newick_list[[vertex_attr(g,"strain_id",spreader)]] = newick(spreader, new_case, t-vertex_attr(g,"lt",spreader), newick_list[[vertex_attr(g,"strain_id",spreader)]])
          g = infect_new_case(g, spreader, new_case, t, names(covid_attrs))
          if (mutation$enabled) {
            if (mutation$mode=="random") {
              mut = mutation$random_params
              mutation_result = mutation_control(g, new_case, mut$p_mutate_infection, mut$p_mutate_death,
                                                 mut$p_increase_infection, mut$p_increase_death,
                                                 mut$infection_step, mut$death_step_lo, mut$death_step_hi,
                                                 mut$p_corr, mut$force_opposite_signs, mut$growth, strain_id)
              g = mutation_result$g
              next_strain_id = mutation_result$next_strain_id
            } else if (mutation$mode=="created" & total_infections == mutation$created_params$n_infections) {
              mut = mutation$created_params
              mutation_result = create_variant(g, new_case, mut$p_infect,
                                               mut$p_death_lo, mut$p_death_hi, strain_id)
              g = mutation_result$g
              next_strain_id = mutation_result$next_strain_id
              print(paste0("created a mutant at t: ",t))
            }
          }
          if (next_strain_id > strain_id) {
            newick_list[[strain_id]] = paste0(' ',new_case,' ;')
            strain_id = next_strain_id
          }
        }
        else { # a recovery or death occurs
          # is this right?
          patient = Sample(infected_patients, 1, prob=vertex_attr(g,"recovery_rate",infected_patients)/rateReco)
          if (rbinom(1,1, vertex_attr(g, "p_death", patient))==1) {
            outcome = "D"
          } else {
            outcome = "R"
          }
          g = set_vertex_attr(g, "typ", patient, outcome)
          g = set_vertex_attr(g, "t_resolved", patient, t)
          # str <- newick(li[ix],F,t-df$lt[li[ix]],str)
          # newick tree can stay here
          newick_tree = newick(patient, FALSE, t-vertex_attr(g,"lt",patient), newick_tree)
          newick_list[[vertex_attr(g,"strain_id",patient)]] = newick(patient, FALSE, t-vertex_attr(g,"lt",patient), newick_list[[vertex_attr(g,"strain_id",patient)]])
        }
      }
      k <- k+1
    }
    hiaux <- c(hiaux,length(which(vertex_attr(g, "typ")!="S"))) # keep track of how many people are not susceptible
    if(total_infections  > min_cases || t>t_max) { # arrive here by having zero new events at this time point or by running out of time. # restart sim if not enough infections happened
      break
    }
  }
  df = as_data_frame(g,"vertices")
  write.table(newick_tree, file = paste0(pre,pop_size,'/newick/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=FALSE)
  for (i in 1:length(newick_list)) {
    if (!is.null(newick_list[[i]]))
    write.table(newick_list[[i]], file = paste0(pre,pop_size,'/newick/',ij,'_',i,'.txt'),
                sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
  }
  write.table(data.frame(hi = hiaux,lvl = pre, N = pop_size), 
              file = paste0(pre,pop_size,'/hiaux/',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=FALSE)
  write.table(df, 
              file = paste0(pre,pop_size,'/dataframe/',ij,'.txt'),
              sep = "\t",append=FALSE)
  write.table(which(df$typ == 'S'), 
              file = paste0(pre,pop_size,'/graph/S',ij,'.txt'),
              sep = "\t",row.names = FALSE,col.names = FALSE,append=FALSE)
  write.table(data.frame(dg = dg,lvl = pre, N = pop_size ),
              file = paste0(pre,pop_size,'/degree/',ij,'.txt'),
              sep ='\t', row.names = F,col.names = F,append = FALSE)
  write.table(igraph::as_edgelist(g), file = paste0(pre,pop_size,'/graph/',ij,'.txt'),sep='\t',
              row.names = F,col.names = F,append = TRUE)
  write.table(dfSIR(df,pop_size), file = paste0(pre,pop_size,'/SIR/',ij,'.txt'),
               sep = "\t",row.names = FALSE,append=FALSE)
  to_return = list(graph=g, df_edges=as_data_frame(g,"edges"), df_vertices=as_data_frame(g,"vertices"))
  return(to_return)
}

#COVID_control <- function(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre){
COVID_control = function(num_sim, num_cores,
                         pop_size, d_connect, lamb,
                         p_death, recovery_rate, p_infect,
                         t_max, f_vulnerable, min_cases, pre,
                         mutation, vaccination){
  dir.create(paste0(pre,pop_size))
  dir.create(paste0(pre,pop_size,'/','hiaux'))
  dir.create(paste0(pre,pop_size,'/','degree'))
  dir.create(paste0(pre,pop_size,'/','graph'))
  dir.create(paste0(pre,pop_size,'/','newick'))
  dir.create(paste0(pre,pop_size,'/','dataframe'))
  dir.create(paste0(pre,pop_size,'/','SIR'))
  if(num_cores == 0){
    #invisible(lapply(1:M, function(i) covid(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,i)))
    invisible(lapply(1:num_sim, function(i) covid(num_sim,num_cores,pop_size,d_connect,lamb,p_death,recovery_rate,p_infect,t_max,f_vulnerable,min_cases,pre,mutation,vaccination,i)))
  }
  else(
    #invisible(parallel::mclapply(1:M,function(i) covid(M,mc,N,d,lamb,pD,rR,pI,t,q,nIs,pre,i),mc.cores = mc))
    invisible(parallel::mclapply(1:num_sim,function(i) covid(num_sim,num_cores,pop_size,d_connect,lamb,p_death,recovery_rate,p_infect,t_max,f_vulnerable,min_cases,pre,mutation,vaccination,i),mc.cores = num_cores))
  )
}

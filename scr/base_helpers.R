# override the behavior of sample() when first argument is a number instead of a vector
# sample() would try to sample from 1:x instead of just returning x
Sample <- function(vec,...){
  if(length(vec) == 1){
    return(vec[1])}
  else{
    return(sample(vec,...))
  }
}


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

# convenience function to calculate infection rate when someone gets sick
calculate_disease_rates = function(g, patient) {
  pI = vertex_attr(g, "p_infect", patient)
  rR = vertex_attr(g, "recovery_rate", patient)
  rI = adjustrI(pI,rR)
  g = set_vertex_attr(g, "infection_rate", patient, rI)
  return(g)
}

# convenience function to determine p_death when someone gets sick based on their risk group
determine_p_death = function(g, patient) {
  p_death_lo = vertex_attr(g, "p_death_lo", patient)
  p_death_hi = vertex_attr(g, "p_death_hi", patient)
  p_death = ifelse(vertex_attr(g, "risk", patient)=="lo", p_death_lo, p_death_hi)
  g = set_vertex_attr(g, "p_death", patient, p_death)
  return(g)
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

# assign attributes to population
initialize_population_attrs = function(g, attr_list) {
  # for things that start the same for the whole population
  # each attr has one value
  pop_size = length(V(g))
  for (attr in names(attr_list)) {
    g = set_vertex_attr(g, attr, value=rep(attr_list[[attr]],pop_size))
  }
  return(g)
}

# assign varying attributes to population
initialize_personal_attrs = function(g, attr_list) {
  # for things that vary across the population
  # each attr is the size of the population
  for (attr in names(attr_list)) {
    g = set_vertex_attr(g, attr, value=attr_list[[attr]])
  }
  return(g)
}

# assign attributes to edges
initialize_edge_attrs = function(g, attr_list) {
  n_edges = length(E(g))
  for (attr in names(attr_list)) {
    g = set_edge_attr(g, attr, value=rep(attr_list[[attr]],n_edges))
  }
  return(g)
}

# assign disease attributes to patient zero
# attributes include characteristics of disease like prob of infection and recovery rate, and modifiers to patient attributes like prob of death and vaccine efficacy
# example: modifier_p_death=1.5, modifier_vacc_efficacy=0.75 --> p_death goes from 0.01 to 0.015, vacc_efficacy goes from 0.80 to 0.60
# attr_list = list(p_infect=p_infect, recovery_rate=recovery_rate, modifier_p_death=1, modifier_vacc_efficacy=1)
initialize_patient_zero = function(g, patient_zero, attr_list) {
  g = set_vertex_attr(g, "typ", patient_zero, "I") %>%
    set_vertex_attr("t_infected", patient_zero, 0) %>%
    set_vertex_attr("parent", patient_zero, 0) %>%
    set_vertex_attr("lt", patient_zero, 0)
  # give patient zero default covid attributes
  for (attr in names(attr_list)) {
    g = set_vertex_attr(g, attr, patient_zero, attr_list[[attr]])
  }
  g = calculate_disease_rates(g, patient_zero)
  g = determine_p_death(g, patient_zero)
  return(g)
}

# determine if and how a random mutation occurs
# p_mutate_x is the probability of a change
# p_increase_x is the probability of in increase, given a change
# p_corr is the probability that infection and death both mutate at once, according to p_mutate_infection
# force_opposite_signs forces infection and death to mutate in opposite directions when the p_corr check succeeds
# mutate_array contains ones and zeros determing whether infection and death will mutate
# sign_array contains 1 and -1 determining whether the mutations are increases or decreases
# finally, calls mutation_function and returns the modified graph
mutation_control = function(g, patient, p_mutate_infection, p_mutate_death, p_increase_infection, p_increase_death,
                            infection_step, death_step_lo, death_step_hi, p_corr, force_opposite_signs,
                            growth=c("additive","multiplicative"), strain_id) {
  if (rbinom(1,1,p_corr)==1) { # if both characteristics are changing at once
    mutate = rbinom(1,1,p_mutate_infection)
    if (mutate==0) { # do nothing
      return(list(g=g, next_strain_id=strain_id))
    } else {
      mutate_array = rep(mutate,2) # set both characteristics to mutate
      if (force_opposite_signs) {
        s = rbinom(1,1,p_increase_infection)
        sign_array = c(s, 1-s)
      } else {
        sign_array = rbinom(2,1,c(p_increase_infection,p_increase_death))
      }
    }
  } else {
    mutate_array = rbinom(2,1,c(p_mutate_infection,p_mutate_death))
    sign_array = rbinom(2,1,c(p_increase_infection,p_increase_death))
    if (sum(mutate_array)==0) { # do nothing
      return(list(g=g, next_strain_id=strain_id))
    }
  }
  sign_array = 2*sign_array-1 # convert 0 to -1
  g = mutation_function2(g, patient, mutate_array, sign_array, infection_step, death_step_lo, death_step_hi, growth, strain_id)
  next_strain_id = strain_id + 1
  return(list(g=g, next_strain_id=next_strain_id))
}

# one possible way to mutate infection and death probabilities
# p_infect and p_death are proportional to modifier coefficients
# for additive mutation, coefficients increase linearly like c=c+x
# for multiplicative mutation, coefficients increase multiplicatively like c=c*(1+x)
# this makes it easier to accommodate multiple death rates, infection probabilities, and recovery rates (if we want) throughout the population
# but it is difficult to bound these coefficients, so the final probabilities have to be bounded after the initial calculation
mutation_function = function(g, patient, mutate_array, sign_array, infection_step, death_step, growth) {
  infection_delta = mutate_array[1]*sign_array[1]*infection_step
  death_delta = mutate_array[2]*sign_array[1]*death_step
  if (growth=="additive") {
    infection_coef_new = vertex_attr(g, "modifier_p_infect", patient) + infection_delta
    infection_coef_new = max(c(infection_step, infection_coef_new)) # not allowed to be 0
    death_coef_new = vertex_attr(g, "modifier_p_death", patient) + death_delta
    death_coef_new = max(c(0, death_coef_new)) # allowed to be 0
  } else if (growth=="multiplicative") {
    infection_coef_new = vertex_attr(g, "modifier_p_infect", patient) * (1+infection_delta) # cannot possibly reach 0
    death_coef_new = vertex_attr(g, "modifier_p_death", patient) * (1+death_delta) # cannot possibly reach 0
  }
  g = set_vertex_attr(g, "modifier_p_infect", patient, infection_coef_new)
  g = set_vertex_attr(g, "modifier_p_death", patient, death_coef_new)
  p_infect_new = vertex_attr(g, "p_infect", patient) * infection_coef_new
  if (growth=="additive") {
    p_infect_new = min(c(p_infect_new, 1-vertex_attr(g, "p_infect", patient)*infection_step)) # must remain at least 1 step below 100% infectious
  } else if (growth=="multiplicative") {
    p_infect_new = min(c(p_infect_new, 1/(1+vertex_attr(g, "p_infect", patient)*infection_step))) # must remain at least 1/(1+step) below 100% infectious
  }
  p_death_new = vertex_attr(g, "p_death", patient) * death_coef_new
  p_death_new = min(c(p_death_new, 1)) # must be <= 100% lethal
  g = set_vertex_attr(g, "p_infect", patient, p_infect_new)
  g = set_vertex_attr(g, "p_death", patient, p_death_new)
  g = set_vertex_attr(g, "infection_rate", patient, adjustrI(p_infect_new,vertex_attr(g, "recovery_rate", patient)))
  # how do I keep the probabilities between 0 and 1 when all I have is the modifier?
  # maybe I let the modifier go as high or low as I want but keep the bounds when I recalculate infection rate and death probability
  # what is the highest or lowest the infection prob can go?
  # additive: [step,1-step]
  # multiplicative: (0,1/(1+step)]
  return(g)
}

# another way to mutate infection and death probabilities
# p_infect and p_death are directly adjusted
# for additive mutation, they are increased by c = c +/- step_size
# for multiplicative mutation, they are increased by c = c * (1 +/- step_size)
# death probabilities:
#   p_death_lo and p_death_hi are death probabilities for low-risk and high-risk individuals
#   they have separate step sizes (they could be equal or different, but they must both be specified)
#   need a new attribute "risk" to determine which value goes into p_death
# this function would have to change if we want to add more infection probabilities or death probabilities in the population
# but the math is much easier to manage
# this meets immediate needs and will be implemented for now
mutation_function2 = function(g, patient, mutate_array, sign_array, infection_step, death_step_lo, death_step_hi, growth, strain_id) {
  infection_delta = mutate_array[1]*sign_array[1]*infection_step
  death_delta_lo = mutate_array[2]*sign_array[2]*death_step_lo
  death_delta_hi = mutate_array[2]*sign_array[2]*death_step_hi
  if (growth=="additive") {
    p_infect_new = vertex_attr(g, "p_infect", patient) + infection_delta
    p_infect_new = sort(c(infection_step, p_infect_new, 1-infection_step))[2]
    p_death_lo_new = vertex_attr(g, "p_death_lo", patient) + death_delta_lo
    p_death_lo_new = sort(c(0, p_death_lo_new, 1))[2]
    p_death_hi_new = vertex_attr(g, "p_death_hi", patient) + death_delta_hi
    p_death_hi_new = sort(c(0, p_death_hi_new, 1))[2]
  } else if (growth=="multiplicative") {
    p_infect_new = vertex_attr(g, "p_infect", patient) * (1+infection_delta)
    p_infect_new = min(c(p_infect_new, 1/(1+infection_step)))
    p_death_lo_new = vertex_attr(g, "p_death_lo", patient) * (1+death_delta_lo)
    p_death_lo_new = min(c(p_death_lo_new, 1))
    p_death_hi_new = vertex_attr(g, "p_death_hi", patient) * (1+death_delta_hi)
    p_death_hi_new = min(c(p_death_hi_new, 1))
  } else {
    stop("Mutation growth mode must be 'additive' or 'multiplicative'")
  }
  g = set_vertex_attr(g, "p_infect", patient, p_infect_new) %>%
      set_vertex_attr("p_death_lo", patient, p_death_lo_new) %>%
      set_vertex_attr("p_death_hi", patient, p_death_hi_new) %>%
      set_vertex_attr("strain_id", patient, strain_id)
  g = calculate_disease_rates(g, patient)
  g = determine_p_death(g, patient)
  return(g)
}

# use divine power to create a super strain in a newly infected patient
create_variant = function(g, patient, p_infect, p_death_lo, p_death_hi, strain_id) {
  g = set_vertex_attr(g, "p_infect", patient, p_infect) %>%
      set_vertex_attr("p_death_lo", patient, p_death_lo) %>%
      set_vertex_attr("p_death_hi", patient, p_death_hi) %>%
      set_vertex_attr("strain_id", patient, strain_id)
  g = calculate_disease_rates(g, patient)
  g = determine_p_death(g, patient)
  next_strain_id = strain_id + 1
  return(list(g=g, next_strain_id=next_strain_id))
}

# manage function calls for different vaccine strategies
# total_doses and available_doses can be used to vaccinate people over time, or they can be set to the same value to give all doses at once
vaccination_control = function(g, total_doses, available_doses, pop_size, strategy, eligible_groups, efficacy) {
  if (strategy=="neighbors") {
    l = vaccinate_neighbors(g, eligible_groups, total_doses, available_doses, efficacy)
  } else if (strategy=="uniform") {
    l = vaccinate_connected(g, eligible_groups, total_doses, available_doses, 1, "uniform", efficacy)
  } else if (strategy=="most_connected") {
    l = vaccinate_connected(g, eligible_groups, total_doses, available_doses, total_doses/pop_size, "most", efficacy)
  } else if (strategy=="among_most_connected") {
    l = vaccinate_connected(g, eligible_groups, total_doses, available_doses, 0.5, "most", efficacy)
  } else if (strategy=="least_connected") {
    l = vaccinate_connected(g, eligible_groups, total_doses, available_doses, total_doses/pop_size, "least", efficacy)
  } else if (strategy=="among_least_connected") {
    l = vaccinate_connected(g, eligible_groups, total_doses, available_doses, 0.5, "most", efficacy)
  }
  return(l)
}

# vaccinate people based on connectivity
# eligible_groups is a vector like c("S") or c("S","I","R")
# fraction_targeted is the proportional of the population targeted for vaccination based on connectivity
# strategy "most" targets highly connected people, strategy "least" targets disconnected people
# returns the graph and the number of remaining doses if we want to keep track of that
vaccinate_connected = function(g, eligible_groups, total_doses, available_doses, fraction_targeted, strategy=c("uniform", "most", "least"), efficacy) {
  eligible_people = which(vertex_attr(g, "typ") %in% eligible_groups)
  if (length(eligible_people) < available_doses) {
    available_doses = length(eligible_people)
  }
  dg = degree(g)
  eligible_dg = dg[eligible_people]
  target_quantile = quantile(dg, fraction_targeted)
  if (strategy=="uniform") {
    targeted_people = eligible_people
  } else if (strategy=="most") {
    targeted_people = eligible_people[which(eligible_dg)>=target_quantile]
  } else if (strategy=="least") {
    targeted_people = eligible_people[which(eligible_dg)<=target_quantile]
  }
  vaxxed_people = Sample(targeted_people, available_doses)
  g = vaccinate(g, vaxxed_people, efficacy)
  doses_remaining = total_doses - available_doses
  return(list(g=g, doses_remaining=doses_remaining))
}

# vaccination strategy based on neighbors
# works like vaccinate_connected()
vaccinate_neighbors = function(g, eligible_groups, total_doses, available_doses, efficacy) {
  eligible_people = which(vertex_attr(g, "typ") %in% eligible_groups)
  nb = adjacent_vertices(g)
  if (length(eligible_people) > available_doses) {
    available_doses = length(eligible_people)
  }
  doses_remaining = total_doses
  while (available_doses>0 & length(eligible_people)>0) {
    eligible_neighbors = sapply(eligible_people, function(i) sum(vertex_attr(g,"typ",nb[[i]]) %in% eligible_groups & vertex_attr(g,"vax",nb[[i]])==F))
    if (sum(eligible_neighbors)==0) {
      vaxxed_people = Sample(eligible_people, available_doses)
      g = vaccinate(g, vaxxed_people, efficacy)
      doses_remaining = doses_remaining - available_doses
      break
    }
    ix = Sample(eligible_people[which(eligible_neighbors>0)],1)
    ni = nb[[ix]][which(vertex_attr(g,"typ",nb[[ix]]) %in% eligible_groups & vertex_attr(g,"vax",nb[[ix]])==F)]
    vaxxed_person = Sample(ni,1)
    g = vaccinate(g, ni, efficacy)
    doses_remaining = doses_remaining - 1
    available_doses = available_doses - 1
    eligible_people = eligible_people[!eligible_people %in% vaxxed_person]
  }
  return(list(g=g, doses_remaining=doses_remaining))
}

# set vaccine attributes for newly vaccinated people
vaccinate = function(g, people, efficacy) {
  g = set_vertex_attr(g, "vax", people, rep(TRUE,length(people))) %>%
      set_vertex_attr("vax_efficacy", people, rep(efficacy,length(people)))
  return(g)
}

# assign disease attributes to neighbor upon infection
infect_new_case = function(g, spreader, new_case, t, covid_attr_names) {
  g = set_vertex_attr(g, "typ", new_case, "I") %>%
    set_vertex_attr("parent", new_case, spreader) %>%
    set_vertex_attr("t_infected", new_case, t) %>%
    set_vertex_attr("lt", c(spreader,new_case), rep(t,2)) %>%
    set_vertex_attr("n_children", spreader, vertex_attr(g, "n_children", spreader)+1)
  # give new_case spreader's covid attributes
  for (attr in covid_attr_names) {
    g = set_vertex_attr(g, attr, new_case, vertex_attr(g, attr, spreader))
    # # and apply multiplicative modifiers to target attributes
    # if startsWith(attr, "modifier_") {
      # target_attr = sub("modifier_", "", attr)
      # g = set_vertex_attr(g, target_attr, new_case, vertex_attr(g, target_attr, new_case)*vertex_attr(g, attr, new_case))
    # }  
  }
  pI = vertex_attr(g, "p_infect", new_case)
  rR = vertex_attr(g, "recovery_rate", new_case)
  g = set_vertex_attr(g, "infection_rate", new_case, adjustrI(pI,rR))
  p_death = ifelse(vertex_attr(g, "risk", new_case)=="lo", vertex_attr(g, "p_death_lo", new_case), vertex_attr(g, "p_death_hi", new_case))
  g = set_vertex_attr(g, "p_death", new_case, p_death)
  return(g)
}

# find all susceptible neighbors of a given individual
susceptible_neighbors = function(g, x, neighbors_list, infection_status) {
  # all_neighbors = igraph::neighbors(g, x)
  all_neighbors = neighbors_list[[x]]
  s_neighbors = all_neighbors[which(infection_status[all_neighbors]=="S")]
  return(as.integer(s_neighbors)) # no need to deal with the 'igraph.vs' class
}

# convert individual and vector of neighbors to c(indiv, n1, indiv, n2, indiv, n3...) to extract edges later
generate_vertex_list = function(x, ns) {
  if (length(ns)==0) {
    return(integer(0))
  }
  vertex_list = integer(0)
  for (i in 1:length(ns)) {
    vertex_list = c(vertex_list, x, ns[i])
  }
  return(vertex_list)
}

# calculate infection rate and recovery rate
# s_neighbor_list: list of each infected patient's susceptible neighbors
# neighbor_susceptbility_list: list of each infected patient's neighbors' susceptibility values
# vaccine_efficacy_list: list of each infected patient's neighbors' vaccine efficacies (0 if not vaccinated)
# vaccine_resistances: 
# s_edge_list: list of each infected patient's edges to susceptible neighbors
# edge_weight_list: list of each infected patient's edge weights to susceptible neighbors
# I would like to set these as edge attributes, but they have to be scalars
calculate_infection_recovery_rates = function(g, infected_patients, neighbors_list) {
  infection_rates = vertex_attr(g, "infection_rate", infected_patients)
  vaccine_resistances = vertex_attr(g, "vax_resistance", infected_patients)
  infection_status = vertex_attr(g, "typ")
  s_neighbor_list = vector("list", length(infected_patients))
  neighbor_susceptibility_list = vector("list", length(infected_patients))
  vaccine_efficacy_list = vector("list", length(infected_patients))
  s_edge_list = vector("list", length(infected_patients))
  edge_weight_list = vector("list", length(infected_patients))
  for (i in 1:length(infected_patients)) {
    # with no susceptible neighbors, these are all integer(0)
    s_neighbor_list[[i]] = susceptible_neighbors(g, infected_patients[i], neighbors_list, infection_status)
    neighbor_susceptibility_list[[i]] = vertex_attr(g, "susceptibility", s_neighbor_list[[i]])
    vaccine_efficacy_list[[i]] = vertex_attr(g, "vax_efficacy", s_neighbor_list[[i]])
    vertex_list = generate_vertex_list(infected_patients[i], s_neighbor_list[[i]])
    s_edge_list[[i]] = get.edge.ids(g, vertex_list, error=TRUE)
    edge_weight_list[[i]] = edge_attr(g, "weight", s_edge_list[[i]])
  }
  # mapply is just a for loop that iterates over the elements of all given arguments
  # coefficients that modify susceptible neighbors' odds of catching disease
  # a list of vectors
  catch_coefficients = mapply(function(w, s, e, r) w*s*(1-e*(1-r)),
                              edge_weight_list, neighbor_susceptibility_list, vaccine_efficacy_list, vaccine_resistances,
                              SIMPLIFY=FALSE)
  # coefficients that modify patients' odds of spreading disease
  # a vector
  spread_coefficients = mapply(function(x,y) x*sum(y), infection_rates, catch_coefficients)
  # spread_coefs could be directly calculated with mapply(function(x,y,z) x*y%*%z, infection_rates, edge_weight_list, neighbor_susceptibility_list)
  # but we need the catch_coefs elsewhere, so it doesn't hurt to do the calculation in 2 steps
  # need to be very careful with mapply, if an element is NULL it just recycles something else
  # empty things should be integer(0), which mapply is happy to use
  # integer(0)*integer(0) = integer(0)
  # sum(integer(0)) = 0
  rateInfect = sum(spread_coefficients)
  recovery_rates = vertex_attr(g, "recovery_rate", infected_patients)
  rateReco = sum(recovery_rates)
  return(list(rateInfect=rateInfect, rateReco=rateReco,
              spread_coefficients=spread_coefficients, catch_coefficients=catch_coefficients))
}

# can't tell what these functions do but I think they have something to do with summarizing results after the simulations
act <- function(i,dfsir,df,ncl){
  if(df[i,1] == 'R'){
    indx <- 3
  }
  else if(df[i,1] == 'D'){
    indx <- 4}
  else if(df[i,1] == 'I'){
    dfsir$S[(df[i,2]+1):ncl] <- dfsir$S[(df[i,2]+1):ncl]-1
    dfsir$I[(df[i,2]+1):ncl] <- dfsir$I[(df[i,2]+1):ncl]+1
    dfsir$pI[(df[i,2]+1):ncl] <- dfsir$pI[(df[i,2]+1):ncl]+df[i,4]
    dfsir$pD[(df[i,2]+1):ncl] <- dfsir$pD[(df[i,2]+1):ncl]+df[i,5]
    return(dfsir)
  }
  else{
    dfsir$S[(df[i,2]+1):ncl] <- dfsir$S[(df[i,2]+1):ncl]-1
    return(dfsir)
  }
  dfsir$S[(df[i,2]+1):ncl] <- dfsir$S[(df[i,2]+1):ncl]-1
  dfsir$I[(df[i,2]+1):(df[i,3]+1)] <- dfsir$I[(df[i,2]+1):(df[i,3]+1)]+1
  dfsir$pI[(df[i,2]+1):(df[i,3]+1)] <- dfsir$pI[(df[i,2]+1):(df[i,3]+1)]+df[i,4]
  dfsir$pD[(df[i,2]+1):(df[i,3]+1)] <- dfsir$pD[(df[i,2]+1):(df[i,3]+1)]+df[i,5]
  dfsir[(df[i,3]+2):ncl,indx] <- dfsir[(df[i,3]+2):ncl,indx]+1
  return(dfsir)
}

# dfSIR <- function(df,N){
#   ncl <- max(c(df[,3],df[,2]),na.rm = T)+2
#   dfsir <- data.frame(S = rep(N,ncl),I = rep(0,ncl),R = rep(0,ncl), 
#                       D = rep(0,ncl))
#   
#   for(i in 1:dim(df)[1] ){dfsir <- act(i,dfsir,df,ncl)}
#   return(dfsir)  
# } 

# index, typ, t_infected, t_resolved
dfSIR = function(df, pop_size) {
  df = df[which(df$typ!="S"),c("typ","t_infected","t_resolved","p_infect","p_death")]
  df$t_infected = ceiling(df$t_infected)
  df$t_resolved = ceiling(df$t_resolved)
  ncl = max(c(df$t_infected, df$t_resolved), na.rm=TRUE)+2
  dfsir = data.frame(S=rep(pop_size, ncl), I=rep(0,ncl), R=rep(0,ncl), D=rep(0,ncl), pI=rep(0,ncl), pD=rep(0,ncl))
  for (i in 1:nrow(df)) {
    dfsir = act(i,dfsir,df,ncl)
  }
  dfsir$pI[1:(ncl-1)] <- dfsir$pI[1:(ncl-1)]/dfsir$I[1:(ncl-1)]
  dfsir$pD[1:(ncl-1)] <- dfsir$pD[1:(ncl-1)]/dfsir$I[1:(ncl-1)]
  return(dfsir)
}

# dfpreSIR <- function(df){
#   df1 <- df[which(df$typ != 'S'),c(2,3,4)]
#   df1[,2] <- ceiling(df1[,2])
#   df1[,3] <- ceiling(df1[,3])
#   return(df1)
# }


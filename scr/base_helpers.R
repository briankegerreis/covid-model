
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

# assign attributes to population
initialize_population_attrs = function(g, attr_list) {
  # for things that start the same for the whole population
  # each attr has one value
  pop_size = len(V(g))
  for (attr in names(attr_list) {
    g = set_vertex_attr(g, attr, value=rep(attr_list[[attr]],pop_size))
  }
  return(g)
}

# assign varying attributes to population
initialize_personal_attrs = function(g, attr_list) {
  # for things that vary across the population
  # each attr is the size of the population
  for (attr in names(attr_list) {
    g = set_vertex_attr(g, attr, value=attr_list[[attr]])
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
  for (attr in names(attr_list) {
    g = set_vertex_attr(g, attr, patient_zero, attr_list[[attr]])
  }
  pI = attr_list[["p_infect"]]
  rR = attr_list[["recovery_rate"]]
  g = set_vertex_attr(g, "infection_rate", patient_zero, adjustrI(pI,rR))
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
    # and apply multiplicative modifiers to target attributes
    if startsWith(attr, "modifier_") {
      target_attr = sub("modifier_", "", attr)
      g = set_vertex_attr(g, target_attr, new_case, vertex_attr(g, target_attr, new_case)*vertex_attr(g, attr, new_case))
  }
  pI = vertex_attr(g, "p_infect", new_case)
  rR = vertex_attr(g, "recovery_rate", new_case)
  g = set_vertex_attr(g, "infection_rate", new_case, adjustrI(pI,rR))
  return(g)
}

# find all susceptible neighbors of a given individual
susceptible_neighbors = function(g, x) {
  all_neighbors = igraph::neighbors(g, x)
  s_neighbors = all_neighbors[which(igraph::vertex_attr(g, "typ", all_neighbors)=="S")]
  return(as.integer(s_neighbors)) # no need to deal with the 'igraph.vs' class
}

# convert individual and vector of neighbors to c(indiv, n1, indiv, n2, indiv, n3...) to extract edges later
generate_vertex_list = function(x, ns) {
  vertex_list = integer(0)
  for (i in 1:length(ns)) {
    vertex_list = c(vertex_list, x, ns[i])
  }
  return(vertex_list)
}

# calculate infection rate and recovery rate
# s_neighbor_list: list of each infected patient's susceptible neighbors
# neighbor_susceptbility_list: list of each infected patient's neighbors' susceptibility values
# s_edge_list: list of each infected patient's edges to susceptible neighbors
# edge_weight_list: list of each infected patient's edge weights to susceptible neighbors
# I would like to set these as edge attributes, but they have to be scalars
calculate_infection_recovery_rates = function(g, infected_patients) {
  infection_rates = vertex_attr(g, "infection_rate", infected_patients)
  s_neighbor_list = vector("list", length(infected_patients))
  neighbor_susceptibility_list = vector("list", length(infected_patients))
  s_edge_list = vector("list", length(infected_patients))
  edge_weight_list = vector("list", length(infected_patients))
  for (i in 1:length(infected_patients)) {
    # with no susceptible neighbors, these are all integer(0)
    s_neighbor_list[[i]] = susceptible_neighbors(g, infected_patients[i])
    neighbor_susceptibility_list[[i]] = vertex_attr(g, "susceptibility", s_neighbor_list[[i]])
    vertex_list = generate_vertex_list(infected_patients[i], s_neighbor_list[[i]])
    s_edge_list[[i]] = get.edge.ids(g, vertex_list, error=TRUE)
    edge_weight_list[[i]] = edge_attr(g, "weight", s_edge_list[[i]])
  }
  # mapply is just a for loop that iterates over the elements of all given arguments
  # coefficients that modify susceptible neighbors' odds of catching disease
  # a list of vectors
  catch_coefficients = mapply(function(x,y) x*y, edge_weight_list, neighbor_susceptibility_list)
  # coefficients that modify patients' odds of spreading disease
  # a vector
  spread_coefficients = mapply(function(x,y), x*sum(y), infection_rates, catch_coefficients)
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


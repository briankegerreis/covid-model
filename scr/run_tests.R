rm(list=ls())
options(stringsAsFactors = FALSE)

setwd("C:/Users/brian/Documents/IU/capstone/covid/covid-model/scr")

library(igraph)
library(yaml)
source("base_sim_slim.R")
source("base_helpers.R")


# num_sim=1
# num_cores=0
# pop_size=1000
# d_connect=20
# lamb=-1
# p_death=c(0.001,0.01)
# recovery_rate=0.10
# p_infect=0.05
# t_max=Inf
# f_vulnerable=0.10
# min_cases=10
# mutation=TRUE
# mutation_params=list(mode="created", p_infect=0.10, p_death_lo=0.05, p_death_hi=0.50)
# # mutation_params=list(mode="random", p_mutate_infection=0.05, p_mutate_death=0.02,
# #                      p_increase_infection=0.75, p_increase_death=0.25,
# #                      infection_step=0.025, death_step_lo=0.001, death_step_hi=0.01,
# #                      p_corr=0.5,force_opposite_signs=FALSE,growth="linear")
# vaccination=TRUE
# vaccination_params=list(fraction_to_vaccinate=0.5, threshold_fraction=0.10,
#                         strategy="uniform", eligible_groups=c("S"), efficacy=0.95)
# pre="test"
config = read_yaml("testconfig.yml.txt")
if (config$general$t_max==0) {
  config$general$t_max = Inf
}

result = COVID_control(config$compute$num_sim, config$compute$num_cores,
                       config$general$pop_size, config$general$d_connect, config$general$lamb,
                       config$general$p_death, config$general$recovery_rate, config$general$p_infect,
                       config$general$t_max, config$general$f_vulnerable, config$general$min_cases, config$general$pre,
                       config$mutation, config$vaccination)
g = result[[1]]$graph
df = result[[1]]$df_vertices
table(df$typ)
View(df[df$typ=="R",])
table(df$typ,df$vax)

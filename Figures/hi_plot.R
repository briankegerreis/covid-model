codedir <- '/N/project/Covid/scr/'
library(ggplot2)

N <- 100
M <- 30

# Simulations folders
lfolders <- c(paste0(c('ER','PLD'),N),
              paste0(c('ER','PLD'),N,'/VAC_SIR_1_10_10'),
              paste0(c('ER','PLD'),N,'/VAC_SIR_4_10_10'))
# Simulations names 
lnames <- rep(c('Erdős–Rényi','Power-Law Distribution'),3)
lsim <- rep(c('Control ','Uniform SIR','Most popular SIR'),2*rep(1,3))

# Data frame Herd Inmiunity

dfhi <- function(dir,M){
  dfhi<-do.call(rbind,lapply(1:M, function(i) 
    read.delim(paste0(dir,'/HI/',i,'.txt'),header = F)))
  return(dfhi)
}

breaks <- c("Uniform S","Uniform SIR","Neighbor S","Neighbor SIR","Among popular S",
            "Among popular SIR","Most popular S","Most popular SIR","Control ",'Vaccination')

values <- c(RColorBrewer::brewer.pal(9,'Paired'),"#E6AB02")

df <- do.call(rbind,lapply(1:6, function(i) cbind(dfhi(lfolders[i],M),
                                                  sim =lsim[i],lvl = lnames[i] ) ) ) 
colnames(df) <- c('hi','graph','N','sim','lvl')
df$hi <- df$hi/N

plt <- ggplot(df, aes(x = hi,group = sim,fill = sim))+
  geom_histogram(aes(x = hi,group = sim, fill = sim,y = stat(count)/sum(count)),position = "identity",
                 alpha = 1, bins =40)+
  scale_fill_manual(breaks = breaks,
                    values = values)+
  labs(x = 'Proportion of infected individuals',
       y = 'Proportion of simulations',fill = 'Strategy')+
  theme_minimal(base_size=15)+
  theme(panel.border = element_rect(color = "black",fill=NA))+
  guides(color = 'none')+
  facet_grid(~lvl,shrink=F,scales="free" )

ggsave(filename = paste0(codedir,'hi_plot.png'), plot = plt,
       device = "png",width = 12, height = 7)
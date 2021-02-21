codedir <- '/N/project/Covid/scr/'
library(ggplot2)

N <- 100
M <- 30
dLD <- 14
# Simulations folders
lfolders <- c(paste0(c('ER','PLD'),N),
              paste0(c('ER','PLD'),N,'/LD_1'),
              paste0(c('ER','PLD'),N,'/LD_2'))
# Simulations names 
lnames <- rep(c('Erdős–Rényi','Power-Law Distribution'),3)
lsim <- rep(c('Control','Lockdown 1','Lockdown 2'),2*rep(1,3))

# Data frame SIR curves

dfsir <- function(dir,N,M){
  lsir <- lapply(1:M, function(i) read.delim(paste0(codedir,'/',dir,'/SIR/',i,'.txt'),header = F))
  lsir <- do.call(rbind,lapply(lsir, function(sir) cbind(t = 1:nrow(sir),sir) ))
  colnames(lsir) <- c('t','S','I','R','D')
  lsir <- cbind(aggregate(lsir$I, list(lsir$t),mean),
              aggregate(lsir$I, list(lsir$t),sd)$x)
  colnames(lsir) <- c('t','mean','sd')
  lsir$mean <-lsir$mean/N
  lsir$sd <- lsir$sd/N
  return(lsir)
}

# Lockdown start and finish time

dfLD <- function(dir,M,dLD){
  btime <- unlist(sapply(1:M,function(i) read.delim(file = paste0(codedir,'/',dir,'/tLD/',i,'.txt'),header = F)[1]))
  return(data.frame(
    bLD = mean(btime),eLD = mean(btime)+dLD,bsd = sd(btime) ))
}

# SIRD curves
df <- do.call(rbind,lapply(1:6, function(i) cbind(dfsir(lfolders[i],N,M), sim = lsim[i], lvl = lnames[i])))

# Begin and finish lockdown time 
dfb <- do.call(rbind,lapply(3:6, function(i) cbind(dfLD(lfolders[i],M,dLD), sim = lsim[i], lvl = lnames[i])))
# Mean by lvl 
dfb <- cbind(aggregate(dfb$bLD,list(dfb$lvl),mean),
      aggregate(dfb$eLD,list(dfb$lvl),mean)$x,
      aggregate(dfb$bsd,list(dfb$lvl),mean)$x)
colnames(dfb) <- c('lvl','bLD','eLD','bsd')

values = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3")

breaks = c('Control',
           'Lockdown 1',
           'Lockdown 2',
           'Lockdown')



q <- quantile(df$mean)[4]
plt <- ggplot(df,aes(x = t, y = mean,color = sim))+
  geom_line()+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color = sim),alpha = 0.3)+
  scale_color_manual(
    breaks = breaks,
    values = values
  )+
  labs(x = 'Time',y = 'Proportion of infected individuals',color = 'Strategy')+
  geom_vline(aes(xintercept = bLD,color = 'Lockdown'),linetype = 'dashed',
             data = dfb)+
  geom_vline(aes(xintercept = eLD,color = 'Lockdown'),linetype = 'dashed',
             data = dfb)+
  geom_segment(aes(y=q,yend =q,x=bLD-bsd,xend=bLD+bsd,color = 'Lockdown'),
               data = dfb)+
  geom_segment(aes(y=q,yend =q,x=eLD-bsd,xend=eLD+bsd,color = 'Lockdown'),
               data = dfb)+
  theme_minimal(base_size=15)+
  theme(panel.spacing.x=unit(10,'mm'),
        strip.text.x = element_text(size = 14),
        panel.border = element_rect(color = "black",fill=NA))+
  #       xlim(0,xli[i])+
  facet_grid(~lvl,shrink=F,
             scales="free")
ggsave(filename = paste0(codedir,'SIRD','.png'), plot = plt,
       device = "png",width = 12, height = 5)  

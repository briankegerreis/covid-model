library(ggplot2)

simdir <- '/N/project/Covid/scr'
setwd(simdir)
# Folder 

lfolder <- c('ER100',
             'ER100/LD_1',
             'ER100/LD_2')

# Simulation names

lnames <- c('Control',
            'Lockdown 1',
            'Lockdown 2')


M <- 30
N <- 100

### Degree plot

# Data frame Rank of infection vs Degree plot 

trim <- function(dir,M){
  ldf <- lapply(1:M, function(i)read.delim(paste0(dir,'/DF/',i,'.txt'),header = F,sep='\t')  )
  ldfch <- lapply(ldf, function(df)  df[!df$V2 %in% c('S','V'),] )
  ldfch <- lapply(ldfch, function(df) data.frame(dg = df$V8[order(df$V3)],n = 1:nrow(df),
                                                 tI = ceiling(df$V4[order(df$V3)])))
  for (i in 1:M) {colnames(ldfch[[i]]) <- c('dg','n','tI')}
  df <- do.call(rbind,ldfch)
  dfn <-  cbind(aggregate(df$dg, list(df$n),mean),
                aggregate(df$dg, list(df$n),sd)$x)
  colnames(dfn) <- c('n','mean','sd')
  return(dfn)
}


df <- lapply(1:length(lnames), function(i) trim(lfolder[i],M) )

for (i in 1:length(lnames)) {
  plt <- ggplot(df[[i]],aes(x = n,y = mean))+
    geom_point(alpha =0.9,shape = 1)+
    labs(x = 'Rank of infection',y = 'Average number of contacts (degree)')+
    guides(color = 'none')+
    theme_minimal(base_size=15)
  
  ggsave(filename = paste0(simdir,'/',lnames[i],'.png'), plot = plt,
         device = "png",width = 12, height = 7)  
}

#### SIRD 

# Data frame SIR curves

dfsir <- function(dir,N,M){
  lsir <- lapply(1:M, function(i) read.delim(paste0(dir,'/SIR/',i,'.txt'),header = F))
  lsir <- do.call(rbind,lapply(lsir, function(sir) cbind(t = 1:nrow(sir),sir) ))
  colnames(lsir) <- c('t','S','I','R','D')
  lsir <- cbind(aggregate(lsir$I, list(lsir$t),mean),
                aggregate(lsir$I, list(lsir$t),sd)$x)
  colnames(lsir) <- c('t','mean','sd')
  lsir$mean <-lsir$mean/N
  lsir$sd <- lsir$sd/N
  return(lsir)
}

df <- do.call(rbind,lapply(1:length(lnames), function(i) cbind(dfsir(lfolder[i],N,M), nam = lnames[i])))

plt <- ggplot(df,aes(x = t, y = mean,color = nam))+
  geom_line()+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color = nam),alpha = 0.3)+
  labs(x = 'Time',y = 'Proportion of infected individuals',color = 'Strategy')+
  theme_minimal(base_size=15)
ggsave(filename = paste0(simdir,'/SIRD','.png'), plot = plt,
       device = "png",width = 12, height = 5)  

#### Herd Inmiunity

# Data frame Herd Inmiunity

dfhi <- function(dir,M){
  dfhi<-do.call(rbind,lapply(1:M, function(i) 
    read.delim(paste0(dir,'/HI/',i,'.txt'),header = F)))
  return(dfhi)
}

df <- do.call(rbind,lapply(1:length(lnames), function(i) cbind(dfhi(lfolder[i],M),
                                                  nam = lnames[i] ) ) ) 
colnames(df) <- c('hi','graph','N','nam')
df$hi <- df$hi/N

plt <- ggplot(df, aes(x = hi,group = nam,fill = nam))+
  geom_histogram(aes(x = hi,group = nam, fill = nam,y = stat(count)/sum(count)),position = "identity",
                 alpha = 1, bins =40)+
  labs(x = 'Proportion of infected individuals',
       y = 'Proportion of simulations',fill = 'Strategy')+
  theme_minimal(base_size=15)

ggsave(filename = paste0(simdir,'hi.png'), plot = plt,
       device = "png",width = 12, height = 7)

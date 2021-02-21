codedir <- '/N/project/Covid/scr/'

# Simulations folders
lfolders <- paste0(c('ER','PLD'),N)
# Simulations names 
lnames <- c('Erdős–Rényi','Power-Law Distribution')

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

library(ggplot2)

df <-do.call(rbind, lapply(1:2, function(i) cbind(trim(lfolders[i],M),lvl = lnames[i])))

plt <- ggplot(df,aes(x = n,y = mean))+
  geom_point(aes(color = lvl),alpha =0.9,shape = 1)+
  labs(x = 'Rank of infection',y = 'Average number of contacts (degree)')+
  scale_color_manual(breaks = lnames,
                     values= c("#41AB5D","#41B6C4"))+
  facet_grid(~lvl,scales = 'free')+ 
  guides(color = 'none')+
  theme_minimal(base_size=15)+
  theme(panel.border = element_rect(color = "black",fill=NA),
        strip.text.x = element_text(size = 14))

ggsave(filename = paste0(codedir,'degree_toy.png'), plot = plt,
       device = "png",width = 12, height = 7)  

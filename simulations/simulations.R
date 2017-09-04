library(ggplot2)
library(dplyr)

data<-read.csv('simulations.csv')
options(scipen=999)
data$time<-data$time.findclades+data$time.findconflicts+data$time.sorting+data$time.bandb
data$time.ms <- round(data$time * 1000)
data$precision = data$tp/(data$tp+data$fp)
data$recall = data$tp/(data$tp+data$fn)

ggplot(data,aes(ntaxa,time,colour=factor(nchar),linetype=factor(nchar)))+
  stat_summary(fun.y=median,geom='point')+
  stat_summary(fun.y=median,geom='line')+
  scale_x_continuous(name='Number of Taxa')+
  scale_y_continuous(name='Median Time (sec)')+
  scale_linetype_discrete(name='Number of Characters')+
  scale_colour_discrete(name='Number of Characters')

ggplot(data,aes(ntaxa,time/60,colour=factor(nchar),linetype=factor(nchar)))+
  stat_summary(fun.y=max,geom='point')+
  stat_summary(fun.y=max,geom='line')+
  coord_cartesian(ylim=c(0,30))+
  scale_x_continuous(name='Number of Taxa')+
  scale_y_continuous(name='Max Time (min)')+
  scale_linetype_discrete(name='Number of Characters')+
  scale_colour_discrete(name='Number of Characters')

ggplot(subset(data,ntaxa>=17),aes(possible_clade_num,log(time)))+
  stat_summary(fun.y=median,geom='line')#+
  #facet_grid(~ntaxa)

library(MASS)
mod <- lm(data=data,log(time)~possible_clade_num+nchar*ntaxa)
bestdown <- stepAIC(mod)


lowerquart<-function(x){quantile(x)[2]}
upperquart<-function(x){quantile(x)[4]}

group_by(data,ntaxa,nchar)%>%summarize(pmin=lowerquart(precision),
                                       pmedian=median(precision),
                                       pmean=mean(precision),
                                       pmax=upperquart(precision),
                                       rmin=lowerquart(recall),
                                       rmedian=median(recall),
                                       rmean=mean(recall),
                                       rmax=upperquart(recall))

ggplot(data,aes(ntaxa,precision,colour=factor(nchar),linetype=factor(nchar)))+
  stat_summary(fun.y=mean,geom='point')+
  stat_summary(fun.y=mean,geom='line')+
  scale_x_continuous(name='Number of Taxa')+
  scale_y_continuous(name='Mean Precision')+
  scale_linetype_discrete(name='Number of Characters')+
  scale_colour_discrete(name='Number of Characters')

ggplot(data,aes(ntaxa,recall,colour=factor(nchar),linetype=factor(nchar)))+
  stat_summary(fun.y=mean,geom='point')+
  stat_summary(fun.y=mean,geom='line')+
  scale_x_continuous(name='Number of Taxa')+
  scale_y_continuous(name='Mean Recall')+
  scale_linetype_discrete(name='Number of Characters')+
  scale_colour_discrete(name='Number of Characters')

data$isConsensus <- factor('Not Consensus')
levels(data$isConsensus) <- c('Not Consensus','Consensus')
data$isConsensus[data$treenum=='consensus'] <- 'Consensus'

ggplot(data,aes(ntaxa,precision,colour=factor(nchar),linetype=factor(nchar)))+
  stat_summary(fun.y=mean,geom='point')+
  stat_summary(fun.y=mean,geom='line')+
  scale_x_continuous(name='Number of Taxa')+
  scale_y_continuous(name='Mean Precision')+
  scale_linetype_discrete(name='Number of Characters')+
  scale_colour_discrete(name='Number of Characters')+
  facet_wrap(~isConsensus)

ggplot(data,aes(ntaxa,recall,colour=factor(nchar),linetype=factor(nchar)))+
  stat_summary(fun.y=mean,geom='point')+
  stat_summary(fun.y=mean,geom='line')+
  scale_x_continuous(name='Number of Taxa')+
  scale_y_continuous(name='Mean Recall')+
  scale_linetype_discrete(name='Number of Characters')+
  scale_colour_discrete(name='Number of Characters')+
  facet_wrap(~isConsensus)
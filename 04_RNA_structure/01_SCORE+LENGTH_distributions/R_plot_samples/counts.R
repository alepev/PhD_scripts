#STANDARD PLOTS
#from table where 1st column is SCORE, others are SAMPLES with N counts for each score
#if 1st column is omitted (i.e. not per-score counts, but individual score observations with repeated values) 
#all samples must have same number of observations! (== enforces same sample size)

c=read.table("counts.txt")
crange=range(0,c) #min & max value in table (0 added to force it as minimum)
crows=row.names(c)

#line
plot(c$BKG,type='l',ylab='COUNTS',xlab='SCORE',ylim=crange,col='blue',axes=F)
#histogram
plot(c$BKG,type='h',ylab='COUNTS',xlab='SCORE',ylim=crange,col='blue',axes=F)

axis(1,at=1:length(crows),lab=crows)
axis(2)
box()
lines(c$LIST_A,col='red')
lines(c$LIST_B,col='green')
title(main='SCORES DISTRIBUTION')

##########################################################################################

#GGPLOT2 PLOTS
#from list where each row is a score observation followed by sample name (i.e. SCORE, SAMPLE pair)
#doesn't require same number of observations for each sample (== allows different sample size)

library(ggplot2)
c2<-read.table("counts_stack.txt", header=T) #if header already included, such as in this case

#NOTE: "values" and "ind" here are the names of column 1 & 2, respectively 
#(can have also more columns, ex. 1=SCORE + 2=LENGTH + 3=SAMPLE)
#kernel density plot - most informative, especially if comparing different sample sizes!
ggplot(c1,aes(x=values)) + geom_density(aes(group=ind,colour=ind,fill=ind),alpha=0.3)
#histogram (overlapping transparent bars; for non-overlapping bars use position='dodge')
ggplot(c2,aes(x=values)) + geom_histogram(aes(group=ind,colour=ind,fill=ind),alpha=0.3, binwidth=1,position='identity')

#want to add observations from an additional sample, or merge lists from different samples directly in R? 
add<-read.table("counts_stack_add-on.txt")
colnames(add)<-colnames(c2) #make sure columns are named the same
c2<- rbind(c2,add)



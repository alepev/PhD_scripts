library(ggplot2)

c<-read.table("allgenesArray.000.out", header=T, sep='\t')
add<-read.table("list_ALL_noLOADneg.000.out",header=T, sep='\t')
c=rbind(c,add)
#colnames(c)<- c('LIST','SCORE','LENGTH','NAME') #add columns names if not already in file

#kernel density plot
g=ggplot(c,aes(x=SCORE)) + geom_density(aes(group=LIST,colour=LIST,fill=LIST),alpha=0.3)
g=ggplot(c,aes(x=LENGTH)) + geom_density(aes(group=LIST,colour=LIST,fill=LIST),alpha=0.3)
#histogram (overlapping transparent bars; for non-overlapping bars use position='dodge')
g=ggplot(c,aes(x=SCORE)) + geom_histogram(aes(group=LIST,colour=LIST,fill=LIST),alpha=0.3, binwidth=1,position='identity')
g=ggplot(c,aes(x=LENGTH)) + geom_histogram(aes(group=LIST,colour=LIST,fill=LIST),alpha=0.3, binwidth=1,position='identity')

#PDF OUTPUT
pdf("mypdf.pdf")
g+theme(legend.position = "bottom") #if long labels name, place legend under plot
dev.off()
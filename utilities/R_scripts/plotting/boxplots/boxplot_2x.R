##2x STANDARD BOXPLOT (ALL TRANSCRIPTS)

library(ggplot2)

setwd("00_Gamm_lists/3UTR_GC/") #change as needed for each panel

#read values (each file will become one boxplot, in the given order)
c1 <- as.character('BKG_Gamm.GC.3UTR.tmp')
c2 <- as.character('UP_Gamm.GC.3UTR.tmp')

#transform into lists
c1 <- read.table(c1, sep='\t') 
l1 <- c1[,1]
c1 <- cbind(c1,rep('BKG',length(l1)))
colnames(c1) <- c("value","sample")
c2 <- read.table(c2, sep='\t') 
l2 <- c2[,1]
c2 <- cbind(c2,rep('UP',length(l2)))
colnames(c2) <- colnames(c1)

#merge lists
cfr <- rbind(c1,c2)

#find plot y axis range without outliers 
#(default: boxplot whiskers extend to the most extreme datapoints 
#within 1.5 times the interquantile range from the box;
#this can be changed setting a different "range" value in the boxplot function) 
cc1=range(boxplot(c1$value, plot=FALSE)$stats)
cc2=range(boxplot(c2$value, plot=FALSE)$stats)
ccc=c(cc1,cc2)
lower=min(ccc)
upper=max(ccc)
diff= upper-lower
ylower = max(c(0,lower-diff/100))
yupper = upper+diff/5
#(more space above to allow for p-value bars)

#BOXPLOT
p = ggplot(cfr, aes(x=sample, y=value)) + 
  geom_boxplot(outlier.shape=NA, fill=c("black","red")) +
  guides(fill=FALSE) + coord_cartesian(ylim = c(ylower, yupper)) + 
  theme(axis.title.x=element_blank())
p = p + stat_summary(fun.y=mean, geom="point", shape=23, size=2,
                     fill=c("white","black"),
                     color=c("white","black"))
p=  p+ stat_summary(geom = "crossbar", width=0.7, fatten=1.5, color=c("white",NA), fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

#add p-value bars
#NOTE: -L at the end = number of lines / asterisks (lines = 3x asterisks)
#coord.s are (x,y) == line.start, (xend,yend) == line.end
#x position refers to category, y refers to actual values 
#(so to plot correctly need to know exactly where plot y axis ends)

unit = (yupper-ylower)/100
up1a = upper + 3*unit
up1b = up1a + 2*unit
up2a = up1a + 6*unit
up2b = up2a + 2*unit

#lines position
lines_df <- structure(
  list(x = c(1,1,1.98), 
       y = c(up1a,up1b,up1a), 
       xend = c(1,1.98,1.98), 
       yend = c(up1b,up1b,up1b)), 
  .Names = c("x", "y", "xend", "yend"),
  row.names = c(NA, -18L), class = "data.frame"
)

#asterisk position
astpos_df <- structure(
  list(x = 1.5, y = up1b+unit), 
  .Names = c("x", "y"), 
  row.names = c(NA, -6L), 
  class = "data.frame"
)
 
#boxplot with significance bars
p <- p + geom_segment(data = lines_df, size = .5, aes(x=x, y=y, xend=xend, yend=yend)) 
p <- p + geom_text(data = astpos_df, aes(x=x, y=y), label="*", size = 6)   

p
##5x STANDARD BOXPLOT (ALL TRANSCRIPTS)

#library(ggplot2)

#setwd('FIGURE_S6_panels/A1_L_CDS/') #change as needed for each panel

#read values (each file will become one boxplot, in the given order)
c1 <- as.character('6UP_FC2_NEW.LEN.tmp')
c2 <- as.character('6DW_FC2_NEW.LEN.tmp')
c3 <- as.character('BKG_Bing.LEN.tmp')
c4 <- as.character('48UP_FC2_NEW.LEN.tmp')
c5 <- as.character('48DW_FC2_NEW.LEN.tmp')

#transform into lists
c1 <- read.table(c1, sep='\t') 
l1 <- c1[,1]
c1 <- cbind(c1,rep('6UP',length(l1)))
colnames(c1) <- c("value","sample")
c2 <- read.table(c2, sep='\t') 
l2 <- c2[,1]
c2 <- cbind(c2,rep('6DW',length(l2)))
colnames(c2) <- colnames(c1)
c3 <- read.table(c3, sep='\t') 
l3 <- c3[,1]
c3 <- cbind(c3,rep('BKG',length(l3)))
colnames(c3) <- colnames(c1)
c4 <- read.table(c4, sep='\t') 
l4 <- c4[,1]
c4 <- cbind(c4,rep('48UP',length(l4)))
colnames(c4) <- colnames(c1)
c5 <- read.table(c5, sep='\t') 
l5 <- c5[,1]
c5 <- cbind(c5,rep('48DW',length(l5)))
colnames(c5) <- colnames(c1)

#merge lists
cfr <- rbind(c1,c2,c3,c4,c5)

#find plot y axis range without outliers 
#(default: boxplot whiskers extend to the most extreme datapoints 
#within 1.5 times the interquantile range from the box;
#this can be changed setting a different "range" value in the boxplot function) 
cc1=range(boxplot(c1$value, plot=FALSE)$stats)
cc2=range(boxplot(c2$value, plot=FALSE)$stats)
cc3=range(boxplot(c3$value, plot=FALSE)$stats)
cc4=range(boxplot(c4$value, plot=FALSE)$stats)
cc5=range(boxplot(c5$value, plot=FALSE)$stats)
ccc=c(cc1,cc2,cc3,cc4,cc5)
lower=min(ccc)
upper=max(ccc)
diff= upper-lower
ylower = max(c(0,lower-diff/100))
yupper = upper+diff/5
#(more space above to allow for p-value bars)

#BOXPLOT
p = ggplot(cfr, aes(x=sample, y=value)) + 
  geom_boxplot(outlier.shape=NA, fill=c("red","blue","black","red","blue")) +
  guides(fill=FALSE) + coord_cartesian(ylim = c(ylower, yupper)) + 
  theme(axis.title.x=element_blank())
p = p + stat_summary(fun.y=mean, geom="point", shape=23, size=2,
                     fill=c("black","black","white","black","black"),
                     color=c("black","black","white","black","black"))
p=  p+ stat_summary(geom = "crossbar", width=0.7, fatten=1.5, color=c(NA,NA,"white",NA,NA), fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

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
  list(x = c(1,1,1.98, 2.02,2.02,2.98, 3.02,3.02,3.98, 4.02,4.02,5, 1,1,2.98, 3.02,3.02,5), 
       y = c(up1a,up1b,up1a, up1a,up1b,up1a, up1a,up1b,up1a, up1a,up1b,up1a, up2a,up2b,up2a, up2a,up2b,up2a), 
       xend = c(1,1.98,1.98, 2.02,2.98,2.98, 3.02,3.98,3.98, 4.02,5,5, 1,2.98,2.98, 3.02,5,5), 
       yend = c(up1b,up1b,up1b, up1b,up1b,up1b, up1b,up1b,up1b, up1b,up1b,up1b, up2b,up2b,up2b, up2b,up2b,up2b)), 
  .Names = c("x", "y", "xend", "yend"),
  row.names = c(NA, -18L), class = "data.frame"
)

#asterisk position
astpos_df <- structure(
  list(x = c(1.5, 2.5, 3.5, 4.5, 2, 4), y = c(rep(up1b+unit,4),rep(up2b+unit,2))), 
  .Names = c("x", "y"), 
  row.names = c(NA, -6L), 
  class = "data.frame"
)
 
#boxplot with significance bars
p <- p + geom_segment(data = lines_df, size = .5, aes(x=x, y=y, xend=xend, yend=yend)) 
p <- p + geom_text(data = astpos_df, aes(x=x, y=y), label=c(rep("*",6)), size = c(rep(6,6)))   

p
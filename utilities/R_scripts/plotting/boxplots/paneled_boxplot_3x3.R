#3x3 PANELED BOXPLOT (CDS with both/no UTRs)

#library(ggplot2)

#setwd("FIGURE_S6_panels/C1_CDS_UTR_LEN") #change as needed for each panel

#read values (each file will become one boxplot, in the given order)
a1 <- as.character('BKG_Bing.LEN.tmp')
a2 <- as.character('48UP_FC2_NEW.LEN.tmp')
a3 <- as.character('48DW_FC2_NEW.LEN.tmp')
b1 <- as.character('BKG_both.LEN.tmp')
b2 <- as.character('48UP_both.LEN.tmp')
b3 <- as.character('48DW_both.LEN.tmp')
c1 <- as.character('BKG_none.LEN.tmp')
c2 <- as.character('48UP_none.LEN.tmp')
c3 <- as.character('48DW_none.LEN.tmp')

#transform into lists
a1 <- read.table(a1, sep='\t') 
la1 <- a1[,1]
a1 <- cbind(a1,rep('BKG',length(la1)),rep('all',length(la1)))
colnames(a1) <- c("value","sample","UTRs")
a2 <- read.table(a2, sep='\t') 
la2 <- a2[,1]
a2 <- cbind(a2,rep('48UP',length(la2)),rep('all',length(la2)))
colnames(a2) <- colnames(a1)
a3 <- read.table(a3, sep='\t') 
la3 <- a3[,1]
a3 <- cbind(a3,rep('48DW',length(la3)),rep('all',length(la3)))
colnames(a3) <- colnames(a1)
b1 <- read.table(b1, sep='\t') 
lb1 <- b1[,1]
b1 <- cbind(b1,rep('BKG',length(lb1)),rep('both',length(lb1)))
colnames(b1) <- colnames(a1)
b2 <- read.table(b2, sep='\t') 
lb2 <- b2[,1]
b2 <- cbind(b2,rep('48UP',length(lb2)),rep('both',length(lb2)))
colnames(b2) <- colnames(a1)
b3 <- read.table(b3, sep='\t') 
lb3 <- b3[,1]
b3 <- cbind(b3,rep('48DW',length(lb3)),rep('both',length(lb3)))
colnames(b3) <- colnames(a1)
c1 <- read.table(c1, sep='\t') 
lc1 <- c1[,1]
c1 <- cbind(c1,rep('BKG',length(lc1)),rep('none',length(lc1)))
colnames(c1) <- colnames(a1)
c2 <- read.table(c2, sep='\t') 
lc2 <- c2[,1]
c2 <- cbind(c2,rep('48UP',length(lc2)),rep('none',length(lc2)))
colnames(c2) <- colnames(a1)
c3 <- read.table(c3, sep='\t') 
lc3 <- c3[,1]
c3 <- cbind(c3,rep('48DW',length(lc3)),rep('none',length(lc3)))
colnames(c3) <- colnames(a1)

#merge lists
cfr <- rbind(a1,a2,a3,b1,b2,b3,c1,c2,c3)

#find plot y axis range without outliers 
#(default: boxplot whiskers extend to the most extreme datapoints 
#within 1.5 times the interquantile range from the box;
#this can be changed setting a different "range" value in the boxplot function) 
aa1=range(boxplot(a1$value, plot=FALSE)$stats)
aa2=range(boxplot(a2$value, plot=FALSE)$stats)
aa3=range(boxplot(a3$value, plot=FALSE)$stats)
bb1=range(boxplot(b1$value, plot=FALSE)$stats)
bb2=range(boxplot(b2$value, plot=FALSE)$stats)
bb3=range(boxplot(b3$value, plot=FALSE)$stats)
cc1=range(boxplot(c1$value, plot=FALSE)$stats)
cc2=range(boxplot(c2$value, plot=FALSE)$stats)
cc3=range(boxplot(c3$value, plot=FALSE)$stats)
abc=c(aa1,aa2,aa3,bb1,bb2,bb3,cc1,cc2,cc3)
lower=min(abc)
upper=max(abc)
diff= upper-lower
ylower = max(c(0,lower-diff/100))
yupper = upper+diff/5
#(more space above to allow for p-value bars)

#BOXPLOT

p = ggplot(cfr, aes(x=sample, y=value)) +
  geom_boxplot(outlier.shape=NA, fill=c("black","red","blue")) +
  guides(fill=FALSE) + coord_cartesian(ylim = c(ylower, yupper)) + 
  theme(axis.title.x=element_blank())
p = p + stat_summary(fun.y=mean, geom="point", shape=23, size=2,
                     fill=c("white","black","black"),
                     color=c("white","black","black"))
p=  p+ stat_summary(geom = "crossbar", width=0.7, fatten=1.5, color=c("white",NA,NA), fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

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
  list(x = c(1,1,1.98, 2.02,2.02,3, 1,1,3), 
       y = c(up1a,up1b,up1a, up1a,up1b,up1a, up2a,up2b,up2a), 
       xend = c(1,1.98,1.98, 2.02,3,3, 1,3,3), 
       yend = c(up1b,up1b,up1b, up1b,up1b,up1b, up2b,up2b,up2b)), 
  .Names = c("x", "y", "xend", "yend"),
  row.names = c(NA, -9L), class = "data.frame"
)

#asterisk position
astpos_df <- structure(
  list(x = c(1.5, 2.5, 2), y = c(rep(up1b+unit,2),up2b+unit)), 
  .Names = c("x", "y"), 
  row.names = c(NA, -3L), 
  class = "data.frame"
)

#boxplot with significance bars
p <- p + geom_segment(data = lines_df, size = .5, aes(x=x, y=y, xend=xend, yend=yend)) 
p <- p + geom_text(data = astpos_df, aes(x=x, y=y), label=c(rep("*",3)), size = c(rep(6,3)))   

# split same plot into panels
p = p + facet_wrap(~ UTRs, ncol = 3)  

p

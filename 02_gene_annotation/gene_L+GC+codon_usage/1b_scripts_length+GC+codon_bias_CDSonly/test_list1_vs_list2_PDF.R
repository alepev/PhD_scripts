args <- commandArgs(TRUE)

f1 <- as.character(args[1]) #file name
f2 <- as.character(args[2]) #file name
d <- as.character(args[3]) #name of density plot
b <- as.character(args[4]) #name of boxplot
VALUE <- as.character(args[5]) #ex. GC%
bkg <- as.character(args[6]) #population, used to compute variance for the Z-test

list1 <- read.table(f1, sep='\t') 
l1 <- list1[,1]
list2 <- read.table(f2, sep='\t') 
l2 <- list2[,1]

m1 <- mean(l1, na.rm=TRUE)
m2 <- mean(l2, na.rm=TRUE)
me1 <- median(l1, na.rm=TRUE)
me2 <- median(l2, na.rm=TRUE)
sd1 <- sd(l1,na.rm=TRUE)
sd2 <- sd(l2,na.rm=TRUE)
iqr1 <- IQR(l1,na.rm=TRUE)
iqr2 <- IQR(l2,na.rm=TRUE)

if (bkg=='-') { 
		var1= var(l1,na.rm=TRUE)
		var2= var(l2,na.rm=TRUE)
     } else {
     	pop <- read.table(bkg, sep='\t') 
		pop <- pop[,1]
		pop.var=var(pop,na.rm=TRUE)
		var1=pop.var
		var2=pop.var}

#2-sample Z-test
zt = (m1 - m2) / (sqrt(var1/length(l1) + var2/length(l2)))
zt.pvalue <- pnorm(-abs(zt))

#2-sample T-test
tt <- t.test(l1,l2, var.equal = TRUE, conf.level = 0.95) #in theory same variance because they come from the same population

#2-sample Wilcoxon test (non-parametric)
wt <- wilcox.test(l1,l2, exact = FALSE, correct = FALSE, conf.int = FALSE, conf.level = 0.95)

if (bkg=='-') { 
	cat(m1,m2,sd1,sd2,me1,me2,iqr1,iqr2,zt.pvalue,tt$p.value,wt$p.value)
   } else {
	cat(m1,m2,sd1,sd2,me1,me2,iqr1,iqr2,zt.pvalue,tt$p.value,wt$p.value,pop.var)}

c1 <- cbind(list1,rep('list1',length(l1)))
colnames(c1) <- c('value','sample')
c2 <- cbind(list2,rep('list2',length(l2)))
colnames(c2) <- colnames(c1)
cfr <- rbind(c1,c2)

library(ggplot2)

#kernel density plot

#pdf(d)
#dp <- ggplot(cfr,aes(x=value)) + geom_density(aes(group=sample,colour=sample,fill=sample),alpha=0.3) + xlab(VALUE)
#dp + scale_fill_manual(values=c("red", "blue")) +
#	  scale_colour_manual(values=c("red", "blue4")) +
#	  geom_vline(xintercept = m1, col='red') + 
#    geom_vline(xintercept = m2, col='blue4') +
#     geom_vline(xintercept = me1, col='red', linetype = 'longdash') +
#     geom_vline(xintercept = me2, col='blue4', linetype = 'longdash') +
#     ggtitle(paste('KERNEL DENSITY PLOT\n mean (continuous line) + median (dashed line)',
#     '\nZ-test p-value (MEAN diff.):',sprintf("%1.2e", zt.pvalue),
#     '\nT-test p-value (MEAN diff.):',sprintf("%1.2e", tt$p.value),
#     '\nWilcoxon p-value (MEDIAN diff.):',sprintf("%1.2e", wt$p.value),'\n')) +
#     theme(plot.title = element_text(size=10))

#boxplot
    
pdf(b)
bp <- ggplot(cfr, aes(x=sample,y=value)) + geom_boxplot() + ylab(VALUE) +
		ggtitle(paste('Z-test p-value (MEAN diff.):',sprintf("%1.2e", zt.pvalue),
     '\nT-test p-value (MEAN diff.):',sprintf("%1.2e", tt$p.value),
     '\nWilcoxon p-value (MEDIAN diff.):',sprintf("%1.2e", wt$p.value),'\n'))
bp + theme(plot.title = element_text(size=10), axis.title.x=element_blank(), axis.text.x = element_text(colour="black",size=12), 
		axis.title.y=element_text(colour="black",size=15))
  
graphics.off() #closes all graphics devices without warnings, differently from dev.off()

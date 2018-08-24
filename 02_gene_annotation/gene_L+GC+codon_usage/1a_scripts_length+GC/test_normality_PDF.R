args <- commandArgs(TRUE)

f <- as.character(args[1]) #file name
pdf1 <- as.character(args[2]) #qqplot name
pdf2 <- as.character(args[3]) #histogram name

x <- read.table(f, sep='\t') 
x <- x[,1]
m <- mean(x,na.rm=TRUE)
s <- sd(x,na.rm=TRUE)

if (length(x) <= 5000) {
	t <- shapiro.test(x)
	testname <- 'Shapiro-Wilk p-value:'
	} else {
	t <- suppressWarnings(ks.test(x,y=pnorm,mean=m, sd=s,alternative='two.sided'))
	testname <- 'Kolmogorov-Smirnov p-value (N>5000):'}

cat(m,s,t$p.value)

nbins=100

pdf(pdf1)
qqnorm(x)
qqline(x)

pdf(pdf2)
hist(x,xlab='Values',main=c('Histogram + Normal Curve (RED)', paste(testname,sprintf("%1.2e", t$p.value))), breaks=nbins,prob=TRUE)
lines(density(x), col = "black", lwd=2)
curve(dnorm(x, mean=m, sd=s), col="red", lwd=2, add=TRUE)

graphics.off() #closes all graphics devices without warnings, differently from dev.off()

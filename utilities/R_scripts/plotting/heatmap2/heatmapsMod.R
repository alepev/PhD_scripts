PCC = 0.5		
Pvalue = 0.1
drop = 0 #if 1, specified experiments are not included in heatmaps
droplist= "drop.txt" #to remove controls instead, you can use new function with ann file (it generates a control sample list)

c <- data.matrix(read.csv("c_fake.txt",header=TRUE,row.names=1,sep="\t"))
p <- data.matrix(read.csv("p_fake.txt",header=TRUE,row.names=1,sep="\t"))
e <- data.matrix(read.csv("e_fake.txt",header=TRUE,row.names=1,sep="\t"))
ann <- read.csv("ann_fake.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE,sep="\t")
cls <- read.csv("col_fake.txt",header=TRUE,row.names=1,stringsAsFactors=FALSE,sep="\t")


######################################################
# DATA FILTERING
######################################################

c2=c
diag(c2)= 0

genes= rownames(p)
for (i in genes) {
	for (j in genes) {
		if (p[i,j]>Pvalue) {
			c2[i,j] <- 0 }
			}
		}		

if (PCC < 0) {
	filter <- apply(c2,1,function(i) any(i<=PCC))
	} else {
	filter <- apply(c2,1,function(i) any(i>=PCC))
	}

OKc= rownames(c2[filter,])
c=c[OKc,OKc]
rm(c2)

#RMA FILE MAPPING e FILTERING

rownames(e) <- m[rownames(e),1]
colnames(e) <- sub(pattern=".gz", replacement="", colnames(e), ignore.case=FALSE, fixed=TRUE) #can substitute '.gz' with whatever undesired extension in CEL names

if (drop == 1) {
	NOe <- scan(droplist, quiet=TRUE ,what="character") #NOTE: to use this you need to specify a "controls" file containing a list of CELs to exclude!
	e=e[,!(colnames(e) %in% NOe)]}
e=e[rownames(e) %in% OKc,]	#keep only genes present in c
e=e[rownames(c),,drop=F]	#sort genes (rows) as in c - in theory not necessary, but without it rows in expression heatmap are not sorted correctly

######################################################
# OUTPUT: HEATMAPS
######################################################

m= "average" #clustering method (best: "average" or "ward", the latter gives more compact clusters)

suppressPackageStartupMessages(library(gplots))

#CORRELATION HEATMAP (symmetrical color mapping of correlations in range [-1;+1])

#NEW: output with adjusted label size!

mycol <- colorpanel(n=20,low="cyan4",mid="grey95",high="magenta4")
rowv <- rev(as.dendrogram(hclust(as.dist(1-c),method=m))) #this allows optimal simmetry visualization when combined with Rowv & Colw settings

x=c
out='c.pdf'

l= max(c(20,nrow(x)*0.1))
if (l<=20) {
	ll= c(1,3)
	} else {
	ll= c(5,l-5)
	}
cl= min(5*l/nrow(x),1)
m=c(10,10)

pdf(file=out,height=l,width=l)
heatmap.2(x, symbreaks=-1.0, col=mycol, trace="none", density.info="none", scale="none", margins=m, Rowv=rowv, Colv=rev(rowv), cexRow=cl, cexCol=cl, lwid=ll, lhei=ll, dendrogram="row")
dev.off()

#EXPRESSION HEATMAP (default heatmap colors based on value range in the set of experiments)

#rowv is the same as before to preserve order between the 2 plots
colv <- as.dendrogram(hclust(as.dist(1-cor(e)),method=m))
nr=nrow(e)
nc=ncol(e)

#NEW: colors for each group of experiments!
groups=unique(ann[,2])
Ngroups=length(groups)

collist=vector()

#assign group colors from mapping file (ensures consistent colors between heatmaps generated at different moments)
for (i in 1:ncol(e)) {
	group= ann[which(rownames(ann)==colnames(e)[i]),2]
	color= cls[which(rownames(cls)==group),1]
	collist=c(collist,color)
	}

#assign group colors on the spot (deprecated)
#important: assumes order of CELs (rows) in ann same as CELs (columns) in expression file. if not, the assignment is messed up!
#cols = sample(colours(), Ngroups)
#for (i in 1:Ngroups) {
#	collist=c(collist,rep(cols[i],sum(ann[,2]==groups[i])))
#	}

#NEW: notation for control experiments!
filter <- apply(ann,1,function(i) i[1]==1)
CTRL= rownames(ann[filter,])
for (i in 1:length(colnames(e))) {
	if (colnames(e)[i] %in% CTRL) {
		colnames(e)[i]= paste(colnames(e)[i], 'C', sep = " ")}
	}

#NEW: output with adjusted label size!

x=e
out='e.pdf'

h= max(c(20,nrow(x)*0.1))
if (h<=20) {
	lh= c(1,3)
	} else {
	lh= c(5,h-5)
	}
cr= min(5*h/nrow(x),1)
w= max(c(20,ncol(x)*0.1))
if (w<=20) {
	lw= c(1,3)
	} else {
	lw= c(5,w-5)
	}
cc= min(5*h/ncol(x),1)
m=c(10,10)

pdf(file=out,height=h,width=w)
heatmap.2(x, trace="none", density.info="none", scale="row", margins=m, Rowv=rowv, Colv=colv, cexRow=cr, cexCol=cc, lwid=lw, lhei=lh, ColSideColors=collist)
dev.off()

#TEST AREA (uncomment to generate matrices of various size, then plot and see effect)

#sample matrices:
#x=matrix(sample(1:100),10,10) 
#x=matrix(sample(1:400),20,20)
#x=matrix(sample(1:2000),20,100) 
#x=matrix(sample(1:2000),100,20) 
#x=matrix(sample(1:10000),100,100)
#x=matrix(sample(1:40000),200,200) 
#x=matrix(sample(1:40000),300,300) 
#x=matrix(sample(1:5000),10,500) 
#x=matrix(sample(1:5000),500,10) 
#x=matrix(sample(1:250000),500,500) 

#test for correlation heatmaps: (only symmetrical matrices!)
#(see above for input parameters)
#pdf(file=out,height=l,width=l)
#heatmap.2(x, symbreaks=-1.0, col=mycol, trace="none", density.info="none", scale="none", margins=m, cexRow=cl, cexCol=cl, lwid=ll, lhei=ll, dendrogram="row")
#dev.off()

#test for expression heatmaps:
#(see above for input parameters)
#pdf(file=out,height=h,width=w)
#heatmap.2(x,trace="none", density.info="none", margins=m, cexRow=cr, cexCol=cc, lwid=lw, lhei=lh)
#dev.off()


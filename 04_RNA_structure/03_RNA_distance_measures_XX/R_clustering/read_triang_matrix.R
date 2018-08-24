tM <- file("fake_matrix.txt",'rt')
h <- scan("fake_headers.txt", quiet=TRUE ,what="character") 

#import LOWER TRIANGULAR MATRIX of distances
a <- readLines(tM)
m <- strsplit(a," ") #insert delimiter here; can use regex
m <- lapply(m,function(x) {
  x <- as.numeric(x)
  length(x) <- max(unlist(lapply(m,length))); 
  return(x)
})
m <- do.call(rbind,m)

#transform into regular matrix
m[is.na(m)] <- 0
na <- c(rep(0,dim(m)[1]))
m <- rbind(na,m)
na <- c(na,0)
m <- cbind(m,na)
m <- m + t(m) - diag(m[1,1],nrow=dim(m)[1],ncol=dim(m)[2])

#assign row/column names
rownames(m) <- h
colnames(m) <- h

#cluster elements based on distances

#SIMPLE DENDROGRAM
out='clustering.pdf'
lw= max(c(20,nrow(m)*0.2))
lh= max(c(10,nrow(m)*0.05))

pdf(file=out,width=lw,height=lh)
dendro <- hclust(as.dist(m),method="ward") #alternative: method="average"
plot(dendro)
dev.off()

#DENDROGRAM + DISTANCE HEATMAP
suppressPackageStartupMessages(library(gplots))
out='clustering+h.pdf'
x=m
l= max(c(20,nrow(x)*0.1))
if (l<=20) {
	ll= c(1,3)
	} else {
	ll= c(5,l-5)
	}
cl= min(5*l/nrow(x),1)
mar=c(20,20)

pdf(file=out,height=l,width=l)
rowv <- rev(as.dendrogram(hclust(as.dist(m),method="ward"))) 		#alternative: method="average"
mycol <- colorpanel(n=20,low="red1",mid="orange1",high="grey95")	#color heatmap
#mycol <- colorpanel(n=20,low="grey1",high="grey95")			#B&W heatmap
heatmap.2(x,col=mycol, trace="none", density.info="none", scale="none", margins=mar, Rowv=rowv, Colv=rev(rowv), cexRow=cl, cexCol=cl, lwid=ll, lhei=ll, dendrogram="row")
dev.off()



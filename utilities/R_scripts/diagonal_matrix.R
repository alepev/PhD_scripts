#DIAGONAL MATRIX CODE
#for symmetrical entries, and value range centered around 0

r <- read.csv(infile,sep="\t")
row.names(r) <- r[,1]
l=length(r[,1])
r <- data.matrix(r[,1:l+1])

#to remove symmetrical values in upper-right corner:
r[upper.tri(r)]=NA
#to remove the diagonal:
#diag(r)=NA
#r=r[2:l,1:l-1]

outfile <- sub(pattern=".txt", replacement=".pdf", infile, ignore.case =FALSE, fixed=FALSE)
mycol <- colorpanel(n=20,low="cyan4",mid="grey95",high="magenta4")
pdf(file=outfile)
heatmap.2(r,symbreaks=-1.0, col=mycol, cexRow=0.2+1/log2(l), cexCol=0.2+1/log2(l), margins=c(5,5), scale="none", cellnote=round(r,digits=2), notecol="black", notecex=0.2+1/log2(l), trace="none", density.info="none",Colv=F,Rowv=F)
dev.off()

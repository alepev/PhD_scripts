#NOTE: t() -> matrix transpose (cor computes correlations between columns of a matrix)

#best clustering methods (at least biological data): Ward e average 
#similar results, Ward finds more defined modules (good for my purpose) but less good when differences in size big vs. small ones
#distance measure: 1-cor works better when direction of correlation matters, otherwise 1-abs(cor) - same weight +/- correlations

infile="rma_fake.txt"
x <- read.csv(infile,header=TRUE,row.names=1,sep="\t")
x <- data.matrix(x)
#if normalized expression data not already in log2:
x=log2(x)


#heatmap correlation gene vs gene (i.e. correlation between rows, otherwise between columns == samples)
mycol <- colorpanel(n=20,low="cyan4",mid="grey95",high="magenta4")
c=cor(t(x)) 
heatmap.2(c,symbreaks=-1.0, symm=TRUE, col=mycol, scale="none", trace="none", density.info="none", hclust=function(c) hclust(c,method="average"), distfun=function(c) as.dist(1-c))

#heatmap expression, clustering based on correlation (computed in place)
heatmap.2(x, scale="none", trace="none", density.info="none", hclust=function(x) hclust(x,method="average"), distfun=function(x) as.dist(1-cor(t(x))))

#heatmap expression, clustering based on correlation (passed independently for rows and columns, i.e. imported from file)
rowv <- as.dendrogram(hclust(as.dist(1-cor(t(x))),method="average"))
colv <- as.dendrogram(hclust(as.dist(1-cor(x)),method="average"))
heatmap.2(x, trace="none", density="none",scale="none",Rowv=rowv, Colv=colv)

#from file: 
#c1 <- read.csv(infile1,sep="\t")
#c2 <- read.csv(infile2,sep="\t")
#rowv <- as.dendrogram(hclust(as.dist(1-c1),method="average"))
#colv <- as.dendrogram(hclust(as.dist(1-c2),method="average"))


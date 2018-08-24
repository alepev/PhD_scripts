
vp=read.table(file="input.txt", header=T, na.strings='-')
out="name"

vp$pval=-log10(vp$pval)
vp$FC=log2(vp$FC)

pdf(paste(out,"output.pdf", sep="."), width=5, height=5)

with(vp, plot(FC, pval, pch=20, main="Volcano plot", 
              xlim=c(min(-2.5,ceiling(min(FC,na.rm=TRUE))-1), max(2,ceiling(max(FC,na.rm=TRUE))+1)), 
              ylab="-log10 pval", xlab="log2 Fold Change"))
plim=8
s=subset(vp, pval>=plim & abs(FC)>=1)
if(length(s$motif_name_list)!=0) {
  with(s, points(FC, pval, pch=20, col='red'))
  with(s, text(FC, pval, labels=motif_name_list, cex=.5))
  }
s=subset(vp, pval<plim & abs(FC)>=1)
if(length(s$motif_name_list)!=0) {
  with(s, points(FC, pval, pch=20, col='blue'))
  with(s, text(FC, pval, labels=motif_name_list, cex=.5))
}
s=subset(vp, pval>=plim & abs(FC)<1)
if(length(s$motif_name_list)!=0) {
  with(s, points(FC, pval, pch=20, col='orange'))
  with(s, text(FC, pval, labels=motif_name_list, cex=.5))
}
s=subset(vp, pval<plim & abs(FC)<1)
if(length(s$motif_name_list)!=0) {
  with(s, text(FC, pval, labels=motif_name_list, cex=.5))
}

graphics.off()


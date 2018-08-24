#Chi-square test for goodness-of-fit to given expected frequencies for 3 categories (A,B,C)

argv <- commandArgs(TRUE)
A1 <- as.numeric(argv[1])
B1 <- as.numeric(argv[2])
C1 <- as.numeric(argv[3])
A2 <- as.numeric(argv[4])
B2 <- as.numeric(argv[5])
C2 <- as.numeric(argv[6])

TOT=A1+B1+C1
pA=A1/TOT
pB=B1/TOT
pC=C1/TOT
ChiSQ=chisq.test(c(A2,B2,C2),p=c(pA,pB,pC))
cat(ChiSQ$p.value)

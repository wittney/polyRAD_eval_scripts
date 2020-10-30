library(polyRAD)
myRADprelim <- readProcessSamMulti("wheat_BfaI_11262019_1_align.csv")
hh <- HindHe(myRADprelim)
TotDepthT <- rowSums(myRADprelim$locDepth)

hhByLoc <- as.data.frame(colMeans(hh, na.rm = TRUE))
colnames(hhByLoc)[1] <- "Hind"
hhByLoc$category <- rep(NA)
#add column name group that breaks up groups


library(ggplot2)
library(gghighlight)
ggplot(hhByLoc, aes(Hind, count)) + geom_histogram(bins = 50) +

<- length(which(hhByLoc$Hind > .12))


#estimate threshold
F <- 0.10
k <- 4
eqn_1 <- log10(((1-F)*(k -1))/k)
eqn_2 <- log10((1-F)*((2*k) -1)/(2*k))
threshold <- exp((eqn_1 + eqn_2)/2)

ggplot(hhByLoc , aes(x= Hind )) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.3, fill="darkcyan") 


rm(list=ls())
source("definedFunc.R")
library(MMGFM)
# Select the number of factors --------------------------------------------
q <- 3
qsvec<-rep(2,3)
N<-100
nvec = c(300, 200, 100) #c(100, 150, 80)#
result<-matrix(NA,N,4)
pveclist =  list('gaussian'=rep(150, 1),'poisson'=rep(50, 2),'binomial'=rep(60, 2))#list('gaussian'=c(50, 150),'poisson'=c(50),'binomial'=c(100,60)) #list('gaussian'=c(250, 400),'poisson'=c(250),'binomial'=c(300,200))#

getratio <- function(B, threshold){
  # B <- reslist1$hA[[2]]
  B_svalues <- svd(B)$d
  B_svalues <- B_svalues[B_svalues>threshold]
  ratio_fac <- B_svalues[-length(B_svalues)] / B_svalues[-1]
  nn <- length(ratio_fac)
  return(which.max(ratio_fac[-nn]))
}

chooseone <- function(x) {
  if(length(unique(x)) == length(x)){
    return(max(x))
  }
  tab <- table(x)
  maxFreq <- max(tab)
  modeValues <- tail(as.numeric(names(tab)[tab == maxFreq]), 1)
  return(modeValues)
}
threshold_vec <- c(1e-2, 1e-3)
for(i in 1:N){
  message("i = ", i)
datlist <- gendata_MM(seed = i,  nvec = nvec, pveclist =pveclist,
                      q = q,  d= 3,qs = qsvec,  rho = rep(3,length(pveclist)), rho_z=0.5,
                      sigmavec=rep(0.5, length(pveclist)),  sigma_eps=1)
XList <- datlist$XList
max(unlist(XList))
str(XList)
ZList <- datlist$ZList
tauList <- datlist$tauList
numvarmat <- datlist$numvarmat

reslist1 <- MMGFM(XList, ZList=ZList, numvarmat, q=6, qsvec = rep(4,3), init='MSFRVI',epsELBO = 1e-20) # LFM
hqvec <- sapply(reslist1$hA, getratio, threshold=threshold_vec[1])
hq <- chooseone(hqvec)
S <- length(reslist1$hB)
reslist1 <- MMGFM(XList, ZList=ZList, numvarmat, q=hq, qsvec = rep(4,3), init='MSFRVI',epsELBO = 1e-20) # LFM
hqsvec <- rep(NA, S)
for(s in 1:S){
  qs <- sapply(reslist1$hB[[s]], getratio, threshold=threshold_vec[2]) ## Use the vote rule.
  print(qs)
  hqsvec[s] <- chooseone(qs)
}

print(hqsvec)
result[i,1]<-hq
result[i,2:4]<-hqsvec
}
colMeans(result, na.rm=TRUE)
colSD(result, na.rm=TRUE)



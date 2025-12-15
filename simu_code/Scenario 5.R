
source("definedFunc.R")
library(MMGFM)
N <- 100
pveclist = list('gaussian'=c(100, 200),'poisson'=c(100, 200),'binomial'=c(60))
beishu<-23
nvec =beishu* c(300, 200)
qsvec<-rep(2, length(nvec))
q<-3
sigma_eps<-1
methodNames <- c("MMGFM", "GFM", "MRRR", "MSFR", "MultiCOAP", "MultiCOAPparti")
n_methods <- length(methodNames)
metricList <- list(F_tr = matrix(NA,N, n_methods), 
                   H_tr = matrix(NA,N, n_methods), 
                   V_tr = matrix(NA,N, n_methods), 
                   A_tr = matrix(NA, N, n_methods),
                   B_tr = matrix(NA, N, n_methods),
                   beta_norm=matrix(NA, N, n_methods),
                   Lambda_norm=matrix(NA, N, n_methods),
                   time = matrix(NA, N, n_methods))
for(ii in seq_along(metricList)) colnames(metricList[[ii]]) <- methodNames

for(i in 1:N){
  message("i = ", i)
  datlist <- gendata_MM(seed = i,  nvec = nvec, pveclist = pveclist ,
                        q = q,  d= 3,qs = qsvec,  rho = rep(2,length(pveclist)), rho_z=0.5,
                        sigmavec=rep(0.5, length(pveclist)), n_bin=1, sigma_eps=sigma_eps)
  XList <- datlist$XList
  max(unlist(XList))
  str(XList)
  ZList <- datlist$ZList
  tauList <- datlist$tauList
  numvarmat <- datlist$numvarmat
  reslist <- MMGFM(XList, ZList=ZList, numvarmat, q=3, qsvec = qsvec, init='MSFRVI', epsELBO = 1e-20) # LFM
  metricList$time[i,1] <- reslist$time.use
  
  ##GFM
  try({
    res_gfm <- gfm_run(XList, numvarmat, q=q, dc_eps  = 1e-20)
    metricList$time[i,2] <- res_gfm$time.use
  },silent=T)
  ### MRRR
  try({
    res_mrrr <- mrrr_run(XList, ZList, numvarmat, q, truncflag=TRUE, trunc=500)
    metricList$time[i,3] <- res_mrrr$time.use
  },silent=T)
  ### MSFR
  try({
    res_msfr <- MSFR_run(XList, ZList, numvarmat, q, qs=qsvec, maxIter=100, log.transform=TRUE)
    metricList$time[i,4] <- res_msfr$time.use
  },silent=T)
  
  ### MultiCOAP
  cc <- length(pveclist)
  dd <- max(sapply(pveclist, length))
  numvarmat1 <- matrix(0,1,length(unlist(pveclist)))
  iij<-1
  for(ic in 1:cc){
    for(id in 1:length(pveclist[[ic]]) ){
      numvarmat1[1,iij] <- unlist(pveclist[[ic]][id])
      iij<-iij+1
    }
  } 
types <- names(pveclist)
  row.names(numvarmat1) <-types[2] 
  XList1<-list()
  Xtmplist<-list()
  for(s in 1: length(nvec)){
    Xtmplist[[1]]<-cbind(XList[[s]][[2]],XList[[s]][[2]],XList[[s]][[3]])
    XList1[[s]]<-Xtmplist
  }
  res_mcoap <- multicoap_run(XcList=XList1, ZList,sum(numvarmat1), q, qsvec,epsELBO = 1e-20)
  metricList$time[i,5] <- res_mcoap$time.use
  
  ### MultiCOAPpartial
  numvarmat2 <- matrix(0,1,length(unlist(pveclist)[3:5]))
  iij<-1
  for(ic in 2:cc){
    for(id in 1:length(pveclist[[ic]]) ){
      numvarmat2[1,iij] <- unlist(pveclist[[ic]][id])
      iij<-iij+1
    }
  } 
  types <- names(pveclist)
  row.names(numvarmat1) <-types[2] 
  XList2<-list()
  Xtmplist<-list()
  for(s in 1: length(nvec)){
    Xtmplist[[1]]<-cbind(XList[[s]][[2]],XList[[s]][[3]])
    XList2[[s]]<-Xtmplist
  }
  res_mcoap <- multicoap_run(XcList=XList2, ZList,sum(numvarmat2), q, qsvec,epsELBO = 1e-20)
  metricList$time[i,6] <- res_mcoap$time.use
  
}

sapply(metricList, colMeans, na.rm=TRUE)
sapply(metricList, colSD, na.rm=TRUE)

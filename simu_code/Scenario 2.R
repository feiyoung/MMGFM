
source("definedFunc.R")
library(MMGFM)
N <- 100
q <- 3
qsvec  <- rep(2,3)
sigma_eps <- 1
methodNames <- c("MMGFM", "GFM", "MRRR", "MSFR")
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
  datlist <- gendata_MM(seed = i,  nvec = c(300, 200, 100),  
                        pveclist = list('gaussian'=c(100, 200),'poisson'=c(50, 150, 200)),
                        q = q,  d= 3,qs = qsvec,  rho = c(3, 2), rho_z=0.5,
                        sigmavec=c(0, 0),  sigma_eps=sigma_eps)
  XList <- datlist$XList
  max(unlist(XList))
  str(XList)
  ZList <- datlist$ZList
  tauList <- datlist$tauList
  numvarmat <- datlist$numvarmat
  reslist <- MMGFM(XList, ZList=ZList, numvarmat, q=3, qsvec = qsvec, init='MSFRVI', epsELBO = 1e-20) # LFM
  
  metricList$F_tr[i,1] <- meanTr(reslist$hF, datlist$F0List)
  metricList$H_tr[i,1] <-meanTr(reslist$hH, datlist$H0List)
  metricList$V_tr[i,1] <-meanTr(lapply(reslist$hv, function(x) Reduce(cbind,x) ), datlist$VList)
  metricList$A_tr[i,1] <-metric_mean(AList=reslist$hA, datlist$A0List, align='unaligned', numvarmat = numvarmat)
  metricList$B_tr[i,1] <- mean(ms_metric_mean(reslist$hB, datlist$B0List, align='unaligned', numvarmat = numvarmat))
  metricList$beta_norm[i,1] <-normvec(Reduce(cbind, reslist$hbeta)- Reduce(cbind,datlist$betaList))
  metricList$Lambda_norm[i,1] <-normvec(1/unlist(reslist$hinvLambda) - sigma_eps)
  metricList$time[i,1] <- reslist$time.use
  
 ##GFM
  try({
  res_gfm <- gfm_run(XList, numvarmat, q=q, dc_eps  = 1e-20)
  metricList$F_tr[i,2] <- meanTr(res_gfm$hF, datlist$F0List)
  metricList$A_tr[i,2] <-metric_mean(AList=res_gfm$hA, datlist$A0List, align='unaligned', numvarmat = numvarmat)
  metricList$time[i,2] <- res_gfm$time.use
  },silent=T)
  ### MRRR
  try({
  res_mrrr <- mrrr_run(XList, ZList, numvarmat, q, truncflag=TRUE, trunc=500)
  metricList$F_tr[i,3] <- meanTr(res_mrrr$hF, datlist$F0List)
  metricList$A_tr[i,3] <- metric_mean(AList=res_mrrr$hA, datlist$A0List, align='unaligned', numvarmat = numvarmat)
  metricList$beta_norm[i,3] <-normvec(res_mrrr$hbeta - Reduce(cbind,datlist$betaList))
  metricList$time[i,3] <- res_mrrr$time.use
},silent=T)
  ### MSFR
  try({
  res_msfr <- MSFR_run(XList, ZList, numvarmat, q, qs=qsvec, maxIter=100, log.transform=TRUE)
  metricList$F_tr[i,4] <- meanTr(res_msfr$hF, datlist$F0List)
  metricList$H_tr[i,4] <- meanTr(res_msfr$hH, datlist$H0List)
  metricList$A_tr[i,4] <- metric_mean(AList=res_msfr$hA, datlist$A0List, align='unaligned', numvarmat = numvarmat)
  metricList$B_tr[i,4] <- mean(ms_metric_mean(res_msfr$hB, datlist$B0List, align='unaligned', numvarmat = numvarmat))
  metricList$beta_norm[i,4] <- normvec(t(res_msfr$hbeta)- Reduce(cbind,datlist$betaList))
  metricList$Lambda_norm[i,4] <-normvec(unlist(res_msfr$hLambda) - sigma_eps)
  metricList$time[i,4] <- res_msfr$time.use
  },silent=T)

}

sapply(metricList, colMeans, na.rm=TRUE)
sapply(metricList, colSD, na.rm=TRUE)

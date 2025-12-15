


# Our method related functions --------------------------------------------

getratio <- function(B, threshold){
  # B <- reslist1$hA[[2]]
  B_svalues <- svd(B)$d
  B_svalues <- B_svalues[B_svalues>threshold]
  ratio_fac <- B_svalues[-length(B_svalues)] / B_svalues[-1]
  nn <- length(ratio_fac)
  return(which.max(ratio_fac[-nn]))
}

chooseone <- function(x) {
  # 使用table()函数统计每个元素出现的次数
  if(length(unique(x)) == length(x)){
    return(max(x))
  }
  tab <- table(x)
  maxFreq <- max(tab)
  modeValues <- tail(as.numeric(names(tab)[tab == maxFreq]), 1)
  return(modeValues)
}

# Metric functions --------------------------------------------------------

colSD <- function(mat, na.rm=TRUE){
  apply(mat, 2, sd, na.rm=na.rm)
}
normvec <- function(vec, norm=c('L1', 'L2')){

  norm <- match.arg(norm)
  if(norm=='L1'){
    val = mean(abs(vec))
  }
  if(norm=='L2'){
    val = mean(sum(vec^2))
  }
  return(val)
}

meanTr <- function(hBlist, Blist,  type='trace_statistic'){
  ###It is noted that the trace statistics is not symmetric, the true value must be in last
  trvec <- sapply(1:length(Blist), function(j) measurefun(hBlist[[j]], Blist[[j]], type = type))
  return(mean(trvec))
}
allTr <- function(hBlist, Blist,  type='trace_statistic'){
  trvec <- sapply(1:length(Blist), function(j) measurefun(Blist[[j]], hBlist[[j]], type = type))
  return(trvec)
}

metric_mean <- function(AList, A0List, type='trace_statistic', align=c('aligned', 'unaligned'), numvarmat=NULL){

  align <- match.arg(align)
  if(align == 'aligned'){
    y <- meanTr(AList, A0List, type=type)
  }else if(align == 'unaligned'){
    cc <- nrow(numvarmat)
    dd <- ncol(numvarmat)
    AtList <- list()
    m <- 1
    for(ic in 1:cc){
      pvec <- rep(0, dd+1)
      for(id in 1:dd){
        if(numvarmat[ic, id]>0){
          pvec[id+1] <- pvec[id] + numvarmat[ic, id]
          AtList[[m]] <- A0List[[ic]][(pvec[id]+1):pvec[id+1],]
          m <- m+1
        }
      }
    }
    y <- meanTr(AList, AtList, type=type)
  }
  return(y)
}
ms_metric_mean <- function(BList, B0List, type='trace_statistic', align=c('aligned', 'unaligned'), numvarmat=NULL){

  align <- match.arg(align)
  S <- length(BList)
  mvec <- rep(NA, S)
  for(s in 1:S){
    mvec[s] <- metric_mean(AList=BList[[s]], B0List[[s]], align=align, numvarmat = numvarmat)

  }
  return(mvec)
}




# Unified function -------------------------------------------------------
library(GFM)
Diag <- function(vec){
  q <- length(vec)
  if(q > 1){
    y <- diag(vec)
  }else{
    y <- matrix(vec, 1,1)
  }
  return(y)
}



mat2list <- function(B, pvec, by_row=TRUE){
  Blist <- list()
  pcum = c(0, cumsum(pvec))
  for(i in 1:length(pvec)){
    if(by_row){
      Blist[[i]] <- B[(pcum[i]+1):pcum[i+1],]
    }else{
      Blist[[i]] <- B[, (pcum[i]+1):pcum[i+1]]
    }
  }
  return(Blist)
}
vec2list <- function(y_int, nvec){

  if(length(y_int) != sum(nvec)) stop("vec2list: Check the argument: nvec!")
  yList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){

    yList_int[[i]] <- y_int[istart: sum(nvec[1:i])]
    istart <- istart + nvec[i]
  }
  return(yList_int)
}

mat2list <- function(B, pvec, by_row=TRUE){
  Blist <- list()
  pcum = c(0, cumsum(pvec))
  for(i in 1:length(pvec)){
    if(by_row){
      Blist[[i]] <- B[(pcum[i]+1):pcum[i+1],]
    }else{
      Blist[[i]] <- B[, (pcum[i]+1):pcum[i+1]]
    }
  }
  return(Blist)
}
vec2list <- function(y_int, nvec){

  if(length(y_int) != sum(nvec)) stop("vec2list: Check the argument: nvec!")
  yList_int <- list()
  istart <- 1
  for(i in 1:length(nvec)){

    yList_int[[i]] <- y_int[istart: sum(nvec[1:i])]
    istart <- istart + nvec[i]
  }
  return(yList_int)
}


gendata_Mgauss <- function (seed = 1, nvec = c(300, 200), pvec = c(50, 150), d = 3,
                           q = 6, qs = rep(2, length(nvec)),  rho = rep(1,length(pvec)), rho_z=1,
                           sigmavec=rep(0.5,length(pvec)), sigma_eps=1){

  # seed = 1; nvec = c(300, 200); pvec = c(50, 150); d = 3;
  # q = 6; qs = rep(2, length(nvec));  rho = rep(1,length(pvec)); rho_z=1;
  # sigmavec=rep(0.5,length(pvec)); sigma_eps=1
  ## rho corresponds to the parameters of signal strength for each modality
  if(length(rho) != length(pvec)) stop("gendata_Mgauss: the length of rho must be equal to that of pvec")
  if(length(sigmavec) != length(pvec)) stop("gendata_Mgauss: the length of sigmavec must be equal to that of pvec")

  S <- length(nvec)
  library(MASS)
  Diag <- GFM:::Diag
  cor.mat <- GFM:::cor.mat
  M <- length(pvec)
  p <- sum(pvec)
  set.seed(1)
  betamat <- matrix(rnorm(d*p, sd=2), d, p)
  glist <- list()
  pcum <- c(0, cumsum(pvec))


  A0List <- list()
  B0List <- list()
  Vlist <- list()
  for(s in 1:S){
    B0tmplist <- list()
    V0 <- matrix(NA, nvec[s], M)
    for(m in 1:M){
      glist[[m]] <- (pcum[m]+1):pcum[m+1]
      if(s == 1){
        Z1 <- matrix(rnorm(pvec[m] * (q+qs[s])), pvec[m], (q+qs[s]))
        svdZ <- svd(Z1)
        AB <-  svdZ$u %*% diag(sqrt(svdZ$d[1:(q+qs[s])])) * rho[m]
        AB <- AB %*% Diag(sign(AB[1,]))
        A0List[[m]] <- AB[,1:q]
        B0tmplist[[m]] <- AB[,(q+1):(q+qs[s])]
      }else{
        Z1 <- matrix(rnorm(pvec[m] * (qs[s])), pvec[m], (qs[s]))
        svdZ <- svd(Z1)
        B <-  svdZ$u %*% diag(sqrt(svdZ$d[1:qs[s]])) * rho[m]
        B <- B %*% Diag(sign(B[1,]))
        B0tmplist[[m]] <- B
      }

      V0[,m] <- rnorm(nvec[s], sd=sqrt(sigmavec[m]))
    }
    Vlist[[s]] <- V0
    B0List[[s]] <- B0tmplist
  }
  XList <- list()
  F0List <- list()
  H0List <- list()
  ZList <- list()
  tauList <- list()
  for(s in 1:S){
    Vmat <- NULL
    for(m in 1:M){
      Vmat <- cbind(Vmat, matrix(Vlist[[s]][,m], nvec[s], pvec[m]))
    }

    set.seed(seed+s)
    F0 <- matrix(rnorm(nvec[s]*q), nvec[s], q)
    F0List[[s]] <- F0
    set.seed(seed + s+100)
    H0 <- matrix(rnorm(nvec[s]*qs), nvec[s], qs)
    H0List[[s]] <- H0
    Z <- cbind(1, matrix(runif(nvec[s]*(d-1), -3, 3), nvec[s], d-1)*rho_z)
    ZList[[s]] <- Z

    B0 <- Reduce(rbind, B0List[[s]])
    A0 <- Reduce(rbind, A0List)
    X <- Z%*% betamat + F0 %*% t(A0) + H0 %*% t(B0) + Vmat +  mvrnorm(nvec[s],rep(0, p), sigma_eps*diag(p))
    XList[[s]] <- list(X)
    tautmpList <- list(matrix(0, nvec[s], length(pvec)))
    tauList[[s]] <- tautmpList
  }


  return(list(XList = XList, ZList=ZList, tauList=tauList, A0List = A0List, B0List = B0List,
              VList=Vlist, F0List=F0List, H0List = H0List,  beta= betamat,
              sigma_eps=sigma_eps,numvarmat= matrix(pvec, nrow=1)))
}


gendata_Mpois <- function (seed = 1, nvec = c(300, 200), pvec = c(50, 150), d = 3,
                            q = 6, qs = rep(2, length(nvec)),  rho = rep(1,length(pvec)), rho_z=1,
                            sigmavec=rep(0.5,length(pvec)), sigma_eps=1){

  # seed = 1; nvec = c(300, 200); pvec = c(50, 150); d = 3;
  # q = 6; qs = rep(2, length(nvec));  rho = rep(1,length(pvec)); rho_z=1;
  # sigmavec=rep(0.5,length(pvec)); sigma_eps = 1
  ## rho corresponds to the parameters of signal strength for each modality
  if(length(rho) != length(pvec)) stop("gendata_Mgauss: the length of rho must be equal to that of pvec")
  if(length(sigmavec) != length(pvec)) stop("gendata_Mgauss: the length of sigmavec must be equal to that of pvec")

  S <- length(nvec)
  library(MASS)
  Diag <- GFM:::Diag
  cor.mat <- GFM:::cor.mat
  M <- length(pvec)
  p <- sum(pvec)
  set.seed(1)
  betamat <- matrix(rnorm(d*p, sd=2), d, p)
  glist <- list()
  pcum <- c(0, cumsum(pvec))


  A0List <- list()
  B0List <- list()
  Vlist <- list()
  for(s in 1:S){
    B0tmplist <- list()
    V0 <- matrix(NA, nvec[s], M)
    for(m in 1:M){
      glist[[m]] <- (pcum[m]+1):pcum[m+1]
      if(s == 1){
        Z1 <- matrix(rnorm(pvec[m] * (q+qs[s])), pvec[m], (q+qs[s]))
        svdZ <- svd(Z1)
        AB <-  svdZ$u %*% diag(sqrt(svdZ$d[1:(q+qs[s])])) * rho[m]
        AB <- AB %*% Diag(sign(AB[1,]))
        A0List[[m]] <- AB[,1:q]
        B0tmplist[[m]] <- AB[,(q+1):(q+qs[s])]
      }else{
        Z1 <- matrix(rnorm(pvec[m] * (qs[s])), pvec[m], (qs[s]))
        svdZ <- svd(Z1)
        B <-  svdZ$u %*% diag(sqrt(svdZ$d[1:qs[s]])) * rho[m]
        B <- B %*% Diag(sign(B[1,]))
        B0tmplist[[m]] <- B
      }

      V0[,m] <- rnorm(nvec[s], sd=sqrt(sigmavec[m]))
    }
    Vlist[[s]] <- V0
    B0List[[s]] <- B0tmplist
  }
  XList <- list()
  F0List <- list()
  H0List <- list()
  ZList <- list()
  tauList <- list()
  for(s in 1:S){
    Vmat <- NULL
    for(m in 1:M){
      Vmat <- cbind(Vmat, matrix(Vlist[[s]][,m], nvec[s], pvec[m]))
    }

    set.seed(seed+s)
    F0 <- matrix(rnorm(nvec[s]*q), nvec[s], q)
    F0List[[s]] <- F0
    set.seed(seed + s+100)
    H0 <- matrix(rnorm(nvec[s]*qs), nvec[s], qs)
    H0List[[s]] <- H0
    Z <- cbind(1, matrix(runif(nvec[s]*(d-1), -3, 3), nvec[s], d-1)*rho_z)
    ZList[[s]] <- Z

    B0 <- Reduce(rbind, B0List[[s]])
    A0 <- Reduce(rbind, A0List)
    Mu <- Z%*% betamat + F0 %*% t(A0) + H0 %*% t(B0) + Vmat +  mvrnorm(nvec[s],rep(0, p), sigma_eps*diag(p))
    XList[[s]] <- list(matrix(rpois(nvec[s] * p, lambda = exp(Mu)), nvec[s], p))
    tautmpList <- list(matrix(0, nvec[s], length(pvec)))
    tauList[[s]] <- tautmpList
  }


  return(list(XList = XList, ZList=ZList, tauList=tauList, A0List = A0List, B0List = B0List,
              VList=Vlist, F0List=F0List, H0List = H0List,  beta= betamat,
              sigma_eps=sigma_eps,numvarmat= matrix(pvec, nrow=1)))
}

gendata_Mbino <- function (seed = 1, nvec = c(300, 200), pvec = c(50, 150), d = 3,
                           q = 6, qs = rep(2, length(nvec)),  rho = rep(1,length(pvec)), rho_z=1,
                           sigmavec=rep(0.5,length(pvec)), sigma_eps=1, n_bin=1){

  # seed = 1; nvec = c(300, 200); pvec = c(50, 150); d = 3;
  # q = 6; qs = rep(2, length(nvec));  rho = rep(1,length(pvec)); rho_z=1;
  # sigmavec=rep(0.5,length(pvec)); sigma_eps = 1
  ## rho corresponds to the parameters of signal strength for each modality
  if(length(rho) != length(pvec)) stop("gendata_Mgauss: the length of rho must be equal to that of pvec")
  if(length(sigmavec) != length(pvec)) stop("gendata_Mgauss: the length of sigmavec must be equal to that of pvec")

  S <- length(nvec)
  library(MASS)
  Diag <- GFM:::Diag
  cor.mat <- GFM:::cor.mat
  M <- length(pvec)
  p <- sum(pvec)
  set.seed(1)
  betamat <- matrix(rnorm(d*p, sd=2), d, p)
  glist <- list()
  pcum <- c(0, cumsum(pvec))


  A0List <- list()
  B0List <- list()
  Vlist <- list()
  for(s in 1:S){
    B0tmplist <- list()
    V0 <- matrix(NA, nvec[s], M)
    for(m in 1:M){
      glist[[m]] <- (pcum[m]+1):pcum[m+1]
      if(s == 1){
        Z1 <- matrix(rnorm(pvec[m] * (q+qs[s])), pvec[m], (q+qs[s]))
        svdZ <- svd(Z1)
        AB <-  svdZ$u %*% diag(sqrt(svdZ$d[1:(q+qs[s])])) * rho[m]
        AB <- AB %*% Diag(sign(AB[1,]))
        A0List[[m]] <- AB[,1:q]
        B0tmplist[[m]] <- AB[,(q+1):(q+qs[s])]
      }else{
        Z1 <- matrix(rnorm(pvec[m] * (qs[s])), pvec[m], (qs[s]))
        svdZ <- svd(Z1)
        B <-  svdZ$u %*% diag(sqrt(svdZ$d[1:qs[s]])) * rho[m]
        B <- B %*% Diag(sign(B[1,]))
        B0tmplist[[m]] <- B
      }

      V0[,m] <- rnorm(nvec[s], sd=sqrt(sigmavec[m]))
    }
    Vlist[[s]] <- V0
    B0List[[s]] <- B0tmplist
  }
  XList <- list()
  F0List <- list()
  H0List <- list()
  ZList <- list()
  tauList <- list()
  for(s in 1:S){
    Vmat <- NULL
    for(m in 1:M){
      Vmat <- cbind(Vmat, matrix(Vlist[[s]][,m], nvec[s], pvec[m]))
    }

    set.seed(seed+s)
    F0 <- matrix(rnorm(nvec[s]*q), nvec[s], q)
    F0List[[s]] <- F0
    set.seed(seed + s+100)
    H0 <- matrix(rnorm(nvec[s]*qs), nvec[s], qs)
    H0List[[s]] <- H0
    Z <- cbind(1, matrix(runif(nvec[s]*(d-1), -3, 3), nvec[s], d-1)*rho_z)
    ZList[[s]] <- Z

    B0 <- Reduce(rbind, B0List[[s]])
    A0 <- Reduce(rbind, A0List)
    Mu <- Z%*% betamat + F0 %*% t(A0) + H0 %*% t(B0) + Vmat +  mvrnorm(nvec[s],rep(0, p), sigma_eps*diag(p))
    XList[[s]] <- list(matrix(rbinom(nvec[s] * p, size=n_bin,prob  = 1/(1+exp(-Mu))), nvec[s], p))
    tautmpList <- list(matrix(0, nvec[s], length(pvec)))
    tauList[[s]] <- tautmpList
  }


  return(list(XList = XList, ZList=ZList, tauList=tauList, A0List = A0List, B0List = B0List,
              VList=Vlist, F0List=F0List, H0List = H0List,  beta= betamat,
              sigma_eps=sigma_eps,numvarmat= matrix(pvec, nrow=1)))
}


gendata_MM <- function (seed = 1,  nvec = c(300, 200),
                        pveclist = list('gaussian'=c(50, 150),'poisson'=c(50),'binomial'=c(100,60)),
                     q = 6,  d= 3,qs = rep(2, length(nvec)),  rho = rep(1,length(pveclist)), rho_z=1,
                     sigmavec=rep(0.5, length(pveclist)), n_bin=1, sigma_eps=1, heter_error=FALSE){

  # seed = 1; nvec = c(300, 200); pveclist = list('gaussian'=c(50, 150),'poisson'=c(50),'binomial'=c(100,60)); d = 3;
  # q = 6; qs = rep(2, length(nvec));  rho = rep(1,length(pveclist)); rho_z=1;
  # sigmavec=rep(0.5,length(pveclist)); sigma_eps = 1;n_bin=1
  ## rho corresponds to the parameters of signal strength for each modality
  if(length(rho) != length(pveclist)) stop("gendata_MM: the length of rho must be equal to that of pveclist")
  if(length(sigmavec) != length(pveclist)) stop("gendata_MM: the length of sigmavec must be equal to that of pveclist")

  types <- names(pveclist)
  S <- length(nvec)
  library(MASS)
  Diag <- GFM:::Diag
  cor.mat <- GFM:::cor.mat
  cc <- length(pveclist)
  dd <- max(sapply(pveclist, length))
  pvec <- sapply(pveclist, sum)
  numvarmat <- matrix(0, cc, dd)
  betaList <- list()
  for(ic in 1:cc){
    ptmp <- pvec[ic]
    for(id in 1:length(pveclist[[ic]]) ){
      numvarmat[ic, id] <- pveclist[[ic]][id]
    }
    set.seed(1)
    betaList[[ic]] <- matrix(rnorm(d*ptmp, sd=2), d, ptmp)
  }
  row.names(numvarmat) <- types
  pcum <- c(0, cumsum(pvec))

  A0List <- list()
  B0List <- list()
  Vlist <- list()
  M <-  sum(numvarmat>0)
  for(s in 1:S){
    B0tmplist <- list()
    V0 <- matrix(NA, nvec[s], M)
    m <- 1
    for(ic in 1:cc){
      #
      Amat_tmp <- NULL
      Bmat_tmp <- NULL
      for(id in 1:length(pveclist[[ic]])){
        p_cd <- numvarmat[ic, id]
        if(s == 1){
          Z1 <- matrix(rnorm(p_cd * (q+qs[s])), p_cd, (q+qs[s]))
          svdZ <- svd(Z1)
          AB <-  svdZ$u %*% diag(sqrt(svdZ$d[1:(q+qs[s])])) * rho[ic]
          AB <- AB %*% Diag(sign(AB[1,]))
          Amat_tmp <- rbind(AB[,1:q], Amat_tmp)
          Bmat_tmp <- rbind(AB[,(q+1):(q+qs[s])], Bmat_tmp)
          # B0tmplist[[m]] <- AB[,(q+1):(q+qs[s])]
        }else{
          Z1 <- matrix(rnorm(p_cd * (qs[s])), p_cd, (qs[s]))
          svdZ <- svd(Z1)
          B <-  svdZ$u %*% diag(sqrt(svdZ$d[1:qs[s]])) * rho[ic]
          B <- B %*% Diag(sign(B[1,]))
          Bmat_tmp <- rbind(B, Bmat_tmp)

        }

        V0[,m] <- rnorm(nvec[s], sd=sqrt(sigmavec[ic]))
        m <- m+1
      }

      if(s == 1){
        A0List[[ic]] <- Amat_tmp
      }
      B0tmplist[[ic]] <- Bmat_tmp
    }
    Vlist[[s]] <- V0
    B0List[[s]] <- B0tmplist
  }
  XList <- list()
  F0List <- list()
  H0List <- list()
  ZList <- list()
  tauList <- list()
  for(s in 1:S){
    Vtmplist <- list()
    m <- 1
    for(ic in 1:cc){
      Vmat <- NULL
      for(id in 1:length(pveclist[[ic]])){
        p_cd <- numvarmat[ic, id]
        Vmat <- cbind(Vmat, matrix(Vlist[[s]][,m], nvec[s], p_cd))
        m <- m+1
      }
      Vtmplist[[ic]] <- Vmat
    } ## Vmat in R^(ns * p_total)

    set.seed(seed+s)
    F0 <- matrix(rnorm(nvec[s]*q), nvec[s], q)
    F0List[[s]] <- F0
    set.seed(seed + s+100)
    H0 <- matrix(rnorm(nvec[s]*qs[s]), nvec[s], qs[s])
    H0List[[s]] <- H0
    Z <- cbind(1, matrix(runif(nvec[s]*(d-1), -3, 3), nvec[s], d-1)*rho_z)
    ZList[[s]] <- Z


    Xtmplist <- list()
    tautmpList <- list()
    for(ic in 1:cc){
      B_tmp <-  B0List[[s]][[ic]]
      A_tmp <- A0List[[ic]]
      beta_tmp <- betaList[[ic]]
      p <- sum(pveclist[[ic]])
      if(heter_error){
        Lam0 <- sigma_eps*diag(runif(p, 0.5, 1.5))
      }else{
        Lam0 <- sigma_eps*diag(p)
      }
      if(types[ic]=='gaussian'){
        Mu <- Z%*% beta_tmp + F0 %*% t(A_tmp) + H0 %*% t(B_tmp) +  Vtmplist[[ic]] +  mvrnorm(nvec[s],rep(0, p), Lam0)
        Xtmplist[[ic]] <- Mu
      }else if(types[ic]=='poisson'){
        Mu <- Z%*% beta_tmp + F0 %*% t(A_tmp) + H0 %*% t(B_tmp) +  Vtmplist[[ic]] +  mvrnorm(nvec[s],rep(0, p), Lam0)
        Xtmplist[[ic]]  <- matrix(rpois(nvec[s] * p, lambda = exp(Mu)), nvec[s], p)
      }else if(types[ic]=='binomial'){
        Mu <- Z%*% beta_tmp + F0 %*% t(A_tmp) + H0 %*% t(B_tmp) +  Vtmplist[[ic]] +  mvrnorm(nvec[s],rep(0, p), Lam0)
        Xtmplist[[ic]]  <- matrix(rbinom(nvec[s] * p, size=n_bin,prob  = 1/(1+exp(-Mu))), nvec[s], p)
      }
      tautmpList[[ic]] <- matrix(0, nvec[s], length(pveclist[[ic]]))
    }
    XList[[s]] <- Xtmplist

    tauList[[s]] <- tautmpList
  }


  return(list(XList = XList, ZList=ZList, tauList=tauList, A0List = A0List, B0List = B0List,
              VList=Vlist, F0List=F0List, H0List = H0List,  betaList= betaList,
              sigma_eps=sigma_eps,numvarmat= numvarmat, types=types, Lam0 = Lam0))
}




# Compared methods --------------------------------------------------------

## Other methods
gfm_run <- function(XList, numvarmat, q, dc_eps=1e-10, ...){
  library(GFM)
  S <- length(XList)
  nvec <- sapply(XList, function(x) nrow(x[[1]]))
  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]
  XnewList <- list()
  cc <- nrow(numvarmat)
  for(ic in 1:cc){
    tmpmat <- NULL
    for(s in 1:S){
      tmpmat <- rbind(tmpmat, XList[[s]][[ic]])
    }
    XnewList[[ic]] <- tmpmat

  }
  types <- row.names(numvarmat)
  tic <- proc.time()
  # res_gfm <- gfm(XnewList, types=types, q=q, dc_eps=dc_eps)
  res_gfm <- gfm(XnewList, types=types, q=q, dc_eps=dc_eps, ...)
  # str(res_gfm)
  hFList <- mat2list(res_gfm$hH, nvec)
  hAList <- mat2list(res_gfm$hB, pvec, by_row = TRUE)
  hmuList <- vec2list(res_gfm$hmu, pvec)
  toc <- proc.time()

  rm(res_gfm)
  return(list(hF=hFList, hA=hAList, hmu=hmuList, time.use= toc[3] - tic[3]))
}
mrrr_run <- function(XList, ZList, numvarmat, q, truncflag=FALSE, trunc=500,lambdaSVD=0.1){
  mrrr_single <- function(Y, Z, numvarmat, rank0,family=list(poisson()),
                          familygroup,  epsilon = 1e-4, sv.tol = 1e-2,
                          lambdaSVD=0.1, maxIter = 2000, trace=TRUE, truncflag=FALSE, trunc=500){
    # epsilon = 1e-4; sv.tol = 1e-2; maxIter = 30; trace=TRUE,lambdaSVD=0.1

    require(rrpack)
    q <- rank0
    n <- nrow(Y); p <- ncol(Y)
    X <- cbind(cbind(1, Z),  diag(n))
    d <- ncol(Z)
    ## Trunction
    if(truncflag){
      Y[Y>trunc] <- trunc
      Y[Y< -trunc] <- -trunc
    }



    pvec <- as.vector(t(numvarmat))
    pvec <- pvec[pvec>0]
    pcums <- cumsum(pvec)
    idxlist <- list()
    idxlist[[1]] <- 1:pcums[1]
    if(length(pvec)>1){
      for(i in 2:length(pvec)){
        idxlist[[i]] <- (pcums[i-1]+1):pcums[i]
      }
    }

    svdX0d1 <- svd(X)$d[1]
    init1 = list(kappaC0 = svdX0d1 * 5)
    offset = NULL
    control = list(epsilon = epsilon, sv.tol = sv.tol, maxit = maxIter,
                   trace = trace, gammaC0 = 1.1, plot.cv = TRUE,
                   conv.obj = TRUE)
    res_mrrr <- mrrr(Y=Y, X=X[,-1], family = family, familygroup = familygroup,
                     penstr = list(penaltySVD = "rankCon", lambdaSVD = lambdaSVD),
                     control = control, init = init1, maxrank = rank0+d) #

    hmu <- res_mrrr$coef[1,]
    hbeta <- t(res_mrrr$coef[2:(d+1),])
    hTheta <- res_mrrr$coef[-c(1:(d+1)),]


    # Matrix::rankMatrix(hTheta)
    svd_Theta <- svd(hTheta, nu=q, nv=q)
    hH <- svd_Theta$u
    hB <- svd_Theta$v %*% Diag(svd_Theta$d[1:q])


    return(list(hH=hH, hB=hB, hmu= hmu, beta=hbeta))
  }

  tic <- proc.time()
  S <- length(XList)
  nvec <- sapply(XList, function(x) nrow(x[[1]]))
  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]
  XnewList <- list()
  cc <- nrow(numvarmat)
  types <- row.names(numvarmat)
  family.use <- list()
  familygroup <- NULL
  for(ic in 1:cc){
    tmpmat <- NULL
    for(s in 1:S){
      tmpmat <- rbind(tmpmat, XList[[s]][[ic]])
    }
    XnewList[[ic]] <- tmpmat
    if(types[ic]=='gaussian'){
      family.use[[ic]] <- gaussian()
    }else if(types[ic] == 'poisson'){
      family.use[[ic]] <- poisson()
    }else if(types[ic] == 'binomial'){
      family.use[[ic]] <- binomial()
    }
    familygroup <- c(familygroup,rep(ic, sum(numvarmat[ic, ])))
  }

  Y <- Reduce(cbind, XnewList)
  Z <- Reduce(rbind, ZList)
  res_mrrr <- mrrr_single(Y, Z[,-1], numvarmat, rank0=q, family = family.use, familygroup = familygroup,
                          epsilon = 1e-4, sv.tol = 1e-2,
                          lambdaSVD=0.1, maxIter = 2000, trace=TRUE, truncflag=truncflag, trunc=trunc)

  toc <- proc.time()
  time_mrrr <- toc[3] - tic[3]

  # str(res_gfm)
  hFList <- mat2list(res_mrrr$hH, nvec)
  hAList <- mat2list(res_mrrr$hB, pvec, by_row = TRUE)
  hmuList <- vec2list(res_mrrr$hmu, pvec)
  return(list(hF=hFList, hA=hAList, hbeta=t(cbind(res_mrrr$hmu, res_mrrr$beta)), hmu=hmuList, time.use=time_mrrr))

}
MSFR_run <- function(XList, ZList, numvarmat, q, qsvec, maxIter=1e4, log.transform=FALSE,load.source=FALSE, dir.source=NULL){

  # require(MSFA)
  require(psych)
  if(!load.source){
    source(paste0(dir.source, "MSFR_main_R_MSFR_V1.R"))
  }
  S <- length(XList)
  nvec <- sapply(XList, function(x) nrow(x[[1]]))
  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]
  XnewList <- list()
  cc <- nrow(numvarmat)
  types <- row.names(numvarmat)
  for(s in 1:S){
    tmpmat <- NULL
    for(ic in 1:cc){
      tmpmat2 <- XList[[s]][[ic]]
      if(types[ic]== 'poisson' && log.transform){
        tmpmat2 <- log(XList[[s]][[ic]]+1)
      }

      tmpmat <- cbind(tmpmat, tmpmat2)
    }
    XnewList[[s]] <- tmpmat

  }
  #fa <- psych::fa
  B_s <- ZList

  t1<- proc.time()
  test <- start_msfa(XnewList, B_s, 5, k=q, j_s=qsvec, constraint = "block_lower2", method = "adhoc")
  # EM_beta <- ecm_msfa(X_s, B_s, start=test, trace = FALSE, nIt=maxIter, constraint = "block_lower1")
  EM_beta <- ecm_msfa(XnewList, B_s, start=test, trace = TRUE, nIt=maxIter)
  t2<- proc.time()
  res_msfr <- list()
  res_msfr$hF <- lapply(EM_beta$E_f, t)
  res_msfr$hH <- lapply(EM_beta$E_l, t)
  res_msfr$time.use <- t2[3] - t1[3]
  res_msfr$hLambda <-EM_beta$psi_s

  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]

  res_msfr$hA <- mat2list(EM_beta$Phi, pvec, by_row = TRUE)
  hB_unalign <- EM_beta$Lambda_s
  tmpList <- list()
  for(s in 1:S){
    tmpList[[s]] <- mat2list(EM_beta$Lambda_s[[s]], pvec, by_row = TRUE)
  }
  res_msfr$hB <- tmpList
  res_msfr$hbeta <- EM_beta$beta

  rm(EM_beta)
  return(res_msfr)
}
multicoap_run <- function(XcList, ZList,numvarmat, q, qsvec,  init = c("MSFRVI", "LFM"),
                          maxIter=30, epsELBO = 1e-20){

  library(MultiCOAP)
  tic <- proc.time()
  XcList <- lapply(XcList, function(x) x[[1]])
  d <- ncol(ZList[[1]])
  res_mc <- MultiCOAP(XcList, ZList, q=q, qs=qsvec, rank_use = d, init=init, maxIter=maxIter, epsELBO=epsELBO)
  toc <- proc.time()
  time.use <- toc[3] - tic[3]


  res_mcoap <- list()
  res_mcoap$hF <- res_mc$F
  res_mcoap$hH <- res_mc$H
  res_mcoap$hLambda <-res_mc$Lambda

  S <- length(XcList)
  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]

  res_mcoap$hA <- mat2list(res_mc$A, pvec, by_row = TRUE)
  tmpList <- list()
  for(s in 1:S){
    tmpList[[s]] <- mat2list(res_mc$B[[s]], pvec, by_row = TRUE)
  }
  res_mcoap$hB <- tmpList
  res_mcoap$hbeta <- res_mc$bbeta
  res_mcoap$time.use <- time.use
  rm(res_mc)
  return(res_mcoap)
}


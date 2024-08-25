#' Generate simulated data
#' @description Generate simulated data from MMGFM models
#' @param seed a postive integer, the random seed for reproducibility of data generation process.
#' @param nvec a  vector with postive integers, specify the sample size in each study/source.
#' @param pveclist a named list, specify the number of modalities for each type and variable dimension in each type of modatlity.
#' @param d a postive integer,  specify the dimension of covariate matrix.
#' @param q a postive integer,  specify the number of study-shared factors.
#' @param qs a  vector with postive integers, specify the number of study-specified factors.
#' @param rho a numeric vector with \code{length(pveclist)} and positive elements, specify the signal strength of loading matrices for each modality type.
#' @param rho_z a positive real, specify the signal strength of covariates.
#' @param sigmavec a positive real vector with \code{length(pveclist)}, specify the variance of study-specified and modality variable-shared factors; default as 0.5 for each element.
#' @param n_bin a positive integer, specify the number of trails when generate Binomial modality matrix; default as 1.
#' @param sigma_eps a positive real, the variance of overdispersion error; default as 1.
#' @param heter_error a logical value, whether to generate the heterogeneous error; default as FALSE.
#' return a list including the following components:
#' \itemize{
#'   \item \code{hbeta} - a M-length list composed by the estimated regression coefficient matrix for each modality;
#'   \item \code{hA} - a M-length list composed by the loading matrix corresponding to study-shared factors for each modality;
#'   \item \code{hB} - a S-length list composed by a M-length loading matrix list corresponding to study-specified factors for each study;
#'   \item \code{hF} - a S-length list composed by the posterior estimation of study-shared factor matrix for each study;
#'   \item \code{hH} - a S-length list composed by the posterior estimation of study-specified factor matrix for each study;
#'   \item \code{hSigma} - a S-length list composed by the estimated posterior variance of the study-shared factor;
#'   \item \code{hPhi} - a S-length list composed by the estimated posterior variance of study-specified factor;
#'   \item \code{hv} - a S-length list composed by a M-length vector list corresponding to the posterior estimation of study-specified and modality variable-shared factor for each study and modality;
#'   \item \code{hzeta} - the estimated posterior variance for study-specified and modality variable-shared factor;
#'   \item \code{hsigma2} - the estimated variance for study-specified and modality variable-shared factor;
#'   \item \code{hinvLambda} - a S-length list composed by a M-length vector list corresponding to the inverse of the estimated variances of error;
#'   \item \code{S} - the approximated posterior covariance for each row of F;
#'   \item \code{ELBO} -  the ELBO value when algorithm stops;
#'   \item \code{ELBO_seq} - the sequence of ELBO values.
#'   \item \code{time_use} - the running time in model fitting of SpaCOAP;
#' }
#' @importFrom stats rnorm rpois rbinom
#' @importFrom MASS  mvrnorm
#' @export
#'
#' @examples
#' q <- 3; qsvec<-rep(2,3)
#' nvec <- c(100, 120, 100)
#' pveclist <-  list('gaussian'=rep(150, 1),'poisson'=rep(50, 2),'binomial'=rep(60, 2))
#' datlist <- gendata_mmgfm(seed = 1,  nvec = nvec, pveclist =pveclist,
#'                          q = q,  d= 3,qs = qsvec,  rho = rep(3,length(pveclist)), rho_z=0.5,
#'                          sigmavec=rep(0.5, length(pveclist)),  sigma_eps=1)
gendata_mmgfm <- function (seed = 1,  nvec = c(300, 200),
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
  #library(MASS)
  Diag <- function(vec){
    q <- length(vec)
    if(q > 1){
      y <- diag(vec)
    }else{
      y <- matrix(vec, 1,1)
    }
    return(y)
  }
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




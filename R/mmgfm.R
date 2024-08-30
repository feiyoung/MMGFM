
# generate man files
# devtools::document()
# R CMD check --as-cran MMGFM_1.1.tar.gz
## usethis::use_data(dat_r2_mac)
# pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_reference()
# pkgdown::build_article("COAPsimu")
# pkgdown::build_article("ProFASTdlpfc2")

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('html_document'))
# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('pdf_document'), clean = F)


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

#' @importFrom  MultiCOAP MSFRVI
initialize_MSFR <- function(XList, ZList, numvarmat,tauList=NULL,  q, qs,
                            log.transform=TRUE, maxIter=30, ...){
  #library(MultiCOAP)
  S <- length(XList)
  nvec <- sapply(XList, function(x) nrow(x[[1]]))
  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]
  d <- ncol(ZList[[1]])
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
  if (is.null(tauList)) {
    tauList <- lapply(nvec, function(n1) rep(1, n1))
  }
  # res_msfrvi <- MSFRVI(XnewList, ZList, q = q, qs = qs, rank_use = d)
  res_msfrvi <- MSFRVI(XnewList, ZList, q = q, qs = qs, rank_use = d,
                       aList = tauList,maxIter = maxIter, ...)

  res_msfr <- list()
  res_msfr$hF <- res_msfrvi$F
  res_msfr$hH <- res_msfrvi$H
  pvec <- as.vector(t(numvarmat))
  pvec <- pvec[pvec>0]

  res_msfr$hA <- res_msfrvi$A
  res_msfr$hB <- res_msfrvi$B
  res_msfr$hbeta <- res_msfrvi$bbeta
  rm(res_msfrvi)
  return(res_msfr)
}


#' Fit the high-dimensional multi-study multi-modality covariate-augmented generalized factor model
#' @description Fit the high-dimensional multi-study multi-modality covariate-augmented generalized factor model via variational inference.
#' @param XList a S-length list with each component a m-length list composed by a combined modality matrix of the same type modalities, which is the observed  matrix from each source/study and each modality, where m is the number of modality types.
#' @param ZList a S-length list with each component a matrix that is the covariate matrix from each study.
#' @param numvarmat a m-by-T matrix with rownames modality types that specifies the variable number for each modality of each modality type, where m is the number of modality types, T is the maximum number of modalities for one of modality types  .
#' @param q an optional string, specify the number of study-shared factors; default as 15.
#' @param qsvec a integer vector with length S, specify the number of study-specifed factors; default as 2.
#' @param tauList an optional S-length list with each component a m-length list correponding the offset term for each combined modality of each study; default as full-zero matrix.
#' @param init an optional string, specify the initialization method, supporting "MSFRVI", "random" and "LFM", default as "MSFRVI".
#' @param epsELBO  an optional positive vlaue, tolerance of relative variation rate of the envidence lower bound value, defualt as '1e-5'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed an optional integer, specify the random seed for reproducibility in initialization.
#' @return return a list including the following components:
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
#' @details If \code{init="MSFRVI"}, it will use the results from multi-study linear factor model in MultiCOAP package as initial values; If \code{init="LFM"}, it will use the results from linear factor model by combing data from all studies as initials.
#' @seealso None
#' @references None
#' @export
#' @useDynLib MMGFM, .registration = TRUE
#' @importFrom  irlba irlba
#' @importFrom  Rcpp evalCpp
#' @importFrom stats  coef lm residuals runif
#' @importFrom utils tail
#' @examples
#' q <- 3; qsvec<-rep(2,3)
#' nvec <- c(100, 120, 100)
#' pveclist <-  list('gaussian'=rep(150, 1),'poisson'=rep(50, 2),'binomial'=rep(60, 2))
#' datlist <- gendata_mmgfm(seed = 1,  nvec = nvec, pveclist =pveclist,
#'                          q = q,  d= 3,qs = qsvec,  rho = rep(3,length(pveclist)), rho_z=0.5,
#'                          sigmavec=rep(0.5, length(pveclist)),  sigma_eps=1)
#' XList <- datlist$XList
#' ZList <- datlist$ZList
#' numvarmat <- datlist$numvarmat
#' ### For illustration, we set maxIter=3. Set maxIter=50 when running formally
#' reslist1 <- MMGFM(XList, ZList=ZList, numvarmat, q=q, qsvec = qsvec, init='MSFRVI',maxIter = 3)
#' str(reslist1)


MMGFM <- function(XList, ZList,  numvarmat, tauList=NULL, q=15, qsvec = rep(2, length(XList)),
                  init=c("MSFRVI", "random", "LFM"), epsELBO=1e-12,
                  maxIter=30,  verbose=TRUE, seed=1){

  # epsELBO=1e-8; maxIter=30; add_IC_iter=FALSE; verbose=TRUE; seed=1
  # init='random'
  # init = 'LFM'; init = 'MSFRVI'
  init <- match.arg(init)
  Diag <- function(vec) {
    q <- length(vec)
    if (q > 1) {
      y <- diag(vec)
    }
    else {
      y <- matrix(vec, 1, 1)
    }
    return(y)
  }
  get_initials <- function(X, q) {
    #library(irlba)
    n <- nrow(X)
    p <- ncol(X)
    mu <- colMeans(X)
    X <- X - matrix(mu, nrow = n, ncol = p, byrow = TRUE)
    svdX <- irlba(A = X, nv = q)
    PCs <- sqrt(n) * svdX$u
    loadings <- svdX$v %*% Diag(svdX$d[1:q])/sqrt(n)
    dX <- PCs %*% t(loadings) - X
    Lam_vec <- colSums(dX^2)/n
    return(list(hH = PCs, hB = loadings, hmu = mu, sigma2vec = Lam_vec))
  }
  add_ident_beta <- function(ZList, betaList, vList, FList, AList, method='remove_fv'){


    S <- length(ZList)
    M <- length(vList[[1]])
    V <- NULL
    for(s in 1:S){
      tmpMat <- NULL
      for(m in 1:M){
        tmpMat <- cbind(tmpMat, vList[[s]][[m]])
      }
      V <- rbind(V,tmpMat)
    }

    mu_V <- colMeans(V)


    nvec <- sapply(ZList, nrow)
    Zmat <- Reduce(rbind, ZList)
    d <- ncol(Zmat)
    Mu_v <- NULL
    for(m in 1:M){
      Mu_v <- c(Mu_v, rep(mu_V[m], ncol(betaList[[m]])))
    }
    shift_beta <- qr.solve(t(Zmat) %*% Zmat) %*% t(Zmat) %*% matrix(Mu_v, nrow(Zmat), length(Mu_v))
    hbeta <- Reduce(cbind,betaList) + shift_beta
    if(method=='remove_v'){
      hbetaList <- mat2list(hbeta, sapply(betaList, ncol), by_row = FALSE)
    }else if(method=='remove_fv'){
      mu_f <- colMeans(Reduce(rbind, FList))
      Amat <- t(Reduce(rbind,AList))
      for(id in 1:d){
        hbeta[id,] <- hbeta[id,] + Amat[id,] * mu_f[id]
      }
      hbetaList <-mat2list(hbeta, sapply(betaList, ncol), by_row = FALSE)
      for(s in 1:S){
        FList[[s]] <- FList[[s]] - matrix(mu_f, nrow=nvec[s], ncol=length(mu_f), byrow = TRUE)
      }
    }
    names(hbetaList) <- names(betaList)
    ## re-organize v
    V_center <- V - matrix(mu_V, nrow=nrow(V), ncol=ncol(V), byrow = TRUE)
    hVList <- lapply(mat2list(V_center,  pvec=nvec), mat2list, pvec=rep(1, M), by_row=FALSE)
    for(s in 1:S){
      for(m in 1:M){
        vList[[s]][[m]][,1] <- hVList[[s]][[m]]
      }
    }

    return(list(hbeta=hbetaList, hV=vList, hF=FList))
  }

  S <- length(XList)
  ### Check input arguments
  flag <- sapply(XList, function(x) inherits(x, "list"))
  if(!all(flag)) stop("MMGFM: each component of XList must be a list consisting of each matrix for each modality type!")
  flag <- sapply(XList[[1]], function(x) inherits(x, "matrix"))
  if(!all(flag)) stop("MMGFM: each component of XList must be a list consisting of each matrix for each modality type!")

  flag <- sapply(ZList, function(x) inherits(x, "matrix"))
  if(!all(flag)) stop("MMGFM: each component of ZList must be a  matrix!")
  if(length(ZList) !=S) stop("MMGFM: length of ZList must be equal to the length of XList!")

  if(!inherits(numvarmat, "matrix") || is.null(row.names(numvarmat)))  stop("MMGFM: numvarmat must be a  matrix with row.names!")
  if(length(qsvec) !=S)  stop("MMGFM: length of qsvec must be equal to the length of ZList and XList!")

  if(q<2) stop("MMGFM: the number of study-shared factors must be greater than 1!")
  flag <- sapply(qsvec, function(x) x<1 )
  if(any(flag)) stop("MMGFM: the number of study-specified factors must be greater than 0!")


  cc <- nrow(numvarmat)
  dd <- ncol(numvarmat)
  num_modals <- rowSums(numvarmat)

  nvec <- sapply(XList, function(x) nrow(x[[1]]))
  s_max <- which.max(nvec)
  d <- ncol(ZList[[1]])
  if(is.null(tauList)){
    tauList <- list()
    for(s in 1:S){
      tmptaulist <- list()
      for(ic in 1:cc){
        tmptaulist[[ic]] <- matrix(0, nrow=nvec[s], ncol=sum(numvarmat[ic,]>0))
      }
      tauList[[s]] <- tmptaulist
    }
  }

  if(!inherits(tauList[[1]], 'list')){
    stop("MMGFM: each component of tauList must be a list consisting of each matrix for each modality type!")
  }
  mapID <- c("gaussian"=1, "poisson"=2, 'binomial'=3)
  if(is.null(row.names(numvarmat))) stop("MMGFM: the row.names of numvarmat must be some of 'gaussian', 'poisson' and 'bionmial'! \n")
  if(!all(row.names(numvarmat) %in% names(mapID))) stop("MMGFM: the row.names of numvarmat must be some of 'gaussian', 'poisson' and 'bionmial'! \n")

  types <- row.names(numvarmat)
  message("Data include ", S, " studies/sources")
  message("variables belong to ", cc, " types: ", paste(types, collapse = ', '))

  set.seed(seed)
  betalist_int <- list()
  Alist_int <- list()
  Blist_int <- list()
  invLambdalist_int <- list()
  Mulist_y_int <- list()
  Slist_y_int <- list()
  Mlist_int <- list()
  Olist_int <- list()
  wlist_int <- list()
  sigma2_int <- array(dim=c(dd, cc, S))
  Sigma_int <- array(dim=c(q, q, S))
  Phi_int <- list()
  zeta_int <- array(dim=c(dd, cc, S))
  typeID <- unname(mapID[types])
  message("Initialization...")
  ## generate the common intials
  for(s in 1:S){
    # s <- 1
    qs <-  qsvec[s]
    ns <- nrow(XList[[s]][[1]])
    Sigma_int[,, s] <- diag(rep(1, q))
    Phi_int[[s]] <- diag(rep(1, qs)) ### Phi is require using field to put the different shape matrix
    invLamlist_tmp <- list()
    Sylist_tmp <- list()
    wlist_tmp <- list()
    Xlist_tmp <- list()
    for(ic in 1:cc){

      pt <- sum(numvarmat[ic, ]) ## the number of variables in type t.
      invLamlist_tmp[[ic]] <- rep(1, pt)
      if(types[ic] == 'gaussian'){
        Sylist_tmp[[ic]] <- matrix(0, ns, pt)
      }else{
        Sylist_tmp[[ic]] <- matrix(1, ns, pt)
      }

      wlist_tmp[[ic]] <-  matrix(0, ns, dd)
      if(typeID[ic] == 2){
        Xlist_tmp[[ic]] <- log(1+XList[[s]][[ic]])
      }else if(typeID[ic]==1){
        Xlist_tmp[[ic]] <- XList[[s]][[ic]]
      }else if(typeID[ic] ==3){
        Xlist_tmp[[ic]] <- XList[[s]][[ic]] # log(1+XList[[s]][[ic]])
      }

      for(id in 1:dd){
        if(numvarmat[ic,id]>0){
          sigma2_int[id, ic, s] <- 1
          zeta_int[id, ic, s] <- 1
        }
      }
    }

    Mulist_y_int[[s]] <- Xlist_tmp
    Slist_y_int[[s]] <- Sylist_tmp
    wlist_int[[s]] <- wlist_tmp
    invLambdalist_int[[s]] <- invLamlist_tmp
  }

  ## generate the initial method-specified initial values of parameters.


  if(init=='MSFRVI'){
    message("Initialization using MSFR...")
    res_msfr <- initialize_MSFR(XList, ZList, numvarmat,tauList=NULL, q, qs=qsvec, maxIter=15, verbose=FALSE)
    Alist_int <- mat2list(res_msfr$hA, num_modals)
    Blist_int <- lapply(res_msfr$hB, function(x) mat2list(x, num_modals))
    # for(s in 1:S){
    #   nc <- length(Blist_int[[s]])
    #   for(ic in 1:nc){
    #     if(!inherits(Blist_int[[s]][[ic]], 'matrix')){
    #       Blist_int[[s]][[ic]] <- matrix(Blist_int[[s]][[ic]], ncol=1)
    #     }
    #   }
    # }
    Mlist_int <- res_msfr$hF
    Olist_int <- res_msfr$hH
    betalist_int <- mat2list(t(res_msfr$hbeta), num_modals, by_row = FALSE)
    if(d==1){ ## special operation for d=1
      betalist_int <- lapply(betalist_int, function(x) matrix(x, nrow=d))
    }

    message("Finish Initialization!")
  }else{
    for(s in 1:S){
      # s <- 1
      qs <-  qsvec[s]
      ns <- nrow(XList[[s]][[1]])
      Blist_tmp <- list()
      for(ic in 1:cc){

        pt <- sum(numvarmat[ic, ]) ## the number of variables in type t.
        if(typeID[ic] == 2){
          Xtmp <- log(1+XList[[s]][[ic]])
        }else if(typeID[ic]==1){
          Xtmp <- XList[[s]][[ic]]
        }else if(typeID[ic] ==3){
          Xtmp <- XList[[s]][[ic]] # log(1+XList[[s]][[ic]])
        }
        if(init=='random'){
          if(s == 1){
            beta_tmp <- matrix(0, nrow=d,ncol=pt)
            betalist_int[[ic]] <- beta_tmp
            Alist_int[[ic]] <- matrix(runif(pt* q), pt, q)
          }
          Blist_tmp[[ic]] <- matrix(runif(pt* qs), pt, qs)
        }else if(init == 'LFM'){

          lm1 <- lm(Xtmp~ZList[[s]]+0)
          res_lfm <- get_initials(residuals(lm1), q=q+qs)
          if(s==s_max){ ## use the data with larger sample size to obtain initial values of common parameters.
            beta_tmp <- coef(lm1) ## Or in the final, use all study to obtain the initial values for better performance.
            betalist_int[[ic]] <- beta_tmp
            Alist_int[[ic]] <- res_lfm$hB[,1:q]
          }

          Blist_tmp[[ic]] <- res_lfm$hB[,(q+1):(q+qs)]
          if(is.element('gaussian', types)){
            if(types[ic] == 'gaussian'){
              Mlist_int[[s]] <- res_lfm$hH[,1:q]
              Olist_int[[s]] <- res_lfm$hH[,(q+1):(q+qs)]
            }
          }else if(is.element('poisson', types)){
            if(types[ic] == 'poisson'){
              Mlist_int[[s]] <- res_lfm$hH[,1:q]
              Olist_int[[s]] <- res_lfm$hH[,(q+1):(q+qs)]
            }
          }else if(is.element('binomial', types)){
            if(types[ic] == 'binomial'){
              Mlist_int[[s]] <- res_lfm$hH[,1:q]
              Olist_int[[s]] <- res_lfm$hH[,(q+1):(q+qs)]
            }
          }

        }


      }
      Blist_int[[s]] <- Blist_tmp
      if(init=='random'){
        Mlist_int[[s]] <- matrix(rnorm(ns * q), ns, q)
        Olist_int[[s]] <- matrix(rnorm(ns * qs), ns, qs)
      }
    }
  }
  for(s in 1:S){
    nc <- length(Blist_int[[s]])
    for(ic in 1:nc){
      if(!inherits(Blist_int[[s]][[ic]], 'matrix')){
        Blist_int[[s]][[ic]] <- matrix(Blist_int[[s]][[ic]], ncol=1)
      }
    }
    if(!inherits(Olist_int[[s]], 'matrix')){
      Olist_int[[s]] <- matrix(Olist_int[[s]], ncol=1)
    }
  }



  tic <- proc.time()
  reslist <- vb_mmgfmcpp(XList, typeID=typeID, numvarmat, tauList, ZList, betalist_int,
                         Alist_int, Blist_int, invLambdalist_int, sigma2_int, Mulist_y_int,
                         Slist_y_int, Mlist_int, Sigma_int, Olist_int, Phi_int, zeta_int,
                         wlist_int, epsELBO, maxIter, verbose, A_fast=TRUE, add_IC_iter = FALSE)
  toc <- proc.time()
  time.use <- toc[3] - tic[3]
  names_modal <- paste0(rep(types,each=dd), rep(1: dd, times=length(types)))
  names_study <- paste0("study", 1:S)
  vars_modal <- c("hbeta", "hA")
  vars_Smodal <- c("hB", "hv", "hinvLambda")
  vars_study <- c(vars_Smodal, "hF", "hH")

  for(vars in vars_study){

    names(reslist[[vars]]) <- names_study
  }
  attr(reslist[['hSigma']], 'dim3names') <- names_study

  for(s in 1:S){
    if(s==1){
      for(vars in vars_modal){
        # vars <- vars_modal[1]
        names(reslist[[vars]]) <- names_modal
      }
    }
    for(vars in vars_Smodal){
      names(reslist[[vars]][[s]]) <- names_modal
    }
  }

  ## remove the redundent modality
  id_remove <- which(t(numvarmat)==0)
  for(s in 1:S){
    if(s==1){
      for(vars in vars_modal){

        reslist[[vars]][id_remove] <- NULL
      }
    }
    for(vars in vars_Smodal){
      reslist[[vars]][[s]][id_remove]  <- NULL
    }
  }
  tmpList <- add_ident_beta(ZList, reslist$hbeta, reslist$hv, FList=reslist$hF, AList=reslist$hA, method='remove_fv')

  reslist$hbeta <- tmpList$hbeta;reslist$hF <- tmpList$hF;reslist$hv <- tmpList$hV

  reslist$time.use <- time.use
  return(reslist)

}



# Determine the number of factors -----------------------------------------

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

#' Select the number of study-shared and study-specified factors for MMGFM
#' @description Select the number of study-shared and study-specified factors for the high-dimensional multi-study multi-modality covariate-augmented generalized factor model.
#' @param XList a S-length list with each component a m-length list composed by a combined modality matrix of the same type modalities, which is the observed  matrix from each source/study and each modality, where m is the number of modality types.
#' @param ZList a S-length list with each component a matrix that is the covariate matrix from each study.
#' @param numvarmat a m-by-T matrix with rownames modality types that specifies the variable number for each modality of each modality type, where m is the number of modality types, T is the maximum number of modalities for one of modality types  .
#' @param q.max an optional integer, specify the upper bound for the number of study-shared factors; default as 15.
#' @param qsvec.max an optional integer vector with length S, specify the upper bound for the number of study-specifed factors; default as 4 for each study.
#' @param threshold.vec an optional real vector with length 2, specify the threshold for the singular values of study-shared loading and study-specified loading matrices, respectively.
#' @param tauList an optional S-length list with each component a m-length list correponding the offset term for each combined modality of each study; default as full-zero matrix.
#' @param init an optional string, specify the initialization method, supporting "MSFRVI", "random" and "LFM", default as "MSFRVI".
#' @param epsELBO  an optional positive vlaue, tolerance of relative variation rate of the envidence lower bound value, defualt as '1e-5'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed an optional integer, specify the random seed for reproducibility in initialization.
#' @return return a list with two components: q and qs.vec.
#'
#' @export
#'
#' @examples
#' q <- 3; qsvec<-rep(2,3)
#' nvec <- c(100, 120, 100)
#' pveclist <-  list('gaussian'=rep(150, 1),'poisson'=rep(50, 2),'binomial'=rep(60, 2))
#' datlist <- gendata_mmgfm(seed = 1,  nvec = nvec, pveclist =pveclist,
#'                          q = q,  d= 3,qs = qsvec,  rho = rep(3,length(pveclist)), rho_z=0.5,
#'                          sigmavec=rep(0.5, length(pveclist)),  sigma_eps=1)
#' XList <- datlist$XList
#' ZList <- datlist$ZList
#' numvarmat <- datlist$numvarmat
#' ### For illustration, we set maxIter=3. Set maxIter=50 when running formally
#' selectFac.MMGFM(XList, ZList=ZList, numvarmat, q.max=6, qsvec.max  = rep(4,3),
#' init='MSFRVI',maxIter = 3)

selectFac.MMGFM <- function(XList, ZList, numvarmat, q.max=15, qsvec.max=rep(4,length(XList)),
                            threshold.vec = c(1e-2, 1e-3), tauList = NULL,init=c("MSFRVI", "random", "LFM"), epsELBO=1e-12,
                            maxIter=30,  verbose=TRUE, seed=1){

  reslist1 <- MMGFM(XList, ZList=ZList, numvarmat, q=q.max, qsvec = qsvec.max,
                    tauList=tauList,init=init, epsELBO=epsELBO,
                    maxIter=maxIter,  verbose=verbose, seed=seed)
  hqvec <- sapply(reslist1$hA, getratio, threshold=threshold.vec[1])
  hq <- chooseone(hqvec)
  S <- length(reslist1$hB)
  if(hq <2) hq <- 2
  reslist1 <- MMGFM(XList, ZList=ZList, numvarmat, q=hq, qsvec = qsvec.max,
                    tauList=tauList, init=init, epsELBO=epsELBO,
                    maxIter=maxIter,  verbose=verbose, seed=seed)
  hqsvec <- rep(NA, S)
  for(s in 1:S){
    qs <- sapply(reslist1$hB[[s]], getratio, threshold=threshold.vec[2]) ## Use the vote rule.
    hqsvec[s] <- chooseone(qs)
  }
  return(list(q=hq, qs.vec=hqsvec))
}



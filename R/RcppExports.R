# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @keywords internal
NULL

vb_mmgfmcpp <- function(XList, typeID, numvarmat, tauList, Zlist, betalist_int, Alist_int, Blist_int, invLambdalist_int, sigma2_int, Mulist_y_int, Slist_y_int, Mlist_int, Sigma_int, Olist_int, Philist_int, zeta_int, wlist_int, epsELBO, maxIter, verbose, A_fast, add_IC_iter = TRUE) {
    .Call(`_MMGFM_vb_mmgfmcpp`, XList, typeID, numvarmat, tauList, Zlist, betalist_int, Alist_int, Blist_int, invLambdalist_int, sigma2_int, Mulist_y_int, Slist_y_int, Mlist_int, Sigma_int, Olist_int, Philist_int, zeta_int, wlist_int, epsELBO, maxIter, verbose, A_fast, add_IC_iter)
}


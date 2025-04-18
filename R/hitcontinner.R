#' Continuous Hitcalls Inner
#'
#' Calculates continuous hitcall using 3 statistical metrics.
#'
#' This function is called either directly from concRespCore or
#' via hitcont. Details of how to compute function input are in
#' concRespCore.
#'
#' @param conc Vector of concentrations.
#' @param resp Vector of responses.
#' @param top Model predicted top, maximal predicted change in response from baseline.
#' @param cutoff Desired cutoff.
#' @param er Model error parameter.
#' @param ps Vector of used model parameters in order: a, tp, b, ga, p, la, q, er.
#' @param fit_method Name of winning fit method (should never be constant).
#' @param caikwt Akaike weight of constant model relative to winning model.
#' @param mll Maximum log-likelihood of winning model.
#' @param errfun Which error distribution to assume for each point, defaults to
#'   "dt4". "dt4" is the original 4 degrees of freedom t-distribution. Another
#'   supported distribution is "dnorm", the normal distribution.
#' @param poly2.biphasic Which fitting method to use for poly2. If poly2.biphasic = TRUE, allows for biphasic polynomial 2
#'   model fits (i.e. both monotonic and non-monotonic). (Defaults to TRUE.)
#' @param verbose If verbose = TRUE, will print status of empirical calculations. (Defaults to FALSE.)
#'
#' @importFrom stats pt
#' @importFrom stats aggregate
#' @importFrom stats pnorm
#'
#' @return Continuous hitcall between 0 and 1.
#' @export
#'
#' @examples
#' conc = c(.03,.1,.3,1,3,10,30,100)
#' resp = c(0,.1,0,.2,.6,.9,1.1,1)
#' top = 1.023239
#' er = -3.295307
#' ps = c(1.033239, 2.453014, 1.592714, er = -3.295307) #tp,ga,p,er
#' fit_method = "hill"
#' caikwt = 1.446966e-08
#' mll = 12.71495
#' hitcontinner(conc,resp,top,cutoff = 0.8, er,ps,fit_method, caikwt, mll)
#' hitcontinner(conc,resp,top,cutoff = 1, er,ps,fit_method, caikwt, mll)
#' hitcontinner(conc,resp,top,cutoff = 1.2, er,ps,fit_method, caikwt, mll)
hitcontinner = function(conc, resp, top, cutoff, er, ps, fit_method, caikwt, mll, errfun = "dt4", poly2.biphasic = TRUE, verbose = FALSE){

  #Each P represents the odds of the curve being a hit according to different criteria; multiply all Ps to get hit odds overall
  if(fit_method == "none") return(0)
  if(fit_method == "hill") fname = "hillfn" else fname = fit_method

  #caikwt is constant model aikaike weight vs all other models. Represents probability that constant model is correct
  P1 = 1-caikwt

  P2 = 1
  med_resp <- aggregate(resp ~ conc, data = data.frame(conc,resp), FUN = median)
  if (errfun == "dt4") {
    for(y in med_resp$resp){
      #multiply odds of each point falling below cutoff to get odds of all falling below
      P2 = P2*pt((y-sign(top)*cutoff)/exp(er),4, lower.tail = top < 0) #use lower tail for positive top and upper tail for neg top
    }
  } else if (errfun == "dnorm"){
    for(y in med_resp$resp){
      P2 = P2*pnorm((y-sign(top)*cutoff)/exp(er), lower.tail = top < 0)
    }
  }
  P2 = 1- P2 #odds of at least one point above cutoff

  # P3 = pnorm((top-cutoff)/topsd) #odds of top above cutoff
  #assume ps may have nas in them
  ps = ps[!is.na(ps)]
  P3 = toplikelihood(fname, cutoff, conc, resp, ps, top, mll, errfun = errfun, poly2.biphasic = poly2.biphasic, verbose = verbose) #odds of top above cutoff

  #multiply three probabilities
  return(P1*P2*P3)

}

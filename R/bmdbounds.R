#' BMD Bounds
#'
#' Uses maximum likelihood method to tune the upper and lower bounds on the BMD (BMDU, BMDL)
#'
#' Takes in concentration response fit details and outputs a bmdu or bmdl, as
#' desired. If bmd is not finite, returns NA. If the objective function doesn't
#' change sign or the root finding otherwise fails, it returns NA. These
#' failures are not uncommon since some curves just don't reach the desired
#' confidence level.
#'
#' @param fit_method Fit method: "exp2", "exp3", "exp4", "exp5", "hill", "gnls",
#'   "poly1", "poly2", or "pow".
#' @param bmr Benchmark response.
#' @param pars Named vector of model parameters: a,b,tp,ga,p,la,q,er output by
#'   httrfit, and in that order.
#' @param conc Vector of concentrations (NOT in log units).
#' @param resp Vector of responses corresponding to given concentrations.
#' @param onesidedp The one-sided p-value. Default of .05 corresponds to 5
#'   percentile BMDL, 95 percentile BMDU, and 90 percent CI.
#' @param bmd Can optionally input the bmd when already known to avoid
#'   unnecessary calculation.
#' @param which.bound Returns BMDU if which.bound = "upper"; returns BMDL if
#'   which.bound = "lower".
#' @param poly2.biphasic If poly2.biphasic = TRUE, constraints are set to allow
#'   for the polynomial 2 model fit to be bi-phasic (i.e. non-monotonic).
#' @param x_v The vertex of the quadratic/parabolic fit.
#'   Only in use when estimating the BMDL and BMDU values for the "poly2" model
#'   when poly2.biphasic = TRUE. No default is set.
#'
#' @importFrom stats uniroot
#' @return Returns either the BMDU or BMDL.
#' @export
#'
#' @examples
#' conc = c(.03, .1, .3, 1, 3, 10, 30, 100)
#' resp = c(.1,-.1,0,1.1,1.9,2,2.1,1.9)
#' pars = c(tp = 1.973356, ga = 0.9401224, p = 3.589397, er =  -2.698579)
#' bmdbounds(fit_method = "hill", bmr = .5, pars, conc, resp)
#' bmdbounds(fit_method = "hill", bmr = .5, pars, conc, resp, which.bound = "upper")

bmdbounds = function(fit_method, bmr, pars, conc, resp, onesidedp = .05, bmd = NULL, which.bound = "lower", poly2.biphasic = TRUE,x_v){

  # check if the if the biphasic poly2 has the vertex passed
  if(poly2.biphasic & fit_method == "poly2" & missing(x_v)){
    stop("Biphasic Polynomial 2 models require the vertex to estimate BMDL & BMDU.")
  }

  #calculate bmd, if necessary
  if(is.null(bmd)) bmd = acy(bmr, as.list(pars), type = fit_method)
  if(!is.finite(bmd)) return(NA_real_)

  # hill model's function name is hillfn, other models are not changed.
  if(fit_method == "hill") fname = paste0(fit_method, "fn") else fname = fit_method
  maxloglik = tcplObj(p = pars, conc = conc, resp = resp, fname = fname)

  #search for bounds to ensure sign change
  if(which.bound == "lower") {
    xs = 10^seq(-5,log10(bmd), length.out = 100)
    ys = sapply(X = xs, FUN = bmdobj, fname = fname, bmr = bmr, conc = conc, resp = resp, ps = pars, mll = maxloglik,
                onesp = onesidedp, partype = 2,poly2.biphasic = poly2.biphasic,x_v = x_v)
    if(!any(ys >= 0, na.rm = TRUE) | !any(ys < 0, na.rm = TRUE)) return(NA_real_)
    bmdrange = c(max(xs[ys >= 0]), bmd)
  }
  if(which.bound == "upper") {
    if(fit_method == "gnls"){
      toploc = acy(bmr, as.list(pars), type = "gnls", returntoploc = TRUE)
      xs = 10^seq(log10(bmd), log10(toploc), length.out = 100)
    } else xs = 10^seq(log10(bmd), 5, length.out = 100)
    ys = sapply(X = xs, FUN = bmdobj, fname = fname, bmr = bmr, conc = conc, resp = resp, ps = pars, mll = maxloglik,
                onesp = onesidedp, partype = 2,poly2.biphasic = poly2.biphasic,x_v = x_v)
    if(!any(ys >= 0, na.rm = TRUE) | !any(ys < 0, na.rm = TRUE)) return(NA_real_)
    bmdrange = c(bmd, min(xs[ys >= 0]))
  }

  #use type 2 param. only
  out = try(uniroot(bmdobj, bmdrange, fname = fname, bmr = bmr, conc = conc, resp = resp, ps = pars, mll = maxloglik,
                onesp = onesidedp, partype = 2,poly2.biphasic = poly2.biphasic,x_v = x_v)$root)

  #sometimes there's no lower/upper bound because the model is such a poor fit, in this case, return NA
  if(is(out,"try-error")) return(NA_real_) else return(out)

}


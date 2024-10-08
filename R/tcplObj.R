#' Concentration Response Objective Function
#'
#' Log-likelihood to be maximized during CR fitting.
#'
#' This function is a generalized version of the log-likelihood estimation
#' functions used in the ToxCast Pipeline (TCPL).
#' Hill model uses fname "loghill" and gnls uses fname "loggnls". Other model
#' functions have the same fname as their model name; i.e. exp2 uses "exp2", etc.
#' errfun = "dnorm" may be better suited to gsva pathway scores than "dt4".
#' Setting err could be used to fix error based on the null data noise
#' distribution instead of fitting the error when maximizing log-likelihood.
#'
#' @param p Vector of parameters, must be in order: a, tp, b, ga, p, la, q, er.
#'   Does not require names.
#' @param conc Vector of concentrations in log10 units for loghill/loggnls, in
#'   regular units otherwise.
#' @param resp Vector of corresponding responses.
#' @param fname Name of model function.
#' @param errfun Which error distribution to assume for each point, defaults to
#'   "dt4". "dt4" is the original 4 degrees of freedom t-distribution. Another
#'   supported distribution is "dnorm", the normal distribution.
#' @param err An optional estimation of error for the given fit.
#'
#' @importFrom stats dt dnorm
#'
#' @return Log-likelihood.
#' @export
#'
#' @examples
#' conc = c(.03,.1 , .3 , 1  , 3 , 10 , 30  , 100)
#' resp = c( 0 , 0 , .1 ,.2 , .5 , 1  , 1.5 , 2  )
#' p = c(tp = 2, ga = 3, p = 4, er = .5)
#' tcplObj(p,conc,resp,"exp5")
#'
#' lconc = log10(conc)
#' tcplObj(p,lconc,resp,"loghill")
tcplObj = function(p, conc, resp, fname, errfun = "dt4", err = NULL) {

  mu = do.call(fname, list(ps = p, x = conc)) #get model values for each conc
  n = length(p)
  if(is.null(err)) err = exp(p[n]) #set error term

  #objective function is sum of log-likelihood of response given the model at each concentration
  #scaled by variance (err)
  if(errfun == "dt4") return( sum( dt((resp - mu)/err, df = 4, log = TRUE) - log(err) ) )
  if(errfun == "dnorm") return( sum( dnorm((resp - mu)/err, log = TRUE) - log(err) ) )

}

#' Constant Model
#'
#' \eqn{f(x) = 0}
#'
#' @param ps Vector of parameters (ignored)
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' cnst(1,1)
#'
cnst = function(ps,x){
  #ignores ps
  return(rep(0,length(x)))
}

#' Exponential 2 Model
#'
#' \eqn{f(x) = a*(e^{(x/b)}- 1)}
#'
#' @param ps Vector of parameters: a,b,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' exp2(c(1,2),1)
#'
exp2 = function(ps,x){
  #a = ps[1], b = ps[2]
  return(ps[1]*(exp(x/ps[2]) - 1)  )
}

#' Exponential 3 Model
#'
#' \eqn{f(x) = a*(e^{(x/b)^p} - 1)}
#'
#' @param ps Vector of parameters: a,b,p,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' exp3(c(1,2,2),1)
#'
exp3 = function(ps,x){
  #a = ps[1], b = ps[2], p = ps[3]
  return(ps[1]*(exp((x/ps[2])^ps[3]) - 1)  )
}

#' Exponential 4 Model
#'
#' \eqn{f(x) = tp*(1-2^{(-x/ga)})}
#'
#' @param ps Vector of parameters: tp,ga,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' exp4(c(1,2),1)
#'
exp4 = function(ps,x){
  #tp = ps[1], ga = ps[2]
  return(ps[1]*(1-2^(-x/ps[2]))  )
}

#' Exponential 5 Model
#'
#' \eqn{f(x) = tp*(1-2^{(-(x/ga)^p)})}
#'
#' @param ps Vector of parameters: tp,ga,p,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' exp5(c(1,2,3),1)
#'
exp5 = function(ps,x){
  #tp = ps[1], ga = ps[2], p = ps[3]
  return(ps[1]*(1-2^(-(x/ps[2])^ps[3]))  )
}

#' Gain-Loss Model
#'
#' \eqn{f(x) = \frac{tp}{[(1 + (ga/x)^p )(1 + (x/la)^q )]}}
#'
#' @param ps Vector of parameters: tp,ga,p,la,q,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' gnls(c(1,2,1,2,2),1)
#'
gnls = function(ps, x){
  #gnls function with regular units
  #tp = ps[1], ga = ps[2], p = ps[3], la = ps[4], q = ps[5]
  gn <- 1/(1 + (ps[2]/x)^ps[3])
  ls <- 1/(1 + (x/ps[4])^ps[5])
  return(ps[1]*gn*ls )
}

#' Log Gain-Loss Model
#'
#' \eqn{f(x) = \frac{tp}{[(1 + 10^{(p*(ga-x))} )(1 + 10^{(q*(x-la))} )]}}
#'
#' @param ps Vector of parameters: tp,ga,p,la,q,er (ga and la are in log10-scale)
#' @param x Vector of concentrations (log10 units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' loggnls(c(1,2,1,2,2),1)
#'
loggnls = function(ps, x){
  #gnls function with log units: x = log10(conc) and ga/la = log10(gain/loss ac50)
  #tp = ps[1], ga = ps[2], p = ps[3], la = ps[4], q = ps[5]
  gn <- 1/(1 + 10^((ps[2] - x)*ps[3]))
  ls <- 1/(1 + 10^((x - ps[4])*ps[5]))
  return(ps[1]*gn*ls )
}

#' Hill Model
#'
#' \eqn{f(x) = \frac{tp}{[(1 + (ga/x)^p )]}}
#'
#' @param ps Vector of parameters: tp,ga,p,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' hillfn(c(1,2,3),1)
#'
hillfn = function(ps,x){
  #hill function with regular units
  #tp = ps[1], ga = ps[2], p = ps[3]
  return(ps[1]/(1 +  (ps[2]/x)^ps[3]) )
}

#' Log Hill Model
#'
#' \eqn{f(x) = \frac{tp}{(1 + 10^{(p*(ga-x))} )}}
#'
#' @param ps Vector of parameters: tp,ga,p,er (ga is in log10-scale)
#' @param x Vector of concentrations (log10 units)
#'
#' @return Vector of model responses.
#' @export
#' @examples
#' loghill(c(1,2,3),1)
#'
loghill = function(ps,x){
  #hill function with log units: x = log10(conc) and ga = log10(ac50)
  #tp = ps[1], ga = ps[2], p = ps[3]
  return(ps[1]/(1 + 10^(ps[3]*(ps[2]-x)) ) )
}

#' Polynomial 1 Model
#'
#' \eqn{f(x) = a*x}
#'
#' @param ps Vector of parameters: a,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' poly1(1,1)
#'
poly1 = function(ps,x){
  #a = ps[1]
  return(ps[1]*x)
}

#' Polynomial 2 Model
#'
#' \eqn{f(x) = a*(\frac{x}{b} + \frac{x^2}{b^2})}
#'
#' @param ps Vector of parameters: a,b,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' poly2(c(1,2),1)
#'
poly2 = function(ps,x){
  #a = ps[1], b = ps[2]
  x0 = x/ps[2]
  return(ps[1]*(x0 + x0*x0))
}

#' Polynomial 2 Model (BMDS)
#'
#' \eqn{f(x) = b1*x + b2*x^{2}}
#'
#' @param ps Vector of parameters: b1,b2,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' poly2bmds(c(1,2),1)
#'
poly2bmds = function(ps,x){
  #b1 = ps[1], b2 = ps[2]
  return(ps[1]*x + ps[2]*(x)^2)
}

#' Power Model
#'
#' \eqn{f(x) = a*x^p}
#'
#' @param ps Vector of parameters: a,p,er
#' @param x Vector of concentrations (regular units)
#'
#' @return Vector of model responses
#' @export
#' @examples
#' pow(c(1,2),1)
#'
pow = function(ps,x){
  #a = ps[1], p = ps[2]
  return(ps[1]*x^ps[2]  )
}

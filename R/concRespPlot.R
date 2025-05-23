#' Concentration Response Plot
#'
#' Plots a concentration response curve for one sample/endpoint combination.
#' This is a generic function and it is expected that users will make their own versions
#'
#' row is one row of data from concRespCore
#'
#' @param row Named list containing:
#'   \itemize{
#'     \item conc - conc string separated by |'s
#'     \item resp - response string separated by |'s
#'     \item method - scoring method determines plot bounds
#'     \item name - chemical name for plot title
#'     \item cutoff - noise cutoff
#'     \item bmr - baseline median response; level at which bmd is calculated
#'     \item er - fitted error term for plotting error bars
#'     \item a, tp, b, ga, p, la, q - other model parameters for fit curve
#'     \item fit_method - curve fit method
#'     \item bmd, bmdl, bmdu - bmd, bmd lower bound, and bmd upper bound
#'     \item ac50, acc - curve value at 50\% of top, curve value at cutoff
#'     \item top - curve top (maximal predicted change in response from baseline)
#'     \item name - name of the chemical
#'     \item assay - name of the assay, signature, or other endpoint
#'     \item other identifiers
#'   }
#'   Other elements are ignored.
#' @param ymin Minimum value of response for the plot
#' @param ymax Maximum value of response for the plot
#' @param draw.error.arrows If TRUE, draw lines representing the uncertainty in the response estimate,
#'   instead of the actual response points
#' @return No output.
#' @export
#'
#' @import stringr
#' @import grDevices
#' @importFrom stats qt
#'
#' @examples
#'
#' conc <- list(.03, .1, .3, 1, 3, 10, 30, 100)
#' resp <- list(0, .2, .1, .4, .7, .9, .6, 1.2)
#' row <- list(conc = conc,
#'             resp = resp,
#'             bmed = 0,
#'             cutoff = 0.25,
#'             onesd = 0.125,
#'             name = "some chemical",
#'             assay = "some assay")
#' res <- concRespCore(row, conthits = TRUE)
#' concRespPlot(res,ymin=-2.5,ymax=2.,5)
#'
#'
concRespPlot <- function(row,ymin=-120,ymax=120,draw.error.arrows=FALSE) {

  #every variable in row goes into the environment to make it easy
  #to update this function to use new row data.
  list2env(row,envir = environment())
  name <- unlist(name)
  assay <- unlist(assay)
  cutoff <- unlist(cutoff)
  bmr <- unlist(bmr)
  er <- unlist(er)
  fit_method <- unlist(fit_method)
  ac50 <- unlist(ac50)
  top <- unlist(top)
  bmd <- unlist(bmd)
  acc <- unlist(acc)
  hitcall <- unlist(hitcall)
  bmdl <- unlist(bmdl)
  bmdu <- unlist(bmdu)

  #reformat conc and resp as vectors
  conc <- as.numeric(str_split(row[1,"conc"],"\\|")[[1]])
  resp <- as.numeric(str_split(row[1,"resp"],"\\|")[[1]])

  # replace untreated controls with pseudo value
  if (any(conc==0)) warning("Data contains untreated controls (conc = 0). A pseudo value replaces -Inf after log-transform.  The pseudo value is set to one log-unit below the lowest experimental `conc`.")
  logconc <- log10(conc)
  # replace the negative infinity with a number that is one log-10 unit
  # less than the second lowest dose (in log).
  logconc <- replace(logconc, logconc == -Inf, sort(unique(logconc))[2]-1)
  conc <- 10**logconc

  #plotting points for curves based on min/max experimental conc
  conc_plot <- 10**seq(from = log10(min(conc)), to = log10(max(conc)), by = 0.05)

  #some deprecated code; later will use j =1 and col.list[j] to mean black
  col.list <- c("black","cyan","red")

  #empty plot to start with
  plotrange = c(min(conc),max(conc))
  plot(c(1,1),type="n",xlab="conc (uM)",ylab="Response",xlim=plotrange,ylim=c(ymin,ymax),
       log="x",main=paste(name,"\n",assay),cex.main=0.9)

  #cutoffs and gray rectangular noise region
  rect(xleft=plotrange[1],ybottom=-cutoff,xright=plotrange[2],ytop=cutoff,col="lightgray")
  lines(plotrange,c(cutoff,cutoff),lwd=1)
  lines(plotrange,c(-cutoff,-cutoff),lwd=1)

  #thick line at 0 and bmrs
  lines(plotrange,c(0,0),lwd=2)
  lines(plotrange,c(bmr,bmr),lwd=1)
  lines(plotrange,c(-bmr,-bmr),lwd=1)

  #height for top labels
  yplot <- ymax*0.95
  #equally spaced labelling points pegged to plotrange
  xplot = 10^(seq(log10(plotrange[1]), log10(plotrange[2]), length.out = 8))[-8]

  #Top label headings
  text(xplot[1],yplot,"mthd",pos=4)
  text(xplot[2],yplot,"AC50",pos=4)
  text(xplot[3],yplot,"Top",pos=4)
  text(xplot[4],yplot,"BMD",pos=4)
  text(xplot[5],yplot,"ACC",pos=4)
  text(xplot[6],yplot,"Hitcall",pos=4)

  j = 1 #j = 1 is black
  conc <- conc[!is.na(resp)]
  resp <- resp[!is.na(resp)]

  #draw error bars according to er, using 95% confidence for t-dist w/ four df
  if(!is.na(er) && draw.error.arrows){
    arrows(conc, resp+exp(er)*qt(.025,4), conc, resp+exp(er)*qt(.975,4), length = .05, angle = 90, code = 3)
  }
  else points(conc, resp, pch=21,bg="black")

  #get model parameters
  parnames = c("a", "tp", "b", "ga", "p", "la", "q")
  modpars = as.list(row[,parnames])
  modpars= modpars[!sapply(modpars, is.na)]

  #gcalculate and plot model curves
  if(fit_method == "hill"){
    resp_plot <- do.call("hillfn",list(ps = unlist(modpars), x = conc_plot))
    lines(resp_plot~conc_plot,col=col.list[j])
  } else if(!fit_method %in% c("cnst","none") ){
    resp_plot <- do.call(fit_method,list(ps = unlist(modpars), x = conc_plot))
    lines(resp_plot~conc_plot,col=col.list[j])
  }

  yplot <- yplot-(ymax-ymin)*0.05 #second row for top labels

  #Fill in top labels second row
  text(xplot[1],yplot,fit_method,pos=4)
  text(xplot[2],yplot,format(ac50,digits=2),pos=4, col = "red")
  text(xplot[3],yplot,format(top,digits=2),pos=4)
  text(xplot[4],yplot,format(bmd,digits=2),pos=4, col = "green")
  text(xplot[5],yplot,format(acc,digits=2),pos=4, col = "blue")

  #color hitcall based on whether it's a hit
  color <- "black"
  font <- 1
  if(hitcall==1) {
    color <- "red"
    font <- 2
  }
  text(xplot[6],yplot,format(hitcall, digits = 2),pos=4,col=color,cex=1,font=font)

  #plot green bmd with range
  if(hitcall>0) {
    lines(c(bmd,bmd),c(ymin/2,ymax/2),col="green",lwd=2, lty = isTRUE(bmd<min(conc)) + 1)
    if(is.na(bmdl)) xleft = plotrange[1]/10 else xleft = bmdl
    if(is.na(bmdu)) xright = plotrange[2]*10 else xright = bmdu

    rect(xleft=xleft,ybottom=ymin/2,xright=xright,ytop=ymax/2,col=rgb(0,1,0, alpha = .5), border = NA)
    lines(c(xleft,xleft),c(ymin/2,ymax/2),col="green",lwd=1)
    lines(c(xright,xright),c(ymin/2,ymax/2),col="green",lwd=1)
  }
}

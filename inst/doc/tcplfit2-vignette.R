## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(tcplfit2)

## ----example1, fig.height = 7, fig.width = 7-----------------------------
  conc <- list(.03,.1,.3,1,3,10,30,100)
  resp <- list(0,.2,.1,.4,.7,.9,.6, 1.2)
  row = list(conc = conc, resp = resp, bmed = 0, cutoff = 1, onesd = .5,name="some chemical")
  res <- concRespCore(row,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                        "exp4", "exp5"),conthits = T, do.plot=T)
  

## ----example2, fig.height = 8, fig.width = 7, eval = FALSE---------------
#    # read in the data
#    file <- "data/mc3.RData"
#    load(file=file)
#  
#    # set up a 3 x 2 grid for the plots
#    oldpar <- par(no.readonly = TRUE)
#    on.exit(par(oldpar))
#    par(mfrow=c(3,2),mar=c(4,4,2,2))
#  
#    # determine the background variation
#    temp <- mc3[mc3$logc<= -2,"resp"]
#    bmad <- mad(temp)
#    onesd <- sd(temp)
#    cutoff <- 3*bmad
#  
#    # select six samples. Note that there may be more than one sample processed for a given chemical
#    spid.list <- unique(mc3$spid)
#    spid.list <- spid.list[1:6]
#  
#    for(spid in spid.list) {
#      # select the data for just this sample
#      temp <- mc3[is.element(mc3$spid,spid),]
#  
#      # The data file has stored concentration in log10 form, so fix that
#      conc <- 10**temp$logc
#      resp <- temp$resp
#  
#      # pull out all of the chemical identifiers and the name of the assay
#      dtxsid <- temp[1,"dtxsid"]
#      casrn <- temp[1,"casrn"]
#      name <- temp[1,"name"]
#      assay <- temp[1,"assay"]
#  
#      # create the row object
#      row <- list(conc = conc, resp = resp, bmed = 0, cutoff = cutoff, onesd = onesd,assay=assay,dtxsid=dtxsid,casrn=casrn,name=name)
#  
#      # run the concentration-response modeling for a single sample
#      res <- concRespCore(row,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
#                                            "exp4", "exp5"),conthits = T, aicc = F,bidirectional=F)
#  
#      # plot the results
#      concRespPlot(res,ymin=-10,ymax=100)
#    }
#  
#  

## ----example3, fig.height = 6, fig.width = 7-----------------------------
  # call additional R packages
  library(stringr)  # string management package

  # read in the file
  data("signatures")
  
  # set up a 3 x 2 grid for the plots
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  par(mfrow=c(3,2),mar=c(4,4,2,2))
    
  # fit 6 observations in signatures
  for(i in 1:nrow(signatures)){
    # set up input data
    row = list(conc=as.numeric(str_split(signatures[i,"conc"],"\\|")[[1]]),
               resp=as.numeric(str_split(signatures[i,"resp"],"\\|")[[1]]),
               bmed=0,
               cutoff=signatures[i,"cutoff"],
               onesd=signatures[i,"onesd"],
               name=signatures[i,"name"],
               assay=signatures[i,"signature"])
    # run concentration-response modeling (1st plotting option)
    out = concRespCore(row,conthits=F,do.plot=T)
    if(i==1){
      res <- out
    }else{
      res <- rbind.data.frame(res,out)
    }
  }

## ----example3_plot2, fig.height = 8, fig.width = 7-----------------------
  # set up a 3 x 2 grid for the plots
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  # plot results using `concRespPlot`(2nd plotting option)
  for(i in 1:nrow(res)){
    concRespPlot(res[i,],ymin=-1,ymax=1)
  }


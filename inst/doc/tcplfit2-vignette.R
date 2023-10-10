## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning = FALSE, echo = FALSE-------------------------------------
library(tcplfit2)

## ----example1, fig.height = 4.55, fig.width = 8-------------------------------
  conc <- list(.03,.1,.3,1,3,10,30,100)
  resp <- list(0,.2,.1,.4,.7,.9,.6, 1.2)
  row = list(conc = conc, resp = resp, bmed = 0, cutoff = 1, onesd = .5,name="some chemical")
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))   
  par(xpd = TRUE)
  res <- concRespCore(row,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                        "exp4", "exp5"),conthits = T, do.plot=T)
  

## ----example1 result, warning=FALSE, echo=FALSE-------------------------------
library(DT)
DT::datatable(res,rownames = FALSE,options = list(scrollX = T))

## ----example2, fig.height = 8, fig.width = 7----------------------------------
  # read in the data
  # Loading in the level 3 example data set from invitrodb
  data("mc3")

  # set up a 3 x 2 grid for the plots
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  # determine the background variation
  temp <- mc3[mc3$logc<= -2,"resp"]
  bmad <- mad(temp)
  onesd <- sd(temp)
  cutoff <- 3*bmad

  # select six samples. Note that there may be more than one sample processed for a given chemical
  spid.list <- unique(mc3$spid)
  spid.list <- spid.list[1:6]

  for(spid in spid.list) {
    # select the data for just this sample
    temp <- mc3[is.element(mc3$spid,spid),]

    # The data file has stored concentration in log10 form, so fix that
    conc <- 10**temp$logc
    resp <- temp$resp

    # pull out all of the chemical identifiers and the name of the assay
    dtxsid <- temp[1,"dtxsid"]
    casrn <- temp[1,"casrn"]
    name <- temp[1,"name"]
    assay <- temp[1,"assay"]

    # create the row object
    row <- list(conc = conc, resp = resp, bmed = 0, cutoff = cutoff, onesd = onesd,assay=assay,dtxsid=dtxsid,casrn=casrn,name=name)

    # run the concentration-response modeling for a single sample
    res <- concRespCore(row,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                          "exp4", "exp5"),conthits = T, aicc = F,bidirectional=F)

    # plot the results
    concRespPlot(res,ymin=-10,ymax=100)
  }



## ----example2 result, warning=FALSE, echo=FALSE-------------------------------
DT::datatable(res,rownames = FALSE,options = list(scrollX = T))

## ----example3, fig.height = 6, fig.width = 7, warning = FALSE-----------------
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

## ----example3_plot2, fig.height = 8, fig.width = 7----------------------------
  # set up a 3 x 2 grid for the plots
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  # plot results using `concRespPlot`(2nd plotting option)
  for(i in 1:nrow(res)){
    concRespPlot(res[i,],ymin=-1,ymax=1)
  }

## ----example4_init, fig.height = 6, fig.width = 7, message=FALSE, warning = FALSE,echo=-4----
# Loading in the level 0 example data set from invitrodb
data("mc0")
library(data.table)
data.table::setDTthreads(2)
dat <- mc0
DT::datatable(head(dat[wllt=='t',]),rownames= FALSE, options = list(scrollX = T))

## ----example4_cndx, fig.height = 6, fig.width = 7-----------------------------
library(tcpl)
## Order by the following columns
setkeyv(dat, c('acid', 'srcf', 'apid', 'coli', 'rowi', 'spid', 'conc'))

## Define replicate id (rpid) column for test compound wells
nconc <- dat[wllt == "t" , ## denotes test well as the well type (wllt)
             list(n = lu(conc)), #total number of unique concentrations
             by = list(acid, apid, spid)][ , list(nconc = min(n)), by = acid]
dat[wllt == "t" & acid %in% nconc[nconc > 1, acid],
    rpid := paste(acid, spid, wllt, srcf, apid, "rep1", conc, sep = "_")]
dat[wllt == "t" & acid %in% nconc[nconc == 1, acid],
    rpid := paste(acid, spid, wllt, srcf, "rep1", conc, sep = "_")]

## Define rpid column for non-test compound wells
dat[wllt != "t",
    rpid := paste(acid, spid, wllt, srcf, apid, "rep1", conc, sep = "_")]

## set repid based on rowid
dat[, dat_rpid := rowid(rpid)]
dat[, rpid := sub("_rep[0-9]+.*", "",rpid, useBytes = TRUE)]
dat[, rpid := paste0(rpid,"_rep",dat_rpid)]

# Define concentration index
indexfunc <- function(x) as.integer(rank(unique(x))[match(x, unique(x))])
dat[ , cndx := indexfunc(conc), by = list(rpid)]

## ----example4_mc2, fig.height = 6, fig.width = 7------------------------------
# If no adjustments are required for the data, the corrected value (cval) should be set as original rval
dat[,cval := rval]

## Poor well quality (wllq) wells should be removed
dat <- dat[!wllq == 0,]

## Fitting generally cannot occur if response values are NA therefore values need to be removed
dat <- dat[!is.na(cval),]

## A column for log10 concentration is added as some of the mc3 methods require logc. Given logging concentration, conc=0 are not allowed therefore a dummy non-zero value should be used
dat[conc == 0 , conc := 0.0001]
dat[ , logc := log10(conc)]

#As a final step to prepare the dataset tcplfit2 processing, a dummy aeid is required if using mc3_mthds from tcpl
dummy_aeid <- 99999
dat[,aeid := dummy_aeid]

## Set aeid as a key
setkey(dat,aeid)

## ----example4_mthdlist, fig.height = 6, fig.width = 7, warning = FALSE--------
mthd_funcs <- tcpl:::mc3_mthds()
DT::datatable(tcpl::tcplMthdList(3),rownames= FALSE, options = list(scrollX = T))

## ----example4_mc3methods, fig.height = 6, fig.width = 7, results = 'hide'-----
# apply level 3 methods
## These methods directly apply the normalization methods from tcpl without the need for a DB connection
lapply(mthd_funcs[["bval.apid.nwlls.med"]](dummy_aeid), eval)
lapply(mthd_funcs[["pval.apid.medncbyconc.min"]](dummy_aeid),eval)
lapply(mthd_funcs[["resp.pc"]](dummy_aeid),eval)

## ----example4_mthdlist_l4, fig.height = 6, fig.width = 7----------------------
mthd_funcs_l4 <- tcpl:::mc4_mthds()
DT::datatable(tcpl::tcplMthdList(4), rownames= FALSE, options = list(scrollX = T))

## ----example4_mc4methods, fig.height = 6, fig.width = 7, results = 'hide'-----
# apply level 4 methods
## These methods directly apply the noise calculation and fitting methods from tcpl without the need for a DB connection
lapply(mthd_funcs_l4[["bmad.aeid.lowconc.twells"]](),eval)
lapply(mthd_funcs_l4[["onesd.aeid.lowconc.twells"]](),eval)
lapply(mthd_funcs_l4[["bidirectional.false"]](),eval)

## ----example4_fitting, fig.height = 6, fig.width = 7--------------------------
#do tcplfit2 fitting
myfun <- function(y) {
  res <- tcplfit2::tcplfit2_core(y$conc,
                          y$resp,
                          cutoff = unique(y$bmad),
                          bidirectional = TRUE,
                          verbose = FALSE,
                          force.fit = TRUE,
                          fitmodels = c("cnst", "hill", "gnls", "poly1",
                                        "poly2", "pow", "exp2", "exp3",
                                        "exp4", "exp5")
                          )
  list(list(res)) #use list twice because data.table uses list(.) to look for values to assign to columns
}

## ----example4_fitting_full, eval=FALSE----------------------------------------
#  # only want to run tcplfit2 for test wells in this case
#  # this chunk doesn't run, fit the curves on the subset below
#  dat[wllt == 't',params:= myfun(.SD), by = .(spid)]

## ----example4_fitting_subset--------------------------------------------------
# create a subset that contains 6 samples and run curve fitting
subdat <- dat[spid %in% unique(spid)[10:15],]
subdat[wllt == 't',params:= myfun(.SD), by = .(spid)]

## ----example4_hitcalling, fig.height = 6, fig.width = 7-----------------------
myfun2 <- function(y) {
  res <- tcplfit2::tcplhit2_core(params = y$params[[1]],
                                 conc = y$conc,
                                 resp = y$resp,
                                 cutoff = 3*unique(y$bmad),
                                 onesd = unique(y$osd)
                                 )
  list(list(res))
}

# continute with hitcalling
res <- subdat[wllt == 't', myfun2(.SD), by = .(spid)]

#pivot wider
res_wide <- rbindlist(Map(cbind, spid = res$spid, res$V1))

DT::datatable(res_wide,options = list(scrollX = T))

## ----example4_plot, fig.height = 8, fig.width = 7-----------------------------
  # set up a 3 x 2 grid for the plots
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))            
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  # plot results using `concRespPlot`(2nd plotting option)
  for(i in 1:nrow(res)){
    concRespPlot(res_wide[i,],ymin=-50,ymax=50)
  }


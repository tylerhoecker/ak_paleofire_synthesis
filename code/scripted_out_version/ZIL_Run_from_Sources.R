zil_fn <- function(lakeList){
  
  ## Import raw charcoal count data, standardize and derive CHAR
  source(file.path('code','ZIL_ImportCHAR.r'))
}  
  # Load functions that define ZIL distribution
  source(file.path('code','ZIL_ZILnorm.r'))
  source(file.path('code','ZIL_Kernels.r'))
  
  # Set up variables and parameters for ZIL estimation procedure
  source(file.path('code','ZIL_Setup.r'))
  
  # Weight CHAR so that all regions contribute equally to 
  if(all(lakeList == alaska)){
    source(file.path('code','ZIL_WeightCHAR.r'))}
  
  # Maxium-likelihood estimatation of ZIL distribution within each window at each timestep.
  source(file.path('code','ZIL_Estimate.r'))
  
  # Summarize output
  output = data.frame(
    yr.bp     = xx,
    composite.mean = meanMLE,
    composite.med  = medMLE,
    n.records			 = n.recs.xx,
    p.hat     = fitparms[,1],
    mu.hat    = fitparms[,2],
    sigma.hat = fitparms[,3],
    n.eff     = n.eff,
    CI.lower  = out.ci.L,
    CI.upper  = out.ci.U,
    mean.se   = out.mean.se)

  return(output)
}


  
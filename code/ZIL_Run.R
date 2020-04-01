# This script defines the objects below, and then creates a function which
# implements the procedure presented in Kelly et al. 2013 for estimating a 
# composite index of biomass burning from multiple lake-sediment charcoal 
# records based on the means of zero-inflated log-normal distributions 
# within 100-year and 500-year wide moving windows, estimated at 
# 10-year timesteps, and 90% bootstrapped confidence intervals. 

# Define list of regions
regions  <- list('noatak'= noatak,'brooks'= brooks,'yukon'= yukon, 'copper'= copper) %>% 
  utils::stack() %>% 
  rename(lake = values, region = ind)
# Alaska-wide "region"
alaska <- c(noatak,brooks,yukon,copper)

zil_fn <- function(lakeList){
  
  # Import raw charcoal count data, standardize and derive CHAR
  #------------------------------------------------------------------------
  
  # Function for standardizing non-zero charcoal accumlations rates
  trans_fn <- function(x) {
    x = ifelse(x > 0, x, NA)
    logX = log(x)
    zX = (logX - mean(logX, na.rm = T)) / sd(logX, na.rm = T)
    expX = exp(zX)
    expX[is.na(expX)] = 0
    # TEST
    #expX[expX > 4] = 4
    
    return(expX)
  }
  
  # Calculate and standardize CHAR
  charData <- map(file.path('data','charData',paste0(lakeList,'_charData.csv')),
                  read_csv) %>%
    `names<-` (lakeList) %>%
    bind_rows(.id = 'lake') %>% 
    group_by(lake) %>% 
    mutate(sedAcc = (cmTop - cmBot) / (ageTop - ageBot),
           rawChar = (charCount / charVol) * sedAcc,
           char = trans_fn(rawChar)) %>% 
    rowwise() %>% 
    mutate(age = round(mean(c(ageTop,ageBot)))) %>% 
    full_join(regions, by = 'lake') %>% 
    select(char, age, lake, region) %>% 
    as.data.frame() 
  
  #------------------------------------------------------------------------

  # Load functions that define ZIL distribution
  #------------------------------------------------------------------------
  # Wrapper function returning kernel weights (K), standardized weights (w), and the effective sample size (n). I made up the formula for the last one and should probably investigate it further someday.
  ry.kernel = function(x, x0, h, kernel.fn) {
    dists = abs(x-x0)
    
    K  = kernel.fn(dists,h)
    #		K[which(K<1e-10)] = 0; # Minimum wernel weight to be considered (think machine precision)
    K0 = kernel.fn(0,h)
    
    if(sum(K)==0) { 
      # The total weight of the data is 0 (e.g. x0 is too far away). Assign manually because can't divide by sum(K).
      w = numeric(length(K));
    } else {
      # Otherwise, just divide by sum(K) so weights add to 1
      w = K/sum(K);
    }
    
    n = sum(K/K0)
    # Derived independently but equivalently to:
    #		Chaudhuri, P., and J. S. Marron. 1999. SiZer for exploration of structures in curves. Journal of the American Statistical Association:807–823. PAGE 812
    
    return(list(w=w,n=n,K=K))
  }
  
  #------------------ Individual kernels
  
  # Uniform kernel with half-window width = h
  kernel.unif = function(dists, h) {
    inWin = dists <= h
    
    K = numeric(length(dists))
    K[inWin] = 1
    
    return(K)
  }
  
  
  # Gaussian kernel with SD=h
  kernel.gauss = function(dists, h) {
    K = dnorm(dists,0,h)
    return(K)
  }
  
  #------------------ Parameter order for all funtions below
  # Zero-inflated lognormal
  # parms[1] = p 			= probability of zero
  # parms[2] = mu 		= mean of lognormal dist
  # parms[3] = sigma  = sd of lognormal dist
  
  #------------------ ZIlnorm.randomdraw
  # Produce a random draw from the zero-inflated lognormal distribution
  ZIlnorm.randomdraw = function(parms,n) {
    (runif(n,0,1)>parms[1]) * rlnorm(n, parms[2], parms[3])
  }
  
  
  #------------------ ZIlnorm.mle
  # Analytical solution for MLE parameters of zero-inflated log-normal distribution
  ZIlnorm.mle = function(sampledata) {
    n=length(sampledata)
    mle.lnorm.p    = sum(sampledata==0)/n
    mle.lnorm.mean = mean(log(sampledata[sampledata>0])) 
    mle.lnorm.sd   = 
      sqrt( sum( (log(sampledata[sampledata>0]) - mle.lnorm.mean)^2 )/sum(sampledata>0) ) 
    mleparms 			 = c(mle.lnorm.p, mle.lnorm.mean, mle.lnorm.sd)
    return(mleparms)
  }
  
  
  #------------------ ZIlnorm.llk
  # Log-likelihood for zero-inflated lognormal
  ZIlnorm.llk = function(parms,x,w) {
    # Weight vector w can be omitted, in which case values are given equal weight
    if( missing(w) ) { w = rep(1,length(x)) }
    
    # This code works in several separated steps:
    #	ind0 = which(x==0)
    #	ind1 = which(x>0)
    #
    #	llike0 = sum( w[ind0] * log(parms[1]) )
    #	llike1 = sum( w[ind1] * (log(1-parms[1]) + dlnorm(x[ind1],parms[2],parms[3], log=TRUE)) )
    #
    #	return( -1* (llike0 + llike1) )
    
    # But this seems a bit faster:
    return( -1* (
      sum( w[x==0] * log(parms[1]) ) +
        sum( w[x>0] * (log(1-parms[1]) + dlnorm(x[x>0],parms[2],parms[3], log=TRUE)) )
    ))
  }
  
  
  ZIlnorm.llk2 = function(parms,x,w) {
    # Weight vector w can be omitted, in which case values are given equal weight
    if( missing(w) ) { w = rep(1,length(x)) }
    
    return( -1*sum( w[x>0] * dlnorm(x[x>0],parms[1],parms[2], log=TRUE) ) )
  }
  
  
  #------------------ ZIlnorm.mean
  # Calculate the mean of a specified ZIlnorm distribution. Works on a 3-parameter vector or an array with 3 parameters per row. 
  ZIlnorm.mean = function(parms) {
    if( length(parms) == 3) { # it's just a vector of three parms
      return(  (1-parms[1]) * (exp(parms[2] + (parms[3]^2)/2))  )
    } else { # it's an array with 3 parms per row
      n.rows = dim(parms)[1]
      output = numeric(n.rows)
      for( i in 1:n.rows) {
        output[i] = (1-parms[i,1]) * (exp(parms[i,2] + (parms[i,3]^2)/2))
      }
      return(output)
    }
  }
  
  ZIlnorm.median = function(parms) {
    if(parms[1]>=0.5) { # at least 50% of points == 0
      return(0)
    } else {		
      return(  qlnorm( ((0.5-parms[1])/(1-parms[1])) ,parms[2],parms[3])  )
    }
  }
  #------------------------------------------------------------------------
  
  # Set up variables and parameters for ZIL estimation procedure
  #------------------------------------------------------------------------
  # Set up data for ZIL procedure
  # Temporal params
  xxStep      <- 10				    # Interval, in years, for moving window
  
  # Use cut off study period for regions with shorter records
  if ( identical(lakeList, copper) ) {
    ageLim = c(-60,7000)
      } else if ( identical(lakeList, noatak) ) {
        ageLim = c(-60,7000)
      } else {
          ageLim = c(-60,10000)}
  
  # Confidence Intervals parameters
  nboot       <- 1000  # Number of bootstrap simulations to calculate C.I.
  alpha       <- 0.1   # Confidence level (0.1 = 90%)
  
  # Index and extract data from large data.frame of CHAR
  # Hard-code smallest binwidth, so that vectors are the same for both 100-yr and 500-yr curves
  ageLim.ind <- which(charData[,'age']>=(ageLim[1]-50) & charData[,'age']<=(ageLim[2]+50))
  
  x <- charData[,'age'][ageLim.ind] 
  y <- charData[,'char'][ageLim.ind]
  
  site <- charData[,'lake'] 
  site <- site[ageLim.ind]
  
  region <- charData[,'region']
  region <- region[ageLim.ind]
  
  sites <- lakeList
  n.sites <- length(lakeList)
  n.regions <- length(unique(regions[,'region']))
  
  # X vals to calculate at, X vector lengths
  xx <- seq(ageLim[1],ageLim[2],xxStep)
  n <- length(y)
  nn <- length(xx)
  
  # Set a tiny number for approximating strictly >0 in optim bounds below. In other words, setting a bound of 0 on one of optim's parameters allows the value 0 to be used, which will kill the optimization if f(0) is undefined. Instead use a number very near 0.
  near0 <- 10^-10
  #------------------------------------------------------------------------
  
  
  # Determine number of records and regions contributing at each timestep
  #------------------------------------------------------------------------
  # Records
  site.ends = data.frame(min=NA,max=NA,site=NA,region=NA)
  for(i in 1:n.sites) {
    x.site = x[site==sites[i]]
    site.ends[i,1:2] = c(min(x.site),max(x.site))
    site.ends[i,3] = sites[i]
    site.ends[i,4] = regions[regions[,1] == sites[i],2] 
  }
  
  n.recs.x = numeric(n)
  for(i in 1:n) {
    n.recs.x[i] = sum(x[i]>=site.ends[,1] & x[i]<=site.ends[,2])
  }
  
  n.recs.xx = numeric(nn)
  for( i in 1:nn) {
    n.recs.xx[i] = sum(xx[i]>=site.ends[,1] & xx[i]<=site.ends[,2])
  }
  
  region.ends = data.frame(min=NA,max=NA,region=NA)
  for( i in 1:n.regions) {
    x.region = x[region==unique(region)[i]]
    region.ends[i,1:2] = c(min(x.region),max(x.region))
    region.ends[i,3] = unique(region)[i]
  }
  
  # Regions
  n.regions.contributing = numeric(nn);
  for( i in 1:n) {
    n.regions.contributing[i] = sum(x[i]>=region.ends[,1] & x[i]<=region.ends[,2])
  }  
  
  #------------------------------------------------------------------------
  
  # Weight CHAR so that all regions contribute equally to the composite record regardless of how many records are in the region 
  #------------------------------------------------------------------------
  if(all(lakeList == alaska)){
    # Save original y's
    # y.uW <- y
    
    # Weight CHAR so that all regions contribute equally to the mean MLE
    w.char = rep(NA,n)
    y.w = rep(NA,n)
    
    for(i in 1:n.regions){
      region.i = which(region == unique(region)[i])
      for(j in 1:nn){
        ind.x  = which( x[region.i] >= (xx[j] - bandWidth) &
                          x[region.i] <= (xx[j] + bandWidth) )
        
        reg.n.sites = length(unique(site[region.i][ind.x]))
        w.char[region.i][ind.x] = 1 / (reg.n.sites * median(n.regions.contributing[region.i][ind.x], na.rm=T)) #
      }
    }
    # Mulitply CHAR by the weight, rescale by multiplying again by the number of total records at that timestep
    y = ( y * w.char  *  n.recs.x ) 
  }
  #------------------------------------------------------------------------
  
  
  # Maxium-likelihood estimatation of ZIL distribution within each window at each timestep.
  #------------------------------------------------------------------------
  # Compute Local Likelihoods
  # Init ouputs
  fitparms <- matrix(data=NA, nrow=nn, ncol=3)
  n.eff		 <- numeric(nn)
  meanMLE  <- numeric(nn)
  medMLE   <- numeric(nn)
  bad.convergence.count <- 0
  kernel.name	<- 'gauss'      # Distribution kernel to use (a function from kernels.r)
  eval(parse(text=paste('kernel.fn = kernel.',kernel.name,sep="")))
  
  # Run local likelihood at every requested X
  for( i in 1:nn) { 
    # Find weights
    x0 = xx[i] # Center-point for kernel
    w = numeric(n)		# Initialize weight vector
    
    for(j in 1:n.sites) { # j=1
      ind.j = which(site==sites[j]) # Index to site j only
      
      # Compute raw kernel weights
      kern.j = ry.kernel(x[ind.j],x0,bandWidth,kernel.fn)			
      k.j = kern.j$K;
      
      if(sum(k.j)==0) { # There are no samples in record j with w>0 for this x0
        w[ind.j] = 0
      } else { # Otherwise,
        # Weight the site by its sample resolution in the kernel around x0, defined by the number of samples with nonzero weight, divided by the total kernel density for this record centered on x0 (may not be 1 because the kernel may be trailing off one end of the record)
        kern.res.j = sum(k.j>0) /
          sum( ry.kernel(min(x[ind.j]):max(x[ind.j]),x0,bandWidth,kernel.fn)$K )	
        
        w[ind.j] = k.j/kern.res.j
      }
    }
    
    # Reweight one last time so sum(w)=1. Might not really be necessary...
    w = w/sum(w)
    
    # Find  n.eff (effective df for the fit dist). This is equivalent to summing up the n.eff for each site j separately, but a little more straightforward.
    kern.out  = ry.kernel(x,x0,bandWidth,kernel.fn) 
    n.eff[i] = kern.out$n # This is (supposed to be) analagous to n in a non-weighted sample
    
    # Fit parameters
    if(all(y[w>0] == 0)) {
      # If all y with nonzero weight are equal to zero, p=1 and mu,sd can't be estimated
      fitparms[i,] = c(1,NA,NA);
    } else {
      # Otherwise, fit distribution using local likelihood. First, find initial parameters from unweighted analytical MLE.
      if(all(y[w>quantile(w[w>0],0.8)] == 0)) {
        # Prefer to use the points with greatest weight (see 'else' below), but if all these are equal to zero, then use all the points with nonzero weight. Some of those must be positive or we wouldn't have gotten past the first 'if' above.
        init.parms = ZIlnorm.mle(y[w>0]);
      } else {
        # In most cases, set initial params as the unweighted analytical MLE of the upper quantile of the values with non-zero weight.
        init.parms = ZIlnorm.mle(y[w>quantile(w[w>0],0.8)])
      }
      
      # Finally, if sd param came pack as zero (can happen if there's only one non-zero point with non-zero weight), replace it with a small fraction of the mean.
      if(init.parms[3]==0) {init.parms[3] =0.01*exp(init.parms[2])}
      
      # Now that the inits are set, calculate parameter scale values to pass to optim. Approximate the partial derivatives of  Zilnorm.llk near the initial parameters, and use those.
      delta1=0.01
      pscale1 = abs(
        (ZIlnorm.llk2(init.parms[2:3],y,w) - ZIlnorm.llk2((init.parms[2:3]+c(delta1,0)),y,w)) / 
          delta1)
      delta2=0.01
      pscale2 = abs(
        (ZIlnorm.llk2(init.parms[2:3],y,w) - ZIlnorm.llk2((init.parms[2:3]+c(0,delta2)),y,w)) /
          delta2)
      
      fit = optim(par=init.parms[2:3], fn=ZIlnorm.llk2, x=y, w=w,
                  method="L-BFGS-B",lower=c(-Inf,near0),upper=c(Inf,Inf),
                  control=list(parscale=c(pscale1,pscale2)))
      
      if(fit$convergence != 0) {
        # If optim didn't converge, then leave mu,sd as NA and print a warning.
        # print(paste(i,fit$message));
        bad.convergence.count = bad.convergence.count+1
      } else {
        # Assuming convergence, store the fitted values.
        fitparms[i,2:3] = fit$par
      }
      
      # Finally, store the estimated of p (which is easy)
      fitparms[i,1] = sum((y==0)*w)
    }
    
    # Compute and save the mean of the distribution specified by the fitted parameters
    meanMLE[i] = ZIlnorm.mean(fitparms[i,]) 
    medMLE[i] = ZIlnorm.median(fitparms[i,])
    # print(paste(i,"of",nn)) Turn on to see progress messages
  }
  #------------------------------------------------------------------------
  
  #  Bootstrap CI
  #------------------------------------------------------------------------
  boot.means = matrix(NA,nrow=nn,ncol=nboot)
  
  for( i in 1:nn ) { #i=1
    parm = fitparms[i,]
    n.eff.i = floor(n.eff[i])
    
    if(n.eff.i>0) {
      boot.dat = matrix(data=ZIlnorm.randomdraw(parm,n.eff.i*nboot),nrow=n.eff.i,ncol=nboot)
      boot.means[i,] = apply(boot.dat,2,mean)
    } else {
      boot.means[i,] = rep(NA,nboot)
    }
    
   # print(paste(i,"of",nn)) Turn this back on to print progress message
  }
  
  out.ci.L = apply(boot.means,1,quantile,probs=alpha/2,na.rm=TRUE)
  out.ci.U = apply(boot.means,1,quantile,probs=1-alpha/2,na.rm=TRUE)
  out.mean.se = apply(boot.means,1,sd,na.rm=T)
  #------------------------------------------------------------------------
  
  
  # Summarize output
  #------------------------------------------------------------------------
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
  #------------------------------------------------------------------------

  return(output)
}


  
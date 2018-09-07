# Set up data for ZIL scripts adapted from Kelly et al. 2013
# Turn multiple dataframes into one long data frame
dat <- do.call(rbind,lapply(dataframes,data.frame,stringsAsFactors=TRUE))
dat <- dat[,c('char','age','lake','region')]

ageLim.ind <- which(dat[,'age']>=(ageLim[1]-2*bandWidth) & dat[,'age']<=(ageLim[2]+2*bandWidth))
x <- dat[,'age'][ageLim.ind] 
y <- dat[,'char'][ageLim.ind]

site <- dat[,'lake'] 
site <- site[ageLim.ind]
site <- site[!is.na(y)]

region <- dat[,'region']
region <- region[ageLim.ind]
region <- region[!is.na(y)]

x    <- x[!is.na(y)]
y    <- y[!is.na(y)]

# Get unique site names
sites <- lakeList
n.sites <- length(sites)

# X vals to calculate at, X vector lengths
xx <- seq(ageLim[1],ageLim[2],xxStep)
n <- length(y)
nn <- length(xx)

n.sites.contributing = NA
n.regions.contributing = NA

for(i in 1:length(xx)){
  n.sites.contributing[i] <- length(unique(dat[ dat[,'age'] >= xx[i] & dat[,'age'] < xx[i+1],'lake']))
  n.regions.contributing[i] <- length(unique(dat[ dat[,'age'] > xx[i] & dat[,'age'] < xx[i+1],'region']))
}

n.records <- read.csv('data/n.sites.regions.csv')
n.sites.contributing <- n.records[,'sites']
n.regions.contributing <- n.records[,'regions']

w.char = 1/(n.regions.contributing * n.sites.contributing)

# Set a tiny number for approximating strictly >0 in optim bounds below. In other words, setting a bound of 0 on one of optim's parameters allows the value 0 to be used, which will kill the optimization if f(0) is undefined. Instead use a number very near 0.
near0 <- 10^-10

# Follow procedure from Kelly et al. 2013
source('code/ZIL_Fn_ZIlnorm.r')
source('code/ZIL_Fn_Kernels.r')
# Compute Local Likelihoods
# Init ouputs
fitparms <- matrix(data=NA, nrow=nn, ncol=3)
n.eff		 <- numeric(nn)
meanMLE  <- numeric(nn)
medMLE   <- numeric(nn)
bad.convergence.count <- 0

# Run local likelihood at every requested X
eval(parse(text=paste('kernel.fn = kernel.',kernel.name,sep="")))


#TH: Now adjust y (CHAR) at each site, at each timestep, for it's regionally-weighted contribution
for (i in 1:n){
  x.ind = which(x < x0 + bandWidth & x > x0 - bandWidth)
  if(y > 0) {
    y = y * w.char
  }
}

for( i in 1:nn) { # i=1
  # Find weights
  x0 = xx[i] # Center-point for kernel
  w = numeric(n)		# Initialize weight vector
  
  for(j in 1:n.sites) { # j=1
    ind.j = which(site==sites[j]) # Index to site j only
    
    # Compute raw kernel weights
    kern.j = ry.kernel(x[ind.j],x0,bandWidth,kernel.fn)			
    k.j = kern.j$K;
    
    if(sum(k.j)==0) { # There are no samples in record j with w>0 for this x0
      w[ind.j] = 0;
    } else { # Otherwise,
      # Weight the site by its sample resolution in the kernel around x0, defined by the number of samples with nonzero weight, divided by the total kernel density for this record centered on x0 (may not be 1 because the kernel may be trailing off one end of the record)
      kern.res.j = sum(k.j>0) /
        sum( ry.kernel(min(x[ind.j]):max(x[ind.j]),x0,bandWidth,kernel.fn)$K )	
      
      w[ind.j] = k.j/kern.res.j
    }
  }
  
  # Reweight one last time so sum(w)=1. Might not really be necessary...
  w = w/sum(w)
  w.char = w.char/sum(w.char)
  
  # Find  n.eff (effective df for the fit dist). This is equivalent to summing up the n.eff for each site j separately, but a little more straightforward.
  kern.out  = ry.kernel(x,x0,bandWidth,kernel.fn) 
  n.eff[i] = kern.out$n # This is (supposed to be) analagous to n in a non-weighted sample
  
  # Fit parameters
  if(all(y[w>0] == 0)) {
    # If all y with nonzero weight are equal to zero, p=1 and mu,sd can't be estimated
    fitparms[i,] = c(1,NA,NA)
  } else {
    # Otherwise, fit distribution using local likelihood. First, find initial parameters from unweighted analytical MLE.
    if(all(y[w>quantile(w[w>0],0.8)] == 0)) {
      # Prefer to use the points with greatest weight (see 'else' below), but if all these are equal to zero, then use all the points with nonzero weight. Some of those must be positive or we wouldn't have gotten past the first 'if' above.
      init.parms = ZIlnorm.mle(y[w>0])
    } else {
      # In most cases, set initial params as the unweighted analytical MLE of the upper quantile of the values with non-zero weight.
      init.parms = ZIlnorm.mle(y[w>quantile(w[w>0],0.8)])
    }
    
    # Finally, if sd param came pack as zero (can happen if there's only one non-zero point with non-zero weight), replace it with a small fraction of the mean.
    if(init.parms[3]==0) {init.parms[3] =0.01*exp(init.parms[2])}
    
    # Now that the inits are set, calculate parameter scale values to pass to optim. Approximate the partial derivatives of  Zilnorm.llk near the initial parameters, and use those.
    delta1=1
    pscale1 = abs(
      (ZIlnorm.llk2(init.parms[2:3],y,w) - ZIlnorm.llk2((init.parms[2:3]+c(delta1,0)),y,w)) / 
        delta1)
    delta2=1
    pscale2 = abs(
      (ZIlnorm.llk2(init.parms[2:3],y,w) - ZIlnorm.llk2((init.parms[2:3]+c(0,delta2)),y,w)) /
        delta2)
    
    fit = optim(
      par=init.parms[2:3], fn=ZIlnorm.llk2, x=y, w=w,
      method="L-BFGS-B",lower=c(-Inf,near0),upper=c(Inf,Inf),
      control=list(parscale=c(pscale1,pscale2)))
    
    if(fit$convergence != 0) {
      # If optim didn't converge, then leave mu,sd as NA and print a warning.
      print(paste(i,fit$message));
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
  print(paste(i,"of",nn))
}


# -------------------- Summarize output
output = data.frame(
  yr.bp     = xx,
  composite.mean = meanMLE,
  composite.med  = medMLE,
  n.records			 = n.recs.contributing,
  p.hat     = fitparms[,1],
  mu.hat    = fitparms[,2],
  sigma.hat = fitparms[,3],
  n.eff     = n.eff
)


# -------------------- Bootstrap CI
source('code/ZIL_BootstrapCI.r')
output$CI.lower       = out.ci.L
output$CI.upper       = out.ci.U		
output$mean.se        = out.mean.se


# -------------------- Save data
# write.csv(output,paste0(output.dir,'ZIL_',reg.name,'_win',bandWidth,'_step',xxStep,iter,'.csv'),
#           row.names=FALSE,na="nan")


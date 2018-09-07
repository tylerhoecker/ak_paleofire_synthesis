# Bootstrap mean with CI for time series of parameters to zero-inflated lognormal distributions

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
	
	print(paste(i,"of",nn))
}

out.ci.L = apply(boot.means,1,quantile,probs=alpha/2,na.rm=TRUE)
out.ci.U = apply(boot.means,1,quantile,probs=1-alpha/2,na.rm=TRUE)
out.mean.se = apply(boot.means,1,sd,na.rm=T)


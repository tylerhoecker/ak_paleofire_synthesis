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


#----------------------------- OLD JUNK
# *** This one works on a whole vector at once. But it can't easily incorporate kernel weights as it's written.
# Log-likelihood for zero-inflated lognormal
#ZIlnorm.llk = function(parms,sampledata) {
#	n0 = sum(sampledata==0)
#	n1 = length(sampledata) - n0
#	return(-1*(
#		n0*log(parms[1]) +
#		n1*log(1-parms[1]) +
#		sum(dlnorm(sampledata[sampledata>0],parms[2],parms[3], log=TRUE))
#	))
#}

# This was a practice, for a single value
# Log-likelihood for zero-inflated lognormal
#ZIlnorm.llk.w = function(parms,x,w) {
#	if( missing(w) ) { w=1 }
#
#	if( x==0 ) {
#		likelihood = log(parms[1])
#	} else {
#		likelihood = log(1-parms[1]) + dlnorm(x,parms[2],parms[3],log=TRUE)
#	}
#	
#	return(-w*likelihood)
#}
#

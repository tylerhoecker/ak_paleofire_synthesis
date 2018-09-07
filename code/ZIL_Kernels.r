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



# Save original y's
# y.uW <- y

# Weight CHAR so that all regions contribute equally to the mean MLE
# Mulitply CHAR by the weight, rescale by multiplying again by the number of total records at that timestep
w.char = rep(NA,n)
y.w = rep(NA,n)

for(i in 1:n.regions){
  region.i = which(region == unique(region)[i])
  for(j in 1:nn){
    ind.x  = which( x[region.i] >= (xx[j] - bandWidth) &
                      x[region.i] <= (xx[j] + bandWidth) )
    
    reg.n.sites = length(unique(site[region.i][ind.x]))
    
    if(reg.n.sites > 1 & median(n.regions.contributing[region.i][ind.x]) > 1){
      w.char[region.i][ind.x] =
        1/( reg.n.sites *
              mean(n.regions.contributing[region.i][ind.x], na.rm=T) )
      y.w[region.i][ind.x] = 
        ( y[region.i][ind.x] * w.char[region.i][ind.x] ) *  n.recs.x[region.i][ind.x]
    }else{
      y.w[region.i][ind.x] = NA
    }
  }
}

# Now need to remove some NAs again...
y    <- y.w[!is.na(y.w)]
x    <- x[!is.na(y.w)]
site <- site[!is.na(y.w)]
region <- region[!is.na(y.w)]
n <- length(y)



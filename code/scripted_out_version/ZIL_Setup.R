# Set up data for ZIL procedure
# Temporal params
xxStep      <- 10				    # Interval, in years, for moving window

if(all(lakeList == copper | lakeList == noatak)){
  ageLim = c(-60,7000)
}else{
  ageLim = c(-60,10000)}

# Confidence Intervals parameters
nboot       <- 1000  # Number of bootstrap simulations to calculate C.I.
alpha       <- 0.1   # Confidence level (0.1 = 90%)

# Index and extract data from large data.frame of CHAR
ageLim.ind <- which(charData[,'age']>=(ageLim[1]-bandWidth) & charData[,'age']<=(ageLim[2]+bandWidth))

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

return()

# Determine number of records and regions contributing at each timestep
site.ends = data.frame(min=NA,max=NA,site=NA,region=NA)
for(i in 1:n.sites) {
  x.site = x[site==sites[i]]
  site.ends[i,1:2] = c(min(x.site),max(x.site))
  site.ends[i,3] = sites[i]
  site.ends[i,4] = regions[regions[,1] == sites[i],2] 
}

n.recs.x = numeric(n);
for(i in 1:n) {
  n.recs.x[i] = sum(x[i]>=site.ends[,1] & x[i]<=site.ends[,2])
}

n.recs.xx = numeric(nn);
for( i in 1:nn) {
  n.recs.xx[i] = sum(xx[i]>=site.ends[,1] & xx[i]<=site.ends[,2])
}

# Determine number of records and regions contributing at each timestep
region.ends = data.frame(min=NA,max=NA,region=NA)
for( i in 1:n.regions) {
  x.region = x[region==unique(region)[i]]
  region.ends[i,1:2] = c(min(x.region),max(x.region))
  region.ends[i,3] = unique(region)[i]
}

n.regions.contributing = numeric(nn);
for( i in 1:n) {
  n.regions.contributing[i] = sum(x[i]>=region.ends[,1] & x[i]<=region.ends[,2])
}  




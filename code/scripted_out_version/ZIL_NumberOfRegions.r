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
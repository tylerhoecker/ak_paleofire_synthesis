# Calculate the number of records at any time

site.ends = data.frame(min=NA,max=NA,site=NA,region=NA)
for( i in 1:n.sites) {
	x.site = x[site==sites[i]]
	site.ends[i,1:2] = c(min(x.site),max(x.site))
	site.ends[i,3] = sites[i]
	site.ends[i,4] = regions[regions[,1] == sites[i],2] 
}

n.recs.x = numeric(n);
for( i in 1:n) {
	n.recs.x[i] = sum(x[i]>=site.ends[,1] & x[i]<=site.ends[,2])
}

n.recs.xx = numeric(nn);
for( i in 1:nn) {
	n.recs.xx[i] = sum(xx[i]>=site.ends[,1] & xx[i]<=site.ends[,2])
}


# n.recs.contributing = data.frame(year=NA,sites=NA,regions=NA)
# for( i in 1:nn) {
#   n.recs.contributing[i,1] = xx[i]
#   n.recs.contributing[i,2] = sum(xx[i]>=site.ends[,1] & xx[i]<=site.ends[,2])
#   n.recs.contributing[i,3] = sum(x[i]>=region.ends[,1] & x[i]<=region.ends[,2]) 
# }


## This script only needs to be run one time to produce the 

# Noatak River Basin (data from Higuera et al. 2011)
noatak <- c('LI','RA','PO','UC')  
# Brooks Range (data from Higuera et al. 2007,2009)
brooks <- c('CO','LC','RP','WK','XI') 
# Yukon Flats (data from Kelly et al. 2014)
yukon  <- c('CP','EP','GA','JA','LD','LT','LU','NR','PI','RE','RO','SL','WC','WI') 
# Copper River Basin (data from Barett et al. 2013)
copper <- c('CR','HD','MP1','SC') 

alaska <- c(noatak,brooks,yukon,copper)
 
regions  <- list('noatak'= noatak,'brooks'= brooks,'yukon'= yukon, 'copper'= copper) %>% 
  stack() %>% 
  rename(lake = values, region = ind)
# Alaska-wide "region"
alaska <- c(noatak,brooks,yukon,copper)

# Import raw charcoal count data, standardize and derive CHAR
#------------------------------------------------------------------------

# Function for standardizing non-zero charcoal accumlations rates
trans_fn <- function(x) {
  x = ifelse(x > 0, x, NA)
  logX = log(x)
  zX = (logX - mean(logX, na.rm = T)) / sd(logX, na.rm = T)
  expX = exp(zX)
  expX[is.na(expX)] = 0
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

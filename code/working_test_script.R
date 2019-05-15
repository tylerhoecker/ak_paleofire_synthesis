# Noatak River Basin (data from Higuera et al. 2011)
noatak <- c('LI','RA','PO','UC')  
# Brooks Range (data from Higuera et al. 2007,2009)
brooks <- c('CO','LC','RP','WK','XI') 
# Yukon Flats (data from Kelly et al. 2014)
yukon  <- c('CP','EP','GA','JA','LD','LT','LU','NR','PI','RE','RO','SL','WC','WI') 
# Copper River Basin (data from Barett et al. 2013)
copper <- c('CR','HD','MP1','SC') 

alaska <- c(noatak,brooks,yukon,copper)

source(file.path('code','ZIL_Run.r'))

bandWidth   <- 50 
output100 <- map(list(alaska), #, noatak, brooks, yukon, copper
                 zil_fn)

compCHAR_100 <- output100 %>% 
  #`names<-` (c('Alaska','Noatak','Kobuk','Yukon Flats','Copper')) %>%
  bind_rows(.id = 'region')

bandWidth   <- 250 	
output500 <- map(list(alaska), #, noatak, brooks, yukon, copper
                 zil_fn) 

compCHAR_500 <- output500 %>% 
  #`names<-` (c('Alaska','Noatak','Kobuk','Yukon Flats','Copper')) %>%
  bind_rows(.id = 'region')


compositeCHAR <- full_join(compCHAR_500, compCHAR_100, 
                           by = c('region', 'yr.bp'),
                           suffix = c('.500','.100'))

ggplot() +
  geom_line(data = compCHAR_100, aes(x = yr.bp, y = composite.mean)) +
  geom_ribbon(data = compCHAR_500, aes(x = yr.bp, ymin = CI.lower, ymax = CI.upper), alpha = 0.3) +
  geom_line(data = compCHAR_500, aes(x = yr.bp, y = composite.mean), size = 1) +
  facet_wrap(~region) +
  scale_x_reverse() +
  coord_cartesian(ylim = c(0,3.5)) +
  theme_bw() +
  labs(x = 'Time (cal years BP)', y = 'Biomass burning (standardized CHAR)')




bandWidth   <- 10 	
output20 <- map(list(alaska, noatak, brooks, yukon, copper),
                zil_fn) 

compCHAR_20 <- output20 %>% 
  `names<-` (c('Alaska','Noatak','Kobuk','Yukon Flats','Copper')) %>%
  bind_rows(.id = 'region')

ggplot(compCHAR_100, aes(x = yr.bp, y = composite.mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), alpha = 0.3) +
  geom_smooth(fill = 'blue', alpha = 0.3) +
  facet_wrap(~region) +
  scale_x_reverse() +
  coord_cartesian(ylim = c(0,3.5)) +
  theme_bw() +
  labs(x = 'Time (cal years BP)', y = 'Biomass burning (standardized CHAR)')





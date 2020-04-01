# Define list of regions
regions  <- list('noatak'= noatak,'brooks'= brooks,'yukon'= yukon, 'copper'= copper) %>% 
  stack() %>% 
  rename(lake = values, region = ind)
# Alaska-wide "region"
alaska <- c(noatak,brooks,yukon,copper)

lakeList = alaska


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


ylimit <- c(-.1,6) 

rawCHAR <-  
  ggplot() + 
  geom_point(data = sample_n(charData,1000), aes(x = age, y = char), 
             fill = 'grey70', shape = 21, color = 'grey50') +
  geom_line(data = filter(compositeCHAR, region == 'Alaska'), aes(x = yr.bp, y = mean_100), 
            color = 'grey20', size = 1) +
  xlab('Time (cal. yr BP)') +
  xlim(10000,-65) +
  coord_cartesian(ylim = ylimit) +
  theme_bw(base_size = 12) +
  theme(plot.margin=unit(c(.5,.5,.5,0), "cm"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) 


zilHIST = 
  ggplot(sample_n(charData,1000)) + 
  geom_histogram(aes(x = char, y = ..density..), 
                 bins = 50, color = 'black', fill = 'grey70', size = 0.4) +
  
  xlim(ylimit) +
  scale_y_continuous(breaks = c(0,1.0,2.0), labels = c(0,1.0,2.0)) +
  labs(x=bquote('Biomass burning (standardized CHAR)'), y= 'Probability density') +
  coord_flip() +
  theme(plot.margin=unit(c(.5,-0.1,.5,.5), "cm"),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line()) +
  theme_bw(base_size = 12)


library(cowplot)
zil_method_plot <- plot_grid(zilHIST, rawCHAR, ncol = 2, rel_widths = c(1,2), align = 'h')


plot_data <- (filter(compositeCHAR, region != 'Alaska')) %>% 
  ungroup() %>% 
  mutate(region = factor(region, levels = c('Noatak River Watershed','Kobuk Ridges and Valleys','Yukon Flats','Copper River Basin')))
#   gather(key, value, -c(region, yr.bp)) %>% 
#   mutate(value = if_else(value > 4, as.numeric(NA), value),
#          value = if_else(value < 0, as.numeric(NA), value)) %>% 
#   spread(key, value)

ggplot(plot_data) + 
  geom_line(aes(x= yr.bp, y = mean_100), color = 'grey40') +
  geom_ribbon(aes(x=yr.bp, ymin = CI_lower_500, ymax = CI_upper_500), alpha = 0.2) +
  geom_line(aes(x=yr.bp,y=mean_500), size = 1) +
  scale_x_reverse(breaks = seq(10000,-50,-2000)) +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7), labels =c(1,2,3,4,5,6,7), expand = c(0,0)) +
  coord_cartesian(ylim = c(0.1,4)) +
  ylab("Biomass burning (standardized CHAR)") +
  xlab('Time (cal years BP)') +
  facet_wrap(~region, ncol = 1, scales = 'free_y') +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        plot.subtitle = element_text(hjust = 0.5, face = 'bold'))
          
            
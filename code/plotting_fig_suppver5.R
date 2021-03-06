# Biomass burning
#--------------------------------------------------------------------------
ylimit = c(0.3,3.1)
xlimit = c(10000,-60)

char.plot <- 
  ggplot(filter(compositeCHAR, region == 'Alaska')) + 
  geom_line(aes(x= yr.bp, y = mean_100), color = 'grey40') +
  geom_ribbon(aes(x=yr.bp, ymin = CI_lower_500, ymax = CI_upper_500), alpha = 0.2) +
  geom_line(aes(x=yr.bp,y=mean_500), size = 1) +
  scale_x_reverse(breaks = seq(10000,-50,-2000)) +
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7), labels =c(1,2,3,4,5,6,7)) +
  coord_cartesian(ylim = c(0.3,3.1)) +
  ylab("Biomass burning\n(standardized CHAR)") +
  xlab('Time (cal years BP)') +
    theme_bw(base_size = 14) + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          #axis.title.y = element_blank(),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          #axis.line.x = element_line(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin=unit(c(0.5,0.5,0,0.5), "cm"))
#--------------------------------------------------------------------------

# FRI 

# FRI
#--------------------------------------------------------------------------
friPlot <- 
  ggplot(FRIsmooth) + theme_bw(base_size = 14) +
  geom_point(data = fireYears_df, aes(x = year, y = FRI), shape = 15, color = 'grey40') +
  geom_line(aes(x = year_BP, y=,mFRI_1000), size = 1) +
  geom_ribbon(aes(x = year_BP, ymin = CI_l, ymax = CI_u), alpha = 0.2) +
  scale_x_reverse(breaks = seq(10000,-50,-2000)) +
  scale_y_continuous(breaks = c(100,150,200), position = 'right') +
  coord_cartesian(ylim = c(75,200)) +
  xlab('Time (cal years BP)') +
  ylab('Fire return interval (yr)') +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0.5,0,0.5),'cm'))
#--------------------------------------------------------------------------

# Ratio of CHAR:FRI 
#--------------------------------------------------------------------------
ratioPlot <- 
  ggplot(ratio, aes(x = yr.bp, y = ratio, ymin = ratio.CIlower, ymax = ratio.CIupper)) + 
  theme_bw(base_size = 14) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.2) +
  scale_x_reverse(limit = xlimit, breaks = seq(10000,-50,-2000)) +
  scale_y_continuous(breaks = c(0.005, 0.015, 0.025)) +
  coord_cartesian(ylim = c(0.000, 0.025)) +
  scale_y_continuous(breaks = c(0.005, 0.015, 0.025)) +
  xlab('Time (cal years BP)') +
  ylab('CHAR : FRI') +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0.5,0,0.5),'cm'))
#--------------------------------------------------------------------------

# Percent burned
#--------------------------------------------------------------------------
#ylimit = c(10,150)

pctBurned_ak <- filter(pctBurned_ak, year <= 8000)


span.100 = (500/10) / (length(!is.na(pctBurned_ak[,'boot.median']))) 


loessPlot <- ggplot(pctBurned_ak, aes(x= year)) +
  geom_smooth(aes(y=boot.lower), se=F, method='loess', 
              span = span.100) +
  geom_smooth(aes(y=boot.upper), se=F, method='loess', 
              span = span.100)
ggloess <- ggplot_build(loessPlot)
dataloess <- data.frame(x= ggloess$data[[1]]$x,
                        ymin = ggloess$data[[1]]$y,
                        ymax = ggloess$data[[2]]$y)


ak.pct.TS.plot <-   
  ggplot(pctBurned_ak) +
  geom_smooth(aes(x= year, y = boot.median), size = 1,
              color = 'grey20', se = F, method = 'loess', span = span.100) +
  # geom_smooth(data = alaska.df, aes(x= year, y = boot.median), size = 0.75,
  #             color = 'red', se = F, method = 'loess', span = span.100) +
  geom_ribbon(data = dataloess, aes(x=x, ymin=ymin, ymax=ymax),alpha = 0.2) +
  #coord_cartesian(ylim = ylimit) +
  scale_x_reverse(limits = xlimit, breaks = seq(10000,-50,-2000)) +
  ylab('Percent\nsites burned') +
  xlab('Time (cal years BP)') +
  scale_y_continuous(breaks = c(50,100), position = 'right') +
  theme_bw(base_size = 14) +
  theme(plot.margin=unit(c(0,0.5,0,0.5), "cm"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.title.y = element_blank(),
        axis.line.y = element_line(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line()) 
#--------------------------------------------------------------------------

# Climate
#--------------------------------------------------------------------------
moisPlot <-  
  ggplot() +   
  geom_col(data = clim.df, aes(x= yr.bp, y = all.moist, fill = moistsign)) +
  scale_fill_manual(values = c('negative' = 'chocolate4','positive'='chartreuse4')) +
  ylab('Moisture\nanomaly (SD)') +
  scale_x_reverse(breaks = seq(10000,-50,-2000)) +
  coord_cartesian(xlim = xlimit) +
  xlab("Time (cal years BP)") +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0.5,0,0.5), "cm")) 

tempPlot <-  
  ggplot() +   
  geom_col(data = clim.df, aes(x= yr.bp, y = all.temp, fill = tempsign)) +
  geom_smooth(data = goa.df, method = 'loess', span = 0.4, se = F, color = 'black', size = .6,
              aes(x= yearBP, y = zscore)) +
  scale_fill_manual(values = c('negative' = 'dodgerblue4','positive'='red3')) +
  ylab('Temperature\nanomaly (SD)') +
  scale_x_reverse(breaks = seq(10000,-50,-2000)) +
  scale_y_continuous(sec.axis = dup_axis(name = 'Temperature\n anomaly (deg. C)')) +
  coord_cartesian(xlim = xlimit) +
  xlab("Time (cal years BP)") +
  theme_bw(base_size = 14) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0.5,0.5,0.5), "cm")) 
#--------------------------------------------------------------------------

# Plot grid
#-------------------------------------------------------------------------- 
fig_5 <- 
  cowplot::plot_grid(char.plot,friPlot,ratioPlot, ak.pct.TS.plot, moisPlot, tempPlot,
                     ncol =1, align = 'v', rel_heights = c(1,1,.75,.75,.75,1),
                     labels = "AUTO", hjust = -10)
#-------------------------------------------------------------------------- 
# 
ggsave(filename = 'synthesis_ak_stack_supplemental.pdf',
       height = 10,
       width = 7,
       units = 'in',
       dpi = 600)

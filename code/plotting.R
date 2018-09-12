# Biomass burning
#--------------------------------------------------------------------------
ylimit = c(0.1,3.1)
xlimit = c(10000,-50)

char.plot <- 
  ggplot(filter(compositeCHAR, region == 'Alaska')) + 
    geom_line(aes(x= yr.bp, y = composite.mean.100), color = 'grey40') +
    geom_ribbon(aes(x=yr.bp, ymin = CI.lower.500, ymax = CI.upper.500), alpha = 0.2) +
    geom_line(aes(x=yr.bp,y=composite.mean.500), size = 1) +
    scale_x_reverse(breaks = seq(10000,-50,-2000)) +
    scale_y_continuous(breaks = c(1,2,3,4,5,6,7), labels =c(1,2,3,4,5,6,7)) +
    coord_cartesian(ylim = ylimit) +
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
          plot.margin=unit(c(0.3,0.3,0,0.3), "cm"))
#--------------------------------------------------------------------------

# FRI 

# FRI
#--------------------------------------------------------------------------
friPlot <- 
  ggplot(FRIsmooth) + theme_bw(base_size = 14) +
  geom_point(data = fireYears, aes(x = year, y = FRI), shape = 15, color = 'grey40') +
  geom_line(aes(x = year_BP, y=,mFRI_1000), size = 1) +
  geom_ribbon(aes(x = year_BP, ymin = CI_l, ymax = CI_u), alpha = 0.2) +
  scale_x_reverse(limits = xlimit, breaks = seq(10000,-50,-2000)) +
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
        plot.margin=unit(c(-0.1,0.3,0,0.3),'cm'))
#--------------------------------------------------------------------------

# Ratio of CHAR:FRI 
#--------------------------------------------------------------------------
ratioPlot <- 
  ggplot(ratio, aes(x = yr.bp, y = ratio, ymin = ratio.CIlower, ymax = ratio.CIupper)) + 
  theme_bw(base_size = 14) +
  geom_line(size = 1) +
  geom_ribbon(alpha = 0.2) +
  scale_x_reverse(limit = xlimit, breaks = seq(10000,-50,-2000)) +
  scale_y_continuous(breaks = c(0.005, 0.010, 0.015)) +
  xlab('Time (cal years BP)') +
  ylab('CHAR : FRI') +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(-0.1,0.3,0,0.3),'cm'))
#--------------------------------------------------------------------------

# Percent burned
#--------------------------------------------------------------------------
ylimit = c(0,150)

pctBurned_ak <- filter(pctBurned_ak, year <= 8000)

samples = length(!is.na(pctBurned.ak[,'win.pct']))
xspan = max(pctBurned.ak[!is.na(pctBurned.ak[,'win.pct']),'year']) -
  min(pctBurned.ak[!is.na(pctBurned.ak[,'win.pct']),'year'])
span.100 = 100/ (xspan/samples) / samples 


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
  ggplot(pctBurned_ak, ) +
  geom_smooth(aes(x= year, y = boot.median), size = 0.75,
              color = 'grey20', se = F, method = 'loess', span = span.100) +
  # geom_smooth(data = alaska.df, aes(x= year, y = boot.median), size = 0.75,
  #             color = 'red', se = F, method = 'loess', span = span.100) +
  geom_ribbon(data = dataloess, aes(x=x, ymin=ymin, ymax=ymax),alpha = 0.2) +
  coord_cartesian(ylim = ylimit) +
  scale_x_reverse(limits = xlimit, breaks = seq(10000,-50,-2000)) +
  ylab('Percent\nsites burned') +
  xlab('Time (cal years BP)') +
  scale_y_continuous(breaks = seq(ylimit[1],ylimit[2],50), position = 'right') +
  theme_bw(base_size = 14) +
  theme(plot.margin=unit(c(-0.1,0.3,0,0.3), "cm"),
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
tempPlot <-  
  ggplot() +   
  geom_col(data = clim.df, aes(x= yr.bp, y = all.temp, fill = tempsign)) +
  geom_smooth(data = goa.df, method = 'loess', span = 0.4, se = F, color = 'black', size = .6,
              aes(x= yearBP, y = zscore)) +
  scale_fill_manual(values = c('negative' = 'dodgerblue4','positive'='red3')) +
  ylab('Composite temp.\nanomaly (deg. C)') +
  scale_x_reverse(breaks = seq(10000,-50,-2000)) +
  scale_y_continuous(sec.axis = dup_axis(name = 'Growing season temp.\n anomaly (deg. C)')) +
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
        plot.margin=unit(c(-0.1,0.3,0.2,0.3), "cm")) 
#--------------------------------------------------------------------------
    
cowplot::plot_grid(char.plot,friPlot,ratioPlot, ak.pct.TS.plot,
          tempPlot,
          ncol =1, align = 'v', rel_heights = c(1,1,.75,.75,1),
          labels = "AUTO", hjust = -10)

ggsave(filename = 'synthesis_ak_stack.pdf',
       height = 12,
       width = 7,
       units = 'in',
       dpi = 600)

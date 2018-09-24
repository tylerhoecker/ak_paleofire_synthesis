# Regional time series
synchPlots = function(df,xlimit,xaxis,ylabel,top){
  
  df <- filter(df, year <= 8000)
  
  samples = length(!is.na(df[,'boot.median']))
  xspan = max(df$year) -
    min(df$year)
  span.100 = 500 / (xspan/samples) / samples 
  
  
  loessPlot <- ggplot(df, aes(x= year)) +
    geom_smooth(aes(y=boot.lower), se=F, method='loess', 
                span = span.100) +
    geom_smooth(aes(y=boot.upper), se=F, method='loess', 
                span = span.100)
  ggloess <- ggplot_build(loessPlot)
  dataloess <- data.frame(x= ggloess$data[[1]]$x,
                          ymin = ggloess$data[[1]]$y,
                          ymax = ggloess$data[[2]]$y)
  synchplot <-   
    ggplot(df) +
    geom_smooth(aes(x= year, y = boot.median), size = 0.75,
                color = 'grey20', se = F, method = 'loess', span = span.100) +
    
    geom_ribbon(data = dataloess, aes(x=x, ymin=ymin, ymax=ymax),alpha = 0.2) +
    scale_x_reverse(limits = xlimit, breaks = seq(10000,-50,-2000)) +
    coord_cartesian(ylim = c(0,150)) +
    ylab('Percent\nsites burned') +
    xlab('Time (cal years BP)') +
    ggtitle(paste0(ylabel)) +
    theme_bw(base_size = 14) +
    theme(plot.margin=unit(c(0.3,0.3,0.3,0.3), "cm"),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line()) 
}

# Brooks
brooks.df <- pctBurned_reg %>% 
  filter(region == 'brooks')
brooks.synch.p = synchPlots(brooks.df,xlimit,xaxis=F,"Kobuk",top=T)
# Copper
copper.df <- pctBurned_reg %>% 
  filter(region == 'copper') 
copper.synch.p = synchPlots(copper.df,xlimit,xaxis=F,"Copper River",top=F)
# Noatak
noatak.df <- pctBurned_reg %>% 
  filter(region == 'noatak') 
noatak.synch.p = synchPlots(noatak.df,xlimit,xaxis=F,"Noatak",top=F)
# Yukon
yukon.df <- pctBurned_reg %>% 
  filter(region == 'yukon') 
yukon.synch.p = synchPlots(yukon.df,xlimit,xaxis=F,"Yukon",top=F)

fig_S2 <- 
cowplot::plot_grid(noatak.synch.p, brooks.synch.p,yukon.synch.p,copper.synch.p,
          ncol = 1, align = 'v', rel_heights = c(.5,.5,.5,.5,.8))

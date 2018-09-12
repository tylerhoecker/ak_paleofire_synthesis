synch_fn <- function(lakeList){
  
  ## Temporal parameters
  # Sites burned per...century
  window <- 100
  # Full study period
  studyPeriod <- c(-60,10000)
  # Timestep to calculate at
  timeStep <- 10
  # Bootstrap iterations
  n.boot <- 99
  # Radiocarbon age error (1 sigma in years)
  error = 25
  # Confidence level for intervals
  alpha = 0.10
  
  # Read in charcoal peak data
  peakData <- map(file.path('data','peakData', paste0(lakeList,'_peakData.csv')),
                  read_csv) %>%
    `names<-` (lakeList) %>%
    bind_rows(.id = 'lake') 
  
  # Comile fire years and Compute FRI
  fireYears <- peakData %>%
    select(lake, age, pkBool) %>%
    filter(age >= studyPeriod[1], age <= studyPeriod[2]) %>%
    group_by(lake) %>%
    filter(pkBool == 1) %>%
    mutate(FRI = lead(age, 1) - age) %>%
    select(lake, year = age, FRI) %>%
    left_join(.,regions, by = 'lake') 
  
  if (regional == 0){
    fireYears[,'region'] <- 'alaska'
  }
  
  # Calculate sites contributing through time
  #------------------------------------------------------------------------------#
  # Vector of years to calculate over
  timePeriod <- seq(studyPeriod[1],studyPeriod[2],timeStep)
  
  # Sites contributing at each timestep, just import from saved (otherwise must
  # determine length based on CHAR records, because peaks do not correspond to
  # full record length).
  
  site_ends <- read_csv(file.path('data','peakData','site_ends_peak.csv'))
  
  if (regional == 0){ 
    site_ends[,'region'] <- 'alaska'
  }
  
  # Iterate through the years and determine the number of sites in each region
  # This is the least efficient part of the code right now. - TH 3/6/17
  sitesByYear <- 
    lapply(timePeriod, function(i){
      site_ends %>%
        group_by(region) %>%
        filter(min <= i, max >= i) %>%
        summarise(year = i, total = n()) 
    }) %>% 
    bind_rows()
  
  # Make a adjustments for CR portion with SNI <3
  if (regional == 0){ 
    sitesByYear <- sitesByYear %>%
      mutate(total = ifelse(year > 2200 & year < 4500,
                            total-1, total))
  }
  sitesByYear <- sitesByYear %>%
    mutate(total = ifelse(region == 'copper' & year > 2200 & year < 4500,
                          total-1, total)) 
  
  
  # Percentage of sites burned calculation 
  #------------------------------------------------------------------------------#
  if (regional == 1){ 
    pctBurned_output <- fireYears %>%
      mutate(year = plyr::round_any(year,timeStep)) %>%
      full_join(sitesByYear,.) %>%
      mutate(fire = ifelse(!is.na(lake),1,0)) %>%
      group_by(year,region) %>%
      mutate(n.burned = ifelse(total >= 2, sum(fire, na.rm = T), NA)) %>%
      mutate(pct = n.burned/total) %>%
      group_by(region,year) %>%
      summarise(total = mean(total),
                n.burned = mean(n.burned),
                pct = mean(pct)) %>%
      # Window
      mutate(win.total = zoo::rollapply(total, window/timeStep, fill= NA,
                                   FUN = mean, na.rm =T)) %>%
      mutate(win.burn = zoo::rollsum(n.burned, window/timeStep, fill= NA)) %>%
      mutate(win.pct = win.burn/win.total*100)}
  
    
  
  # Percent burned with age uncertainty 
  if (regional == 0){ 
    
    pctBurned_output <- fireYears %>% 
      group_by(lake) %>% 
      bootstrap(n.boot, by_group = T) %>% 
      do(lapply(seq(1:n.boot), function(i){
        mutate(.,year = round(truncnorm::rtruncnorm(1, mean=year, sd=error, 
                                         a = timePeriod[1]))) %>%
          mutate(year = plyr::round_any(year,timeStep)) %>%
          full_join(sitesByYear,.) %>%
          mutate(fire = ifelse(!is.na(lake),1,0)) %>%
          group_by(region,year) %>%
          mutate(n.burned = sum(fire)) %>%
          mutate(pct = ifelse(total > 3, n.burned/total*100, NA)) %>%
          # Window 
          group_by(region) %>%
          mutate(win.total = zoo::rollapply(total, window/timeStep, fill= c(100,200,300),
                                       FUN = mean, na.rm =T)) %>%
          mutate(win.total = ifelse(win.total == 100, 
                                    max(total, na.rm=T), win.total)) %>%
          mutate(win.total = ifelse(win.total == 300, total, win.total)) %>%
          mutate(win.burn = zoo::rollsum(n.burned, window/timeStep, 
                                    fill= c(100,200,300))) %>%
          mutate(win.burn = ifelse(win.burn == 100, n.burned, win.burn)) %>%
          mutate(win.burn = ifelse(win.burn == 300, n.burned, win.burn)) %>%
          # Minimum # of records criterion (4)
          mutate(win.pct = ifelse(win.total > 3, win.burn/win.total*100, NA))}) %>% 
          bind_rows()) 
    pctBurned_output <- pctBurned_output %>% 
        ungroup(replicate) %>% 
        group_by(year) %>%
        summarise(boot.total = mean(total, na.rm = T),
                  boot.burned = mean(n.burned, na.rm = T),
                  boot.win.total = mean(win.total, na.rm = T),
                  boot.win.burn = mean(win.burn, na.rm = T),
                  boot.mean = mean(win.pct, na.rm = T),
                  boot.median = quantile(win.pct, 0.50, na.rm = T),
                  boot.lower = quantile(win.pct, alpha/2, na.rm = T),
                  boot.upper = quantile(win.pct, 1-(alpha/2), na.rm = T))
  }
  
  return(pctBurned_output)
}

  
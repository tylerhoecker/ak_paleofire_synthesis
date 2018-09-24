fireyears_fn <- function(lakeList){
  
  ## Temporal parameters
  studyPeriod <- c(-60,10000)
  
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
    mutate(FRI = as.integer(lead(age, 1) - age)) %>%
    select(lake, year = age, FRI) %>%
    left_join(.,regions, by = 'lake') 

  return(fireYears)
}

synch_fn <- function(fireYears){
  
  # Calculate sites contributing through time
  #------------------------------------------------------------------------------#
  ## Temporal parameters
  # Sites burned per...century
  window <- 100
  # Full study period
  studyPeriod <- c(-60,10000)
  # Timestep to calculate at
  timeStep <- 10
  # Bootstrap iterations
  n.boot <- 1000
  # Radiocarbon age error (1 sigma in years)
  error = 25
  # Confidence level for intervals
  alpha = 0.10
  # Vector of years to calculate over
  timePeriod <- seq(studyPeriod[1],studyPeriod[2],timeStep)
  
  # Sites contributing at each timestep, just import from saved (otherwise must
  # determine length based on CHAR records, because peaks do not correspond to
  # full record length).
  
  site_ends <- read_csv(file.path('data','peakData','site_ends_peak.csv'))
  
  # Group regions together if Alaska-wide analysis
  if (regional == 0){ 
    fireYears[,'region'] <- 'alaska'
    site_ends[,'region'] <- 'alaska'
  }
  
  # Iterate through the years and determine the number of sites in each region
  # This is the least efficient part of the code right now. - TH 3/6/17
  sitesByYear <- 
    lapply(timePeriod, function(i){
      site_ends %>%
        mutate(region = as.factor(region)) %>% 
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
  
  if (regional == 1){ 
    sitesByYear <- sitesByYear %>%
      mutate(total = ifelse(region == 'copper' & year > 2200 & year < 4500,
                            total-1, total)) 
  }
  
  # Calculate percentage of sites burned, by iteratively drawing from a distribution of possible fire dates define by the standard deviation of chronological error (25 yr) 
    ageUncrty <- 
      # Age uncertainty component (at individual lake-level)
      lapply(seq(1:n.boot), function(i){
        fireYears %>% 
          mutate(year = round(rtruncnorm(1, mean = year, sd = error, a = timePeriod[1]))) %>%     
          mutate(year = plyr::round_any(year,timeStep)) %>%
          group_by(lake) %>%
          mutate(FRI = lead(year, 1) - year) %>%
          mutate(FRI = ifelse(FRI < 0, NA, ifelse(FRI > 2000, NA, FRI))) %>%
          full_join(sitesByYear,.) %>%
          mutate(fire = ifelse(!is.na(lake),1,0)) %>%
          group_by(region,year) %>%
          mutate(n.burned = ifelse(total >= 2, sum(fire, na.rm = T), NA)) %>%
          mutate(pct = n.burned/total*100) %>%
          summarise(total = mean(total),
                    n.burned = mean(n.burned),
                    pct = mean(pct)) %>%
          # Moving-window component (100-yr window) (at regional level)
          group_by(region) %>%
          mutate(win.total = rollapply(total, window/timeStep, fill= NA,
                                       FUN = mean, na.rm =T)) %>%
          mutate(win.burn = ifelse(win.total >= 2, 
                                   rollapply(n.burned, window/timeStep, fill = NA,FUN = sum, na.rm = T),
                                   NA)) %>%
          mutate(win.pct = win.burn/win.total*100)}) %>%
      bind_rows(.id = 'rep') %>% 
      ungroup(region) 
    
    # If analysis is regional, group by region and year and summarize replicates
    if (regional == 1){
      pctBurned_output <- ageUncrty %>%
        group_by(region,year) %>%
        summarise(boot.total = mean(total, na.rm = T),
                  boot.burned = mean(n.burned, na.rm = T),
                  boot.win.total = mean(win.total, na.rm = T),
                  boot.win.burn = mean(win.burn, na.rm = T),
                  boot.mean = mean(win.pct, na.rm = T),
                  boot.median = quantile(win.pct, 0.50, na.rm = T),
                  boot.lower = quantile(win.pct, alpha/2, na.rm = T),
                  boot.upper = quantile(win.pct, 1-(alpha/2), na.rm = T))
    }
    
    # If analysis is Alaska-wide, group by year and summarize replicates
    if (regional == 0){
      pctBurned_output <- ageUncrty %>%
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

  
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

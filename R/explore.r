require(tidyverse)
require(magrittr)
require(ggplot2)


d <- read.table('./../sdtab020_clean', header = TRUE)
str(d)

dd <- d %>% mutate(across(everything(), as.numeric)
                   , dose = ceil(ID/192))
str(dd)

dd %>% group_by(ID, TIME, CKDEPI, LOGPRO) %>% 
  summarise(n=n())

require(Hmisc)
Hmisc::describe(dd)



dd %>% group_by(TIME, CKDEPI, LOGPRO) %>% 
  dplyr::filter(CKDEPI %in% c(34,40,50) & LOGPRO %in% c(3,3.2,3.4) & ID <= 192) %>% 
  ggplot(aes(x = DV, y=C_VENT, fill = as.factor(TIME))) +
  geom_point(alpha = 1, shape = 21, color = 'white') +
  facet_grid(LOGPRO ~ CKDEPI)

dd %>% select(ID,CKDEPI,LOGPRO) %>% 
  group_by(ID, CKDEPI, LOGPRO)%>%
  ggplot(aes(x = ID, y = CKDEPI, color = as.factor(LOGPRO))) +
  geom_point()


dd %>% group_by(TIME, CKDEPI, LOGPRO) %>% 
  dplyr::filter(CKDEPI %in% c(10) & LOGPRO %in% c(2.2) & ID == 1) %>% 
  ggplot(aes(x = TIME, y = C_VENT, group=ID)) +
  geom_point()







  
temp %>% nest() %>% 
    mutate(tps = map(data, ~fields::Tps(.[c('CKDEPI','LOGPRO')], .$p))
           , tps.xy = map(tps, ~expand_grid(x = seq(0,160,1), y = seq(2,4.5,.01)))
           , tps.z = map(tps, ~predict(., expand_grid(seq(0,160,1), seq(2,4.5,.01)), lambda = 0.1)[,1])) %>% 
    unnest(c(tps.xy, tps.z)) %>% 
  ungroup() %>% group_by(x, y) %>% 
  filter(tps.z>=0.9) %>% slice(which.min(dose)) %>% 
  ggplot(aes(x = x, y = y, fill = dose)) +
  geom_tile()
  



dd %>% group_by(ID, CKDEPI, LOGPRO, dose) %>% 
  filter(TIME==96) %>% 
  summarise(p = sum(C_VENT>=2)/n()) %>% 
  ungroup() %>% group_by(CKDEPI, LOGPRO) %>% 
  filter(p>=0.9) %>% slice(which.min(dose)) %>% 
  select(CKDEPI, LOGPRO, dose) %>% 
  write_delim('mic2.dat')

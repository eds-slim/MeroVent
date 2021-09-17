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
  summarise(meanc = mean(C_VENT)) %>% 
  filter(TIME == 96) %>% 
  ggplot(aes(x = CKDEPI, y = LOGPRO, fill = meanc)) +
  geom_tile()


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


dd %>% group_by(ID, CKDEPI, LOGPRO, dose) %>% 
  #filter(TIME==96) %>% 
  summarise(p = sum(C_VENT>=8)/n()) %>% 
  ggplot(aes(x = CKDEPI, y = LOGPRO, fill = p)) +
  geom_tile() +
  facet_wrap(~dose) +
  scale_fill_gradient2(low = 'red', mid = 'yellow', high = 'darkgreen', midpoint = .5)

d.plot1 <- dd %>% group_by(ID, CKDEPI, LOGPRO, dose) %>% 
  #filter(TIME==96) %>% 
  filter(CKDEPI %in% c(50,60,70) & LOGPRO %in% c(2.3, 2.4, 2.6)) %>% 
  summarise(p = sum(C_VENT>=2)/n()) %>% 
  ungroup() %>% group_by(CKDEPI, LOGPRO) %>% nest() %>% 
  mutate(tps = map(data, ~fields::Tps(.$dose, .$p, lambda = 0))
         , tps.r = map(tps, ~.$fitted.values)) 

d.plot.2 <- d.plot1 %>% 
  mutate(tps.x = map(tps, ~seq(0,6,.1))
         , tps.y = map(tps, ~predict(., seq(0,6,.1))[,1])) %>% 
  unnest(c(tps.x, tps.y))

d.plot1 %>% 
  unnest(c(tps.r, data)) %>% 
  ggplot(aes(x = dose, y = p)) +
  geom_point() +
  geom_point(aes(y = `tps.r`)) +
  geom_line(data = d.plot.2, inherit.aes = FALSE, aes(x=tps.x, y=tps.y), color = 'orange') +
  ggalt::geom_xspline(spline_shape = -0.5, color = 'darkgreen') +
  facet_grid(CKDEPI ~ LOGPRO, scales = 'free')

tps <- dd %>% group_by(ID, CKDEPI, LOGPRO, dose) %>% 
  #filter(TIME==96) %>% 
  filter(CKDEPI %in% c(60) & LOGPRO %in% c(2.4)) %>% 
  summarise(p = sum(C_VENT>=2)/n()) %$% 
  fields::Tps(.$dose, .$p, lambda = 0)

temp <- dd %>% group_by(ID, CKDEPI, LOGPRO, dose) %>% 
  #filter(TIME==96) %>% 
  summarise(p = sum(C_VENT>=2)/n()) %>% 
  ungroup() %>% group_by(dose)

temp %>% 
  ggplot(aes(x = CKDEPI, y = LOGPRO, fill = p)) +
  geom_tile() +
  facet_wrap(~dose) +
  scale_fill_gradient2(low = 'red', mid = 'yellow', high = 'darkgreen', midpoint = -.5)

  
temp %>% nest() %>% 
  mutate(tps = map(data, ~fields::Tps(.[c('CKDEPI','LOGPRO')], .$p))
         , tps.xy = map(tps, ~expand_grid(x = seq(0,160,1), y = seq(2,4.5,.01)))
         , tps.z = map(tps, ~predict(., expand_grid(seq(0,160,1), seq(2,4.5,.01)), lambda = 0.1)[,1])) %>% 
  unnest(c(tps.xy, tps.z)) %>% 
  ggplot(aes(x = x, y = y, fill = tps.z)) +
  geom_tile() +
  facet_wrap(~dose) +
  scale_fill_gradient2(low = 'red', mid = 'yellow', high = 'blue', midpoint = .5)
  
temp %>% nest() %>% 
    mutate(tps = map(data, ~fields::Tps(.[c('CKDEPI','LOGPRO')], .$p))
           , tps.xy = map(tps, ~expand_grid(x = seq(0,160,1), y = seq(2,4.5,.01)))
           , tps.z = map(tps, ~predict(., expand_grid(seq(0,160,1), seq(2,4.5,.01)), lambda = 0.1)[,1])) %>% 
    unnest(c(tps.xy, tps.z)) %>% 
  ungroup() %>% group_by(x, y) %>% 
  filter(tps.z>=0.9) %>% slice(which.min(dose)) %>% 
  ggplot(aes(x = x, y = y, fill = as.factor(dose))) +
  geom_tile()
  
dd %>% group_by(ID, CKDEPI, LOGPRO, dose) %>% 
  filter(TIME==96) %>% 
  summarise(p = sum(C_VENT>=2)/n()) %>% 
  ungroup() %>% group_by(CKDEPI, LOGPRO) %>% 
  filter(p>=0.9) %>% slice(which.min(dose)) %>% 
  ggplot(aes(x = CKDEPI, y = LOGPRO, fill = as.factor(dose))) +
  geom_tile()


dd %>% group_by(ID, CKDEPI, LOGPRO, dose) %>% 
  filter(TIME==96) %>% 
  summarise(p = sum(C_VENT>=2)/n()) %>% 
  ungroup() %>% group_by(CKDEPI, LOGPRO) %>% 
  filter(p>=0.9) %>% slice(which.min(dose)) %>% 
  select(CKDEPI, LOGPRO, dose) %>% 
  write_delim('mic2.dat')

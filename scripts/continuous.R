setwd("~/Desktop/TPC-a-palooza-main/")

library(tidyverse)
library(DBI)
library(RSQLite)
library(janitor)
library(mgcv)
library(gratia)
library(nls.multstart)
library(broom)
library(minpack.lm)
library(rTPC)
library(purrr)
library(ggrepel)
library(cowplot)
################################################### DATA SETUP ################# #####
files <- list.files(path = "./data/Continuous", pattern = "*PM.csv", full.names = TRUE)
meta<-read_csv("./data/Continuous/metadata_continuous.csv")%>%clean_names()%>%filter(is.na(notes))%>%
  select(sample_id,chamber,chamber_volume_m_l)

stepwise_results<-read_csv("./output/Stepwise/Stepwise_TPC_metrics.csv")%>%
  clean_names()%>%
  select(sample_id,trait,estimate)%>%
  filter(trait=="Topt_C")%>%
  select(-trait)%>%
  rename(topt_stepwise=estimate)%>%
  mutate(id=parse_number(sample_id))

################################################### PROCESSING ################# #####
bucket<-list()
for (f in files){
  bucket[[f]]<-read_csv(f)%>%
    clean_names()%>%
    select(date,temperature,oxygen)%>%
    mutate(date=dmy_hms(date),
           file=f,
           chamber = parse_number(str_extract(file, "Ch\\d+(?=_)")))
}

raw_output<-bind_rows(bucket)%>%
  inner_join(.,meta,by="chamber")%>%
  select(-file)%>%
  filter(temperature>19)
  #filter(date > ymd_hms("2026-10-15 12:40:00", tz = tz(date)))

ggplot(raw_output)+
  geom_point(aes(date,oxygen))+
  facet_wrap(~sample_id)

################################################### BLANK FIT ################## #####
blank_data<-raw_output%>%filter(sample_id=="blank")
blank_gam <- gam(oxygen ~ s(temperature,k=5),data = blank_data)
ref_level <- mean(fitted(blank_gam))

################################################### CORRECTED FITS ############# #####
working<-raw_output%>%
  mutate(oxygen_c = oxygen - predict(blank_gam, newdata = cur_data()) + ref_level)%>%
  group_by(sample_id)%>%
  nest()%>%
  #mutate(data = map(data, ~ mutate(.x, well = droplevels(factor(well))))) %>%
  mutate(model = map(data, ~ gam(oxygen_c ~ s(temperature, k = 5), data = .x)),
         derivs = map2(model, data, ~ derivatives(.x, data = .y, type = "central", select = "s(temperature)")),
         data = map2(data, derivs, ~ .x %>% mutate(slope = .y$.derivative, slope_se = .y$.se)),
         data = map2(data, model, ~ .x %>% mutate(fitted = predict(.y, newdata = .x))))%>%
  select(sample_id,data) %>%
  unnest(data)%>%
  mutate(rate=slope*-1,
         lower = rate - 1.96 * slope_se,
         upper = rate + 1.96 * slope_se)
  
ggplot(working)+
  geom_point(aes(temperature,oxygen),alpha=0.005)+
  geom_point(aes(temperature,oxygen_c))+
  geom_line(aes(temperature,fitted),color="blue")+
  facet_wrap(~sample_id,scales="free_y")+
  theme_bw(base_size=8)+
  scale_x_continuous(breaks=seq(15,40,2))+
  ylab("Oxygen Concentration")+
  xlab("Temperature")
quartz()
################################################### TPC ######################## #####
working_tpc<-working%>%rename(temp=temperature)%>%filter(sample_id!="blank")%>%
  group_by(sample_id)%>%
  arrange(sample_id,date)

ggplot(working_tpc)+geom_point(aes(temp,rate))

aic_comparison<-data.frame(model=character(),aic=numeric())
for (i in c(1:49)){
  mod<-get_model_names()[i]
  params<-  setdiff(names(formals(get(mod, asNamespace("rTPC")))), "temp");params
  start_vals <- get_start_vals(working_tpc$temp, working_tpc$rate, model_name = mod)
  low_lims <- get_lower_lims(working_tpc$temp, working_tpc$rate, model_name = mod)
  upper_lims <- get_upper_lims(working_tpc$temp, working_tpc$rate, model_name = mod)
  formula <- as.formula(paste("rate ~",sprintf("%s(temp = temp, %s)",mod,paste(params, collapse = ", "))));formula
  
  fit <- try(nls_multstart(formula,
                           data = working_tpc,
                           iter = 500,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = low_lims,
                           upper = upper_lims,
                           supp_errors = 'Y'),    silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit)) {
    next
  }
  
  aic_comparison <- aic_comparison %>%
    bind_rows(tibble(model = mod, aic = AIC(fit))) %>%
    arrange(aic) %>%
    distinct()
};beepr::beep()

mod<-aic_comparison%>%slice_min(aic)%>%pull(model);head(aic_comparison)
mod<-get_model_names()[40]
params <- setdiff(names(formals(get(mod, asNamespace("rTPC")))), c("temp", "tref"))
start_vals <- get_start_vals(working_tpc$temp, working_tpc$rate, model_name = mod)
low_lims <- get_lower_lims(working_tpc$temp, working_tpc$rate, model_name = mod)
upper_lims <- get_upper_lims(working_tpc$temp, working_tpc$rate, model_name = mod)
formula <- as.formula(paste("rate ~",sprintf("%s(temp = temp, %s)",mod,paste(params, collapse = ", "))));formula

fits <- working_tpc %>%
  group_by(sample_id) %>%
  group_nest() %>%
  mutate(
    start = map(data, ~ get_start_vals(.x$temp, .x$rate, model_name = mod)),
    low   = map(data, ~ get_lower_lims(.x$temp, .x$rate, model_name = mod)),
    high  = map(data, ~ get_upper_lims(.x$temp, .x$rate, model_name = mod)),
    formula = list(as.formula(paste("rate ~", sprintf("%s(temp = temp, %s)", mod, paste(params, collapse = ", "))))),
    fit = pmap(list(data, start, low, high, formula),
               ~ nls_multstart(..5,
                               data = ..1,
                               iter = 500,
                               start_lower = ..2 - 10,
                               start_upper = ..2 + 10,
                               lower = ..3,
                               upper = ..4,
                               supp_errors = "Y",
                               lhstype = "random")),
    params_all = map(fit, ~ calc_params(.x) %>% as_tibble()),
    topt_obs = map_dbl(params_all, "topt"))

param_table <- fits %>%
  select(sample_id, params_all)%>%
  unnest(params_all)%>%
  mutate(adj=topt+2)%>%
  mutate(id=parse_number(sample_id));param_table

predicted <-fits %>%mutate(newdata = map(data, ~ tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 500))),
                                 aug = map2(fit, newdata, ~ augment(.x, newdata = .y))) %>%
  select(sample_id, aug) %>%
  unnest(aug)%>%rename(rate=3)

a<-ggplot(working_tpc)+
  geom_hline(yintercept=0,color="red")+
  geom_point(aes(temp,rate),alpha=0.005)+
  geom_line(aes(temp,rate),data=predicted)+
  geom_point(aes(topt,rmax,fill=sample_id),data=param_table,pch=21,size=3)+
  facet_wrap(~sample_id,nrow=2)+
  theme_bw(base_size=8)+
  theme(legend.position="none")+
  scale_x_continuous(breaks=seq(15,40,5))+
  xlab("Temperature (°C)")+ylab("Respiration Rate")

results<-stepwise_results%>%left_join(.,param_table,by="id")
summary(lm(topt_stepwise~topt,data=results))
b<-ggplot(results)+
  geom_abline()+
  geom_point(aes(topt_stepwise,topt))+
  geom_text_repel(aes(topt_stepwise,topt,label=id),size=3)+
  geom_smooth(aes(topt_stepwise,topt),method="lm",se=FALSE)+
  theme_bw(base_size=8)+
  coord_fixed()+
  scale_x_continuous(breaks=seq(20,35,1))+
  scale_y_continuous(breaks=seq(20,35,1))+
  ylab(expression(T[OPT]~"(Continuous)"))+
  xlab(expression(T[OPT]~"(Stepwise)"))

c<-ggplot(results)+
  geom_abline()+
  geom_point(aes(topt_stepwise,adj))+
  geom_text_repel(aes(topt_stepwise,adj,label=id),size=3)+
  geom_smooth(aes(topt_stepwise,adj),method="lm",se=FALSE)+
  theme_bw(base_size=8)+
  coord_fixed()+
  scale_x_continuous(breaks=seq(20,35,1))+
  scale_y_continuous(breaks=seq(20,35,1))+
  ylab(expression(T[OPT]~"(Adjusted, Continuous)"))+
  xlab(expression(T[OPT]~"(Stepwise)"))

quartz(w=7.2,h=3)
cowplot::plot_grid(a,b,c,nrow=1,labels="AUTO",label_size=8)
ggsave("~/Desktop/comparision_fig.pdf")


##### SHARPE ######

tref_val <- 27   # reference temp; pick what's sensible for your system
params <- setdiff(names(formals(get(mod, asNamespace("rTPC")))), c("temp", "tref"))
formula <- as.formula(
  sprintf("rate ~ %s(temp = temp, %s, tref = %s)",
          mod,
          paste(sprintf("%s = %s", params, params), collapse = ", "),
          tref_val)
)

fits <- working_tpc %>%
  group_by(sample_id) %>%
  group_nest() %>%
  mutate(
    start = map(data, ~ get_start_vals(.x$temp, .x$rate, model_name = mod)),
    low   = map(data, ~ get_lower_lims(.x$temp, .x$rate, model_name = mod)),
    high  = map(data, ~ get_upper_lims(.x$temp, .x$rate, model_name = mod)),
    fit = pmap(list(data, start, low, high),
               ~ nls_multstart(formula,
                               data = ..1,
                               iter = 500,
                               start_lower = ..2 - 10,
                               start_upper = ..2 + 10,
                               lower = ..3,
                               upper = ..4,
                               supp_errors = "Y",
                               lhstype = "random")),
    params_all = map(fit, ~ calc_params(.x) %>% as_tibble()),
    topt_obs = map_dbl(params_all, "topt"))

param_table <- fits %>%
  select(sample_id, params_all)%>%
  unnest(params_all)%>%
  mutate(adj=topt+1.5)%>%
  mutate(id=parse_number(sample_id));param_table

predicted <-fits %>%mutate(newdata = map(data, ~ tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 500))),
                           aug = map2(fit, newdata, ~ augment(.x, newdata = .y))) %>%
  select(sample_id, aug) %>%
  unnest(aug)%>%rename(rate=3)

a<-ggplot(working_tpc)+
  geom_hline(yintercept=0,color="red")+
  geom_point(aes(temp,rate),alpha=0.005)+
  geom_line(aes(temp,rate),data=predicted)+
  geom_point(aes(topt,rmax,fill=sample_id),data=param_table,pch=21,size=3)+
  facet_wrap(~sample_id,nrow=2)+
  theme_bw(base_size=8)+
  theme(legend.position="none")+
  scale_x_continuous(breaks=seq(15,40,5))+
  xlab("Temperature (°C)")+ylab("Respiration Rate")

results<-stepwise_results%>%left_join(.,param_table,by="id")
summary(lm(topt_stepwise~topt,data=results))
b<-ggplot(results)+
  geom_abline()+
  geom_point(aes(topt_stepwise,topt))+
  geom_text_repel(aes(topt_stepwise,topt,label=id),size=3)+
  geom_smooth(aes(topt_stepwise,topt),method="lm",se=FALSE)+
  theme_bw(base_size=8)+
  coord_fixed()+
  scale_x_continuous(breaks=seq(20,35,1))+
  scale_y_continuous(breaks=seq(20,35,1))+
  ylab(expression(T[OPT]~"(Continuous)"))+
  xlab(expression(T[OPT]~"(Stepwise)"))

c<-ggplot(results)+
  geom_abline()+
  geom_point(aes(topt_stepwise,adj))+
  geom_text_repel(aes(topt_stepwise,adj,label=id),size=3)+
  geom_smooth(aes(topt_stepwise,adj),method="lm",se=FALSE)+
  theme_bw(base_size=8)+
  coord_fixed()+
  scale_x_continuous(breaks=seq(20,35,1))+
  scale_y_continuous(breaks=seq(20,35,1))+
  ylab(expression(T[OPT]~"(Adjusted, Continuous)"))+
  xlab(expression(T[OPT]~"(Stepwise)"))

#quartz(w=7.2,h=3)
cowplot::plot_grid(a,b,c,nrow=1,labels="AUTO",label_size=8)
ggsave("~/Desktop/comparision_fig.pdf")

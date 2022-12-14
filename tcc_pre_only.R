## LSF_DOCKER_VOLUMES='/storage1/fs1/bafritz/Active/EICU\ Hot\ List/:/research/ /storage1/fs1/christopherking/Active/tcc/:/export/ /home/christopherking/gitdir/tcc_analysis/:/script_home/ ' bsub -G 'compute-christopherking' -n 7 -R 'rusage[mem=32GB] span[hosts=1]' -M 32GB -q general -a 'docker(cryanking/rstan:1.1)' R -f /script_home/tcc_pre_only.R
## LSF_DOCKER_VOLUMES='/storage1/fs1/bafritz/Active/EICU\ Hot\ List/:/research/ /storage1/fs1/christopherking/Active/tcc/:/export/ /home/christopherking/gitdir/tcc_analysis/:/script_home/' bsub -G 'compute-christopherking' -n 6 -R 'rusage[mem=32GB] span[hosts=1]' -M 32GB -q general-interactive -Is -a 'docker(cryanking/rstan:1.1)' /bin/bash

setwd("/research")
library(magrittr)
library(dplyr)
library(cmdstanr)
library(posterior)
set_cmdstan_path("/root/.cmdstan/cmdstan-2.30.1")
rawdata <- read.csv('hot list full data set with organ systems.csv',strip.white=T)
################################################
#Filter to the correct cohort
################################################
pts <- rawdata[!rawdata$actHospMort=='NULL',]

pts <- pts[pts$wardName %in% c('CCU','CTICU','SICU'),]

#For this analysis, I believe you will also want to filter by DCyr and DCqtr:
# pts <- pts[(pts$DCyr %in% c(2016,2017)) & (pts$DCqtr %in% c(1,2)),]

grouped_data <- pts %>% as_tibble %>% filter( (DCyr  > 2016) | (HotListPeriod1_6v6months != "NULL")  ) %>% filter(DCyr < 2018 )  %>% mutate(intergroup = HotListPeriod1_6v6months != "Pre")  %>% group_by(wardName, intergroup, DiseaseGrp) %>% summarize( apachel = qlogis(mean(as.numeric(predictedHospitalMortality))) , N=n(), y=sum(as.numeric(actHospMort), na.rm=T) ) %>% ungroup

grouped_data <- grouped_data %>%mutate(intergroup=as.numeric(intergroup)) %>% tidyr::pivot_wider(names_from='intergroup', values_from=one_of('apachel', 'N', 'y'),names_sep="" ) %>% filter(!is.na(N0) & !is.na(N1) )

setwd("/export")




## a sample transformation function

stan_code <- '
data {
  int<lower=0> J;          // number of cluster 
  array[J] int<lower=0> y0;      // outcomes in initial time
  array[J] int<lower=0> N0;      // at risk in pre time
  vector[J] apachel0;         // log of average apache for that group
}

parameters {
  // real drift; //parameters for the change in quality logit noise
  vector[J] eta0; // time0 quality logits
  real<lower=0> qual_sd; //parameters for the initial quality logit noise
}

model {
  eta0 ~ cauchy(0, qual_sd );
  qual_sd ~ cauchy(0,2);
  y0 ~ binomial_logit(N0, eta0 + apachel0);
}

' ## this set of modifications (high temperature, moderate drift_sd) allows rapid but very inefficent sampling
  ## it is basically a hierarchical model where each cluster-specific drift compensates for the difference betweeen q1 and q0, but it is "cheaper" for the prior to move mu than each drift [the cauchy prior imposes relatively little shrinkage]
  ## allowing individual group shocks to have high sd -> very low precision on mu (transformation function) [sd drift and sd mu are tightly co-identified]
  ## low temp -> slowish sampling and high aliasing of mu (every point is affected by every mu)

writeLines(stan_code, "simu_pre.stan")


mod <- cmdstan_model("simu_pre.stan")


fit <- mod$sample(
  data =  c( grouped_data%>% select(apachel0, y0, N0) %>% as.list , list(J=nrow(grouped_data)  ) ), 
  seed = 123, 
  chains = 12, 
  parallel_chains = 6,
  iter_warmup = 5000,
  iter_sampling = 10000,
  refresh = 1000 # print update every 500 iters
)

fit$save_object("tcc_pre_only.obj")

my_fit_summary <- fit$summary()

saveRDS(my_fit_summary, "tcc_pre_summary.rda")


# tcc_samples %>% select(starts_with("eta0")) %>% summarize_all(mean) %>% unlist %>% sort %>% data.frame
# fit <- readRDS("tcc_pre_only.obj")
# my_fit_summary <- readRDS("tcc_pre_summary.rda")

tcc_samples<-fit$draws(format="df")
tcc_samples <- tcc_samples[seq(from=1, to=nrow(tcc_samples), by=1000) ,]



new_N0 <- grouped_data %>% mutate(id=row_number()) %>% select(apachel1, y1, N1, id) %>%  mutate(expected_out = NA_real_, sd_out = NA_real_)

for( i in new_N0$id) {
sim_new <- rbinom(n=nrow(tcc_samples), size= rep(as.integer(new_N0[i,"y1"]), times=nrow(tcc_samples) ),  prob= plogis( as.numeric(new_N0[i, "apachel1"]) + as.numeric(tcc_samples[[paste0("eta0[", i,"]")]] ) ) )

new_N0[i, "expected_out"] <- mean(sim_new)
new_N0[i, "sd_out"] <- sd(sim_new)
}

fitted_ <- (tcc_samples %>% select(starts_with("eta0")) %>% as.matrix ) 



## match N1 and apache1 -> logit simulation







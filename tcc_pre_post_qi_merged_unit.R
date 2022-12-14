## LSF_DOCKER_VOLUMES='/storage1/fs1/bafritz/Active/EICU\ Hot\ List/:/research/ /storage1/fs1/christopherking/Active/tcc/:/export/ /home/christopherking/misc_applied_in_dev/tcc/:/script_home/ ' bsub -G 'compute-christopherking' -n 6 -R 'rusage[mem=32GB] span[hosts=1]' -M 32GB -q general -a 'docker(cryanking/rstan:1.1)' R -f /script_home/tcc_pre_post_qi.R
## LSF_DOCKER_VOLUMES='/storage1/fs1/bafritz/Active/EICU\ Hot\ List/:/research/ /storage1/fs1/christopherking/Active/tcc/:/export/ /home/christopherking/misc_applied_in_dev/tcc/:/script_home/' bsub -G 'compute-christopherking' -n 6 -R 'rusage[mem=32GB] span[hosts=1]' -M 32GB -q general-interactive -Is -a 'docker(cryanking/rstan:1.1)' /bin/bash

setwd("/research")
library(magrittr)
library(dplyr)
library(cmdstanr)
set_cmdstan_path("/root/.cmdstan/cmdstan-2.30.1")
library(posterior)

rawdata <- read.csv('hot list full data set with organ systems.csv',strip.white=T)
################################################
#Filter to the correct cohort
################################################
pts <- rawdata[!rawdata$actHospMort=='NULL',]

pts <- pts[pts$wardName %in% c('CCU','CTICU','SICU'),]

#For this analysis, I believe you will also want to filter by DCyr and DCqtr:
# pts <- pts[(pts$DCyr %in% c(2016,2017)) & (pts$DCqtr %in% c(1,2)),]

grouped_data <- pts %>% as_tibble %>% filter( (DCyr  > 2016) | (HotListPeriod1_6v6months != "NULL")  ) %>% filter(DCyr < 2018 )  %>% mutate(intergroup = HotListPeriod1_6v6months != "Pre")  %>% group_by(intergroup, DiseaseGrp) %>% summarize( apachel = qlogis(mean(as.numeric(predictedHospitalMortality))) , N=n(), y=sum(as.numeric(actHospMort), na.rm=T) ) %>% ungroup

grouped_data <- grouped_data %>%mutate(intergroup=as.numeric(intergroup)) %>% tidyr::pivot_wider(names_from='intergroup', values_from=one_of('apachel', 'N', 'y'),names_sep="" ) %>% filter(!is.na(N0) & !is.na(N1) )

setwd("/export")


## a sample transformation function

stan_code <- '
data {
  int<lower=0> J;          // number of cluster 
  int<lower=0> K;          // number of kernel points
  vector[K] kernPoints;      //location of kernel regression points
  array[J] int<lower=0> y0;      // outcomes in initial time
  array[J] int<lower=0> N0;      // at risk in pre time
  array[J] int<lower=0> y1;      // outcomes in post time
  array[J] int<lower=0> N1;      // at risk in post time
  vector[J] apachel0;         // log of average apache for that group
  vector[J] apachel1;         // log of average apache for that group
  // real<lower=0> drift_sd; //parameters for the change in quality logit noise
  real<lower=0> mu_sd; //parameters for the change in quality logit noise
  real<lower=0> temper; // at the short absolute distances, unscaled kernel regression smooths out over all points
}

parameters {
  // real drift; //parameters for the change in quality logit noise
  real<lower=0> drift_sd; //parameters for the change in quality logit noise
  real<lower=0> qual_sd; //parameters for the initial quality logit noise
  vector[J] drifts; // quality logits noise
  vector[J] eta0; // time0 quality logits
  vector[K] mu; // kernel regression values
}
transformed parameters {
  vector[J] eta1; //time1 quality logits
  for (n in 1:J) {
    //eta1[n] = sum(eta0[n]-kernPoints ); //todo add simplex intermediate?
     eta1[n] = sum(softmax( (-temper)* (eta0[n] - kernPoints) .* (eta0[n] - kernPoints)  )   .* mu); // distance of eta from each knot -> softmax standardized
  }
}
model {
  drifts ~ normal(0, drift_sd); // fixed vs prior on sd
  // drifts ~ normal(drift, drift_sd);
  drift_sd ~ cauchy(0,.2); // allowing high variance -> no penalty on shocks -> aliases away effect
  // drift ~ normal(0, .25);
  
  mu ~ cauchy(0,mu_sd); // vs normal with mu_sd ~ cauchy
  
  eta0 ~ cauchy(0, qual_sd );
  qual_sd ~ cauchy(0,2);
  y0 ~ binomial_logit(N0, eta0 + apachel0);
  y1 ~ binomial_logit(N1, eta0 + apachel1 + eta1 + drifts);
}

' ## this set of modifications (high temperature, moderate drift_sd) allows rapid but very inefficent sampling
  ## it is basically a hierarchical model where each cluster-specific drift compensates for the difference betweeen q1 and q0, but it is "cheaper" for the prior to move mu than each drift [the cauchy prior imposes relatively little shrinkage]
  ## allowing individual group shocks to have high sd -> very low precision on mu (transformation function) [sd drift and sd mu are tightly co-identified]
  ## low temp -> slowish sampling and high aliasing of mu (every point is affected by every mu)

writeLines(stan_code, "simu_drift.stan")


mod <- cmdstan_model("simu_drift.stan")

temp_used <- 50.

fit <- mod$sample(
  data =  c( grouped_data%>% select(apachel0, apachel1, y0, y1, N0, N1) %>% as.list , list(J=nrow(grouped_data) , drift_sd = .1, mu_sd = 2., K=5, kernPoints=seq(from=-.25, to=.25, length.out=5), temper=temp_used ) ), 
  seed = 123, 
  chains = 12, 
  parallel_chains = 6,
  iter_warmup = 15000,
  iter_sampling = 65000,
  refresh = 500 # print update every 500 iters
)

fit$save_object("tcc_pre_post_merge2.obj")



my_fit_summary <- fit$summary()

my_fit_summary %>% filter(grepl(variable, pattern="mu[", fixed=T))

## create the fitted curve

test_points <- seq(from=-.25, to=.25, length.out=100 )
kernPoints <- seq(from=-.25, to=.25, length.out=5) %>% rep(each=length(test_points)) %>% matrix(nrow=length(test_points))  %>% subtract(test_points ) %>% raise_to_power(2) %>% multiply_by(-temp_used) %>% exp %>% sweep(x=., MARGIN=1, STATS=rowSums(.), FUN="/" )

kernPoints %>% round(1) %>%set_colnames(paste0('x', seq(from=-.25, to=.25, length.out=5)) ) %>% cbind(test_points)

# fit <- readRDS("tcc_pre_post_first2.obj")

tcc_samples<-fit$draws(format="df")
fitted_mu <- (tcc_samples %>% select(starts_with("mu[")) %>% as.matrix ) %*% t(kernPoints )

fitted_mu %>% apply(2, mean) %>% exp -> mu_mean
fitted_mu %>% apply(2, quantile, prob=0.05) %>% exp -> mu_lower
fitted_mu %>% apply(2, quantile, prob=0.95) %>% exp-> mu_upper

pdf(file="eta_trend_qi_m.pdf")

# plot(cbind(test_points%>% exp  ,  mu_mean), ylim=c(min(mu_lower), max(mu_upper) )  , xlab="initial surplus mort ratio", ylab="change in surplus mort (ratio)", log="xy")
plot(cbind(test_points%>% exp  ,  mu_mean), ylim=c(1/3, 3 )  , xlab="initial surplus mort ratio", ylab="change in surplus mort (ratio)", log="xy")
points(x=test_points%>% exp , y=mu_lower, col="red", pch=19 )
points(x=test_points%>% exp , y=mu_upper, col="red", pch=19 )

abline(h=1)
abline(v=1)

dev.off()

## other plots: SMR pre-post with CI box, eta0 vs eta1

tcc_samples %>% select(starts_with("eta0")) %>% summarize_all(mean) %>% unlist %>% sort %>% data.frame

## qual_sd is extimated about 0.3 and almost always less than 0.5. -> almost all eta0 are <.6 and most are < 0.3 (exp = 1.35)
##  -> no real info on "far out" SMR, and most high SMR are illusions

## TODO: create a metric on the fitted curve

## estimate the SMR with a ci

pdf(file="eta_change_qi_m.pdf")
my_fit_summary %>% filter(grepl(variable, pattern="eta", fixed=T)) %>% select(variable, mean, sd) %>% mutate(period=grepl(variable, pattern="eta1")) %>% mutate(variable=substring(variable, first=5)) %>% tidyr::pivot_wider(names_from='period', values_from=c('mean', 'sd' ),names_sep="_" ) %>% select(mean_FALSE, mean_TRUE) %>% mutate_all(exp) %>% plot(xlab="initial quality factor", ylab="tcc multiplier", xlim=c(0.5,2), log="xy")
abline(h=1)
abline(v=1)
dev.off()

pdf(file="eta_density_qi_m.pdf")
my_fit_summary %>% filter(grepl(variable, pattern="eta0", fixed=T)) %>% mutate( mean = case_when(mean< -.6~-.65, mean>.6~.65, TRUE~mean) ) %>% select(mean) %>% unlist %>% hist(main="histogram of quality factor" ,axes=F, breaks = c( -.7, seq(from=-.6, to=.6, length.out=10 ), .7 ) )
axis(1, at= seq(from=-.5, to=.5, length.out=7), labels=seq(from=-.5, to=.5, length.out=7) %>% exp %>% round(1) )
dev.off()

pdf(file="eta_pre_post_qi_m.pdf")
my_fit_summary %>% 
  filter(grepl(variable, pattern="eta", fixed=T) | grepl(variable, pattern="drifts", fixed=T)) %>% 
  select(variable, mean,sd) %>% 
  mutate(period=case_when(
    grepl(variable, pattern="eta0") ~ 0 ,
    grepl(variable, pattern="eta1") ~ 1 ,
    grepl(variable, pattern="drift") ~ 2 )) %>% 
  mutate(group=sub(variable, pattern=".*\\[", replacement="") %>% sub(pattern=']', fixed=T, replacement='') ) %>% 
  tidyr::pivot_wider(names_from='period', values_from=c('mean', 'sd' ),names_sep="_" ,id_cols="group", names_prefix="v")  %>% 
  mutate(post=mean_v1+mean_v2+mean_v0, sd_post=sqrt(sd_v1^2+sd_v2^2) ) -> templ
  
  templ %>% select(mean_v0, post) %>% 
  mutate_all(exp) %>% plot(xlab="initial quality factor", ylab="post quality factor", xlim=c(0.5,2), log="xy", ylim=c(.4,1.5))
#   with(templ, arrows(x0=mean_v0 %>% exp, y0=exp(post+1.96*sd_post), y1=exp(post-1.96*sd_post) , length=0) )
abline(h=1)
abline(v=1)
abline(a=0, b=1)
dev.off()



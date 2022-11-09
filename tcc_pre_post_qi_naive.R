## LSF_DOCKER_VOLUMES='/storage1/fs1/bafritz/Active/EICU\ Hot\ List/:/research/ /storage1/fs1/christopherking/Active/tcc/:/export/ /home/christopherking/misc_applied_in_dev/tcc/:/script_home/' bsub -G 'compute-christopherking' -n 6 -R 'rusage[mem=32GB] span[hosts=1]' -M 32GB -q general-interactive -Is -a 'docker(cryanking/rstan:1.1)' /bin/bash

setwd("/research")
library(magrittr)
library(dplyr)
library(cmdstanr)
library(posterior)

rawdata <- read.csv('hot list full data set with organ systems.csv',strip.white=T)
################################################
#Filter to the correct cohort
################################################
pts <- rawdata[!rawdata$actHospMort=='NULL',]

pts <- pts[pts$wardName %in% c('CCU','CTICU','SICU'),]


raw_glm <- pts %>% as_tibble %>% filter( (DCyr  > 2016) | (HotListPeriod1_6v6months != "NULL")  ) %>% filter(DCyr < 2018 )  %>% mutate(postTCC = as.numeric(HotListPeriod1_6v6months != "Pre")) %>% mutate(y=as.numeric(actHospMort), apachel=qlogis(as.numeric(predictedHospitalMortality)  ) ) %>%
glm( y~ offset(apachel) + postTCC , data=. , family=binomial())

raw_glm %>% summary


#For this analysis, I believe you will also want to filter by DCyr and DCqtr:
# pts <- pts[(pts$DCyr %in% c(2016,2017)) & (pts$DCqtr %in% c(1,2)),]

grouped_data <- pts %>% as_tibble %>% filter( (DCyr  > 2016) | (HotListPeriod1_6v6months != "NULL")  ) %>% filter(DCyr < 2018 )  %>% mutate(intergroup = HotListPeriod1_6v6months != "Pre")  %>% group_by(wardName, intergroup, DiseaseGrp) %>% summarize( apachel = qlogis(mean(as.numeric(predictedHospitalMortality))) , N=n(), y=sum(as.numeric(actHospMort), na.rm=T) ) %>% ungroup

grouped_data <- grouped_data %>%mutate(intergroup=as.numeric(intergroup)) %>% tidyr::pivot_wider(names_from='intergroup', values_from=one_of('apachel', 'N', 'y'),names_sep="" ) %>% filter(!is.na(N0) & !is.na(N1) )

setwd("/export")

## is there evidence that there are any poorly performing groups?
## fixed effects doesn't really suggest so
smr_glm <- grouped_data %>% filter(y0>0) %>% filter(y0<N0) %>% filter(y0<N0) %>% mutate(grp=interaction(wardName, DiseaseGrp) ) %>% glm(cbind(y0,N0-y0) ~ offset(apachel0) + grp , data=. , family=binomial())



smr_glm %>% summary %>% extract2("coefficients") %>% (function(x) { extract(x, -1, -3) + cbind(rep(extract(x,1,1),nrow(x)-1) , rep(0, nrow(x)-1) ,rep(0,nrow(x)-1 ) ) }) %>%set_colnames(c("lOR","se","p" ) ) %>% as_tibble(rownames="grp") %>% mutate(OR = exp(lOR),lOR = round(lOR, 2), se=round(se, 2), p= round(p,4) ) %>% arrange(p) 


## these are sorted p-values
## there are maybe 1-2 likely poor-performing diagnoses

if(FALSE) {
## what about an lme4 fit (before pulling in any big guns)
library(lme4)
gm0 <- glmer(cbind(y0,N0-y0) ~ offset(apachel0) + (1|grp),
           family = binomial, data = grouped_data  %>% mutate(grp=interaction(wardName, DiseaseGrp) ))

gm00 <- glm(cbind(y0,N0-y0) ~ offset(apachel0) ,
           family = binomial, data = grouped_data %>% mutate(grp=interaction(wardName, DiseaseGrp) ))
                     
library(varTestnlme)
varCompTest(gm0,gm00)         
# p-value from exact weights: 0.01402405
## the glmer fit is convinced that there are non-null quality factors

gm1 <- glmer(cbind(y1,N1-y1) ~ offset(apachel1) + (1|grp),
           family = binomial, data = grouped_data  %>% mutate(grp=interaction(wardName, DiseaseGrp) ))


## shrinkage point estimates of the quality factor
gm0 %>% coef %>% extract2("grp") %>% as_tibble(rownames="rname") %>% arrange(-`(Intercept)`) %>% rename(qual_fact=`(Intercept)`) %>% mutate(qual_fact_e=exp(qual_fact) ) %>% head



gm0 %>% coef %>% extract2("grp") %>% as_tibble(rownames="rname") %>% arrange(`(Intercept)`) %>% rename(qual_fact=`(Intercept)`) %>% mutate(qual_fact_e=exp(qual_fact) ) %>% head


compare_blups <- inner_join( 
  gm0 %>% coef %>% extract2("grp") %>% as_tibble(rownames="rname") %>% arrange(-`(Intercept)`) %>% rename(qual_fact0=`(Intercept)`) %>% mutate(qual_fact_e0=exp(qual_fact0) ) , 
  gm1 %>% coef %>% extract2("grp") %>% as_tibble(rownames="rname") %>% arrange(-`(Intercept)`) %>% rename(qual_fact1=`(Intercept)`) %>% mutate(qual_fact_e1=exp(qual_fact1) ) , by="rname" )

pdf(file="blup_change.pdf")
compare_blups %>% select(qual_fact_e0, qual_fact_e1) %>% plot( xlab="initial quality factor", ylab="post quality factor")
abline(h=1)
abline(v=1)
abline(a=0, b=1)
dev.off()


}


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
  vector[J] etaprime; // time0 quality logits naive
  vector[K] mu; // kernel regression values
}
transformed parameters {
  vector[J] eta1; //time1 quality logits
  for (n in 1:J) {
    //eta1[n] = sum(etaprime[n]-kernPoints ); //todo add simplex intermediate?
     eta1[n] = sum(softmax( (-temper)* (etaprime[n] - kernPoints) .* (etaprime[n] - kernPoints)  )   .* mu); // distance of eta from each knot -> softmax standardized
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
  etaprime ~ cauchy(0, 5 );
  y0 ~ binomial_logit(N0, eta0 + apachel0);
  y0 ~ binomial_logit(N0, etaprime + apachel0);
  y1 ~ binomial_logit(N1, eta0 + apachel1 + eta1 + drifts);
}

' ## this set of modifications (high temperature, moderate drift_sd) allows rapid but very inefficent sampling
  ## it is basically a hierarchical model where each cluster-specific drift compensates for the difference betweeen q1 and q0, but it is "cheaper" for the prior to move mu than each drift [the cauchy prior imposes relatively little shrinkage]
  ## allowing individual group shocks to have high sd -> very low precision on mu (transformation function) [sd drift and sd mu are tightly co-identified]
  ## low temp -> slowish sampling and high aliasing of mu (every point is affected by every mu)

writeLines(stan_code, "simu_drift2.stan")


mod <- cmdstan_model("simu_drift2.stan")

temp_used <- 70.

fit2 <- mod$sample(
  data =  c( grouped_data%>% select(apachel0, apachel1, y0, y1, N0, N1) %>% as.list , list(J=nrow(grouped_data) , drift_sd = .1, mu_sd = 2., K=5, kernPoints=seq(from=-.75, to=.75, length.out=5), temper=temp_used ) ), 
  seed = 123, 
  chains = 12, 
  parallel_chains = 6,
  iter_warmup = 15000,
  iter_sampling = 65000,
  refresh = 500 # print update every 500 iters
)

fit2$save_object("tcc_pre_post_first_naive.obj")

# fit2 <- readRDS("tcc_pre_post_first_naive.obj")


my_fit2_summary <- fit2$summary()

my_fit2_summary %>% filter(grepl(variable, pattern="mu[", fixed=T))

## create the fit2ted curve

test_points <- seq(from=-.75, to=.75, length.out=100 )
kernPoints <- seq(from=-.75, to=.75, length.out=5) %>% rep(each=length(test_points)) %>% matrix(nrow=length(test_points))  %>% subtract(test_points ) %>% raise_to_power(2) %>% multiply_by(-temp_used) %>% exp %>% sweep(x=., MARGIN=1, STATS=rowSums(.), FUN="/" )

kernPoints %>% round(1) %>%set_colnames(paste0('x', seq(from=-.25, to=.25, length.out=5)) ) %>% cbind(test_points)


tcc_samples<-fit2$draws(format="df")
fitted_mu <- (tcc_samples %>% select(starts_with("mu[")) %>% as.matrix ) %*% t(kernPoints )

fitted_mu %>% apply(2, mean) %>% exp -> mu_mean
fitted_mu %>% apply(2, quantile, prob=0.05) %>% exp -> mu_lower
fitted_mu %>% apply(2, quantile, prob=0.95) %>% exp-> mu_upper

pdf(file="eta_trend_naive.pdf")

# plot(cbind(test_points%>% exp  ,  mu_mean), ylim=c(min(mu_lower), max(mu_upper) )  , xlab="initial surplus mort ratio", ylab="change in surplus mort (ratio)", log="xy")
plot(cbind(test_points%>% exp  ,  mu_mean), ylim=c(1/3, 3 )  , xlab="initial surplus mort ratio", ylab="change in surplus mort (ratio)", log="xy")
points(x=test_points%>% exp , y=mu_lower, col="red", pch=19 )
points(x=test_points%>% exp , y=mu_upper, col="red", pch=19 )

abline(h=1)
abline(v=1)

dev.off()

## other plots: SMR pre-post with CI box, eta0 vs eta1

tcc_samples %>% select(starts_with("eta0")) %>% summarize_all(mean) %>% unlist %>% sort %>% data.frame

## create a metric on the fitted curve

## estimate the SMR with a ci

pdf(file="eta_change_naive.pdf")
my_fit2_summary %>% filter(grepl(variable, pattern="eta\\d", fixed=F)) %>% select(variable, mean, sd) %>% mutate(period=grepl(variable, pattern="eta1")) %>% mutate(variable=substring(variable, first=5)) %>% tidyr::pivot_wider(names_from='period', values_from=c('mean', 'sd' ),names_sep="_" ) %>% select(mean_FALSE, mean_TRUE) %>% mutate_all(exp) %>% plot(xlab="initial quality factor", ylab="tcc multiplier", xlim=c(0.5,2), log="xy")
abline(h=1)
abline(v=1)
dev.off()

pdf(file="eta_prime_change_naive.pdf")
my_fit2_summary %>% filter(grepl(variable, pattern="eta1", fixed=T) | grepl(variable, pattern="etaprime", fixed=T)) %>% select(variable, mean, sd) %>% mutate(period=grepl(variable, pattern="eta1")) %>% mutate(variable=sub(variable, pattern=".*\\[", replacement="")) %>% tidyr::pivot_wider(names_from='period', values_from=c('mean', 'sd' ),names_sep="_" ) %>% select(mean_FALSE, mean_TRUE) %>% mutate_all(exp) %>% plot(xlab="initial SMR", ylab="tcc multiplier", xlim=c(0.5,3.5), log="xy")
abline(h=1)
abline(v=1)
dev.off()


pdf(file="eta_density_naive.pdf")
my_fit2_summary %>% filter(grepl(variable, pattern="eta0", fixed=T)) %>% mutate( mean = case_when(mean< -.6~-.65, mean>.6~.65, TRUE~mean) ) %>% select(mean) %>% unlist %>% hist(main="histogram of quality factor" ,axes=F, breaks = c( -.7, seq(from=-.6, to=.6, length.out=10 ), .7 ) )
axis(1, at= seq(from=-.5, to=.5, length.out=7), labels=seq(from=-.5, to=.5, length.out=7) %>% exp %>% round(1) )
dev.off()

pdf(file="eta_prme_density_naive.pdf")
my_fit2_summary %>% filter(grepl(variable, pattern="etaprime", fixed=T)) %>% mutate( mean = case_when(mean< -.6~-.65, mean>.6~.65, TRUE~mean) ) %>% select(mean) %>% unlist %>% hist(main="histogram of quality factor" ,axes=F, breaks = c( -.7, seq(from=-.6, to=.6, length.out=10 ), .7 ) )
axis(1, at= seq(from=-.5, to=.5, length.out=7), labels=seq(from=-.5, to=.5, length.out=7) %>% exp %>% round(1) )
dev.off()



pdf(file="eta_pre_post_naive.pdf")
my_fit2_summary %>% 
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

my_fit2_summary %>% 
  filter(grepl(variable, pattern="eta", fixed=T) ) %>% 
  select(variable, mean) %>% 
  mutate(group=sub(variable, pattern=".*\\[", replacement="") %>% sub(pattern=']', fixed=T, replacement='') ) %>% 
  mutate(type = sub(variable, pattern="\\[.*", replacement="")) %>%
  tidyr::pivot_wider(names_from='type', values_from=c('mean' ),names_sep="_" ,id_cols="group") %>% 
  inner_join(grouped_data %>% mutate(group=as.character(row_number())), by="group"  ) %>%
  arrange(eta0) %>% select(-wardName, -group) %>% mutate(SMR0= (y0)/(N0)/plogis(apachel0) ) %>% select(-starts_with("apache"))   %>% mutate( across(starts_with("eta"), exp ) ) %>% select(-one_of(c('N1','y1', 'eta1')))  %>% print(n=50)
  

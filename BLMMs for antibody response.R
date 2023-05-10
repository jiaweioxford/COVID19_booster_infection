####################################################################################################
# Example of code written by Jia Wei (Nuffield Department of Medicine, University of Oxford)       #
# for analyses of antibody waning following booster vaccination or breakthrough infection.         #
# Accompanying paper: Protection against SARS-CoV-2 Omicron BA.4/5 variant following booster       #
# vaccination or breakthrough infection in the UK                                                  #
####################################################################################################

library(tidyverse)
library(brms)


###piecewise model examining antibody responses from 21 days post-second vaccination to 14 days post-third/booster vaccination or infection 

data1=data1 %>% 
  filter(time<=14) %>% 
  mutate(knot1=0) %>% 
  mutate(T1=(time<knot1)*time7+(time>=knot1)*knot1,
         T2=(time>=knot1)*(time-knot1)+(time<knot1)*knot1) %>% 
  group_by(id)%>%
  mutate(ID=cur_group_id()) %>% 
  ungroup() %>% 
  mutate(ID=as.factor(ID))

priors=c(
  set_prior("normal(6,2)",class="Intercept"),
  set_prior("normal(0,1)",coef="T1",class="b"),
  set_prior("normal(0,1)",coef="T2",class="b"),
  set_prior("normal(0,1)",coef="age10",class="b"),
  set_prior("normal(0,1)",coef="Male1",class="b"),
  set_prior("normal(0,1)",coef="ethnicity1",class="b"),
  set_prior("normal(0,1)",coef="lthc1",class="b"),
  set_prior("normal(0,1)",coef="hcw1",class="b"),
  set_prior("normal(0,1)",coef="dur30",class="b"),
  
  set_prior("normal(0,0.1)",coef="age10:T1",class="b"),
  set_prior("normal(0,0.1)",coef="T2:age1",class="b"),
  set_prior("normal(0,0.1)",coef="T2:age2",class="b"),
  
  set_prior("normal(0,0.1)",coef="T1:Male1",class="b"),
  set_prior("normal(0,0.1)",coef="T2:Male1",class="b"),
  
  set_prior("normal(0,0.1)",coef="T1:ethnicity1",class="b"),
  set_prior("normal(0,0.1)",coef="T2:ethnicity1",class="b"),
  
  set_prior("normal(0,0.1)",coef="T1:lthc1",class="b"),
  set_prior("normal(0,0.1)",coef="T2:lthc1",class="b"),
  
  set_prior("normal(0,0.1)",coef="T1:hcw1",class="b"),
  set_prior("normal(0,0.1)",coef="T2:hcw1",class="b"),
  
  set_prior("normal(0,0.1)",coef="T1:dur30",class="b"),
  set_prior("normal(0,0.1)",coef="dur30:T2",class="b"),
  
  set_prior("cauchy(0,1)",coef="T1",class="sd",group="ID"),
  set_prior("cauchy(0,1)",coef="T2",class="sd",group="ID"),
  set_prior("cauchy(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("cauchy(0,1)",class="sigma"),
  set_prior("lkj_corr_cholesky(1)",class="L")
)


set_inits=function(seed=1){
  set.seed(seed)
  list(
    Intercept=rnorm(1,6,2),
    b=c(rnorm(1,0,1),#T1
        rnorm(1,0,1),#T2
        rnorm(1,0,1),#age
        rnorm(1,0,1),#sex
        rnorm(1,0,1),#ethnicity
        rnorm(1,0,1),#lthc
        rnorm(1,0,1),#hcw
        rnorm(1,0,1),#dur
        
        rnorm(1,0,1),#age:T1
        rnorm(1,0,1),#age:T2
        rnorm(1,0,1),#age:T2
        rnorm(1,0,1),#sex:T1
        rnorm(1,0,1),#sex:T2
        rnorm(1,0,1),#ethnicity:T1
        rnorm(1,0,1),#ethnicity:T2
        rnorm(1,0,1),#lthc:T1
        rnorm(1,0,1),#lthc:T2
        rnorm(1,0,1),#hcw:T1
        rnorm(1,0,1),#hcw:T2
        rnorm(1,0,1),#dur:T1
        rnorm(1,0,1) #prior:T2
    ),
    sigma=runif(1,0,1),
    sd_1=c(runif(1,0,1),runif(1,0,1),runif(1,0,1)),
    z_1=matrix(rep(c(6,-0.1,1),21448),3,21448)
  )
}

inits_list=list(
  set_inits(1),
  set_inits(2),
  set_inits(3),
  set_inits(4),
  set_inits(5),
  set_inits(6)
)

fit1 <- brm(formula=log2(assay)|cens(cen1)~1+T1:age10+age10+T2:age1+T2:age2+
                       T1*dur30+T2*dur30+T1*Male+T2*Male+
                       T1*ethnicity+T2*ethnicity+
                       T1*lthc+T2*lthc+T1*hcw+T2*hcw+(1+T1+T2|ID),
                     data=data1,cores=6,family=gaussian(),
                     prior = priors,inits=inits_list,
                     chains = 6, iter=4000, warmup = 2000, seed=42,
                     control = list(adapt_delta=0.95))


## linear model examining antibody decline from 42 days post third/booster vaccination or infection

priors=c(
  set_prior("normal(10,2)",class="Intercept"),
  set_prior("normal(0,0.5)",coef="T3",class="b"),
  set_prior("normal(0,1)",coef="age1",class="b"),
  set_prior("normal(0,1)",coef="age2",class="b"),
  set_prior("normal(0,1)",coef="Male1",class="b"),
  set_prior("normal(0,1)",coef="ethnicity1",class="b"),
  set_prior("normal(0,1)",coef="lthc1",class="b"),
  set_prior("normal(0,1)",coef="hcw1",class="b"),
  set_prior("normal(0,1)",coef="dur30",class="b"),
  
  set_prior("normal(0,0.1)",coef="T3:age1",class="b"),
  set_prior("normal(0,0.1)",coef="T3:age2",class="b"),
  set_prior("normal(0,0.1)",coef="T3:Male1",class="b"),
  set_prior("normal(0,0.1)",coef="T3:ethnicity1",class="b"),
  set_prior("normal(0,0.1)",coef="T3:lthc1",class="b"),
  set_prior("normal(0,0.1)",coef="T3:hcw1",class="b"),
  set_prior("normal(0,0.1)",coef="T3:dur30",class="b"),
  
  set_prior("normal(0,0.1)",coef="T3",class="sd",group="ID"),
  set_prior("normal(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("normal(0,1)",class="sigma"),
  set_prior("lkj_corr_cholesky(1)",class="L")
)

set_inits=function(seed=1){
  set.seed(seed)
  list(
    Intercept=rnorm(1,10,2),
    b=c(rnorm(1,0,0.5),#time
        rnorm(1,0,1),#age1
        rnorm(1,0,1),#age2
        rnorm(1,0,1),#sex
        rnorm(1,0,1),#ethnicity
        rnorm(1,0,1),#lthc
        rnorm(1,0,1),#hcw
        rnorm(1,0,1),#dur
        
        rnorm(1,0,0.1),#age1
        rnorm(1,0,0.1),#age2
        rnorm(1,0,0.1),#sex
        rnorm(1,0,0.1),#ethnicity
        rnorm(1,0,0.1),#lthc
        rnorm(1,0,0.1),#hcw
        rnorm(1,0,0.1),#dur
    ),
    sigma=runif(1,0,1),
    sd_1=c(runif(1,0,1),runif(1,0,0.1)),
    z_1=matrix(rep(c(10,-0.1),14748),2,14748)
  )
}

inits_list=list(
  set_inits(1),
  set_inits(2),
  set_inits(3),
  set_inits(4),
  set_inits(5),
  set_inits(6)
)

fit2 <- brm(formula=log2(assay)|cens(cen1)~1+T3*age1+T3*age2+T3*Male+
              T3*ethnicity+T3*lthc+T3*hcw+T3*dur30+(1+T3|ID),
                          data=data2,cores=6,family=gaussian(),
                          prior = priors,inits=inits_list,
                          chains = 6, iter=4000, warmup = 2000, seed=42,
                          control = list(adapt_delta=0.95))

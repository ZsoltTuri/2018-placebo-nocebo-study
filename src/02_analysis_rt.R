# This script is part of 
#
# Turi, Z., Bjørkedal, E., Gunkel, L., Antal, A., Paulus, W. & Mittner, M. (2018).
# Evidence for Cognitive Placebo and Nocebo Effects in Healthy Individuals.
#
# Analysis of reaction time data.
#
#
library(ProjectTemplate)
load.project()

options(mc.cores=parallel::detectCores())
theme_set(theme_bw())

bname<-tools::file_path_sans_ext(basename(this.file.name()))
#stop()

#======================== 
## data preparation
#========================
learn %>%
  mutate(group = fct_relevel(group, "nhg"),
         symbl_pair = factor(symbl_pair, levels = c("AB", "CD", "EF")),
         symbl_pair = fct_relevel(symbl_pair, "AB"),
         day = fct_relevel(day, "1")) %>%
  group_by(day, group, participant, symbl_pair) %>%
  mutate(ztrial = scale(trial)) %>% 
  na.omit %>% ungroup  %>%
  group_by(participant) %>% mutate(zrt=(reaction_time-mean(reaction_time))/sd(reaction_time)) %>% ungroup -> d
str(d)


#========================
## function applied to every model for fitting and plotting
#========================

fit_and_plot <- function(mod.name,frm){
  #mod.name = formula.name.fname(frm)
  
  if(!is.cached.var(mod.name, base=bname)){
    mod <- brm(frm, data = d, family = gaussian())
    assign(mod.name, mod, envir=.GlobalEnv)
    cache.var(mod.name, bname)
  } else {
    mod <- load.cache.var(mod.name,bname)
  }
  
  pdf(plot.filename(sprintf("diag_%s.pdf", mod.name), base=bname), width=5,height=10)
  mcmc_rhat(rhat(mod)) %>% print
  mcmc_neff(neff_ratio(mod)) %>% print
  dev.off()
  
  if(dim(attr(terms(frm), "factors"))[2]>1){
    ## needs to catch error because null model does not have effects
    pdf(plot.filename(sprintf("marginal_%s.pdf", mod.name), base=bname), width=4,height=3)
    map(plot(marginal_effects(mod), ask=F,plot=F), function(p) p+ylim(0,1)) %>% print
    plot(marginal_effects(mod, probs = c(0.05, 0.95)), ask=F)  
    dev.off()
  } else {
    message("Not plotting marginal effects for ", frm);
  }
  
  pdf(plot.filename(sprintf("ppred_%s.pdf", mod.name), base=bname), width=15,height=10)
  d %>% add_predicted_samples(mod, n=100) %>%
    group_by(day, group, participant, symbl_pair, .iteration) %>%
    summarize(pred=mean(pred),reaction_time=mean(reaction_time)) %>%
    ggplot(aes(x=symbl_pair, y=reaction_time))+
    stat_pointinterval(aes(y = pred), .prob = c(.99, .95, .8, .5))+
    #geom_point(size=2, color="red")+
    stat_summary(fun.y=mean, size=2, color="red", geom="point")+
    facet_wrap(~participant) -> p
  print(p)
  
  d %>% add_predicted_samples(mod, n=100) %>%
    group_by(day, participant, group, symbl_pair, .iteration) %>%
    summarize(pred=mean(pred),reaction_time=mean(reaction_time)) %>%
    ggplot(aes(x=symbl_pair, y=reaction_time))+
    stat_pointinterval(aes(y = pred), .prob = c(.99, .95, .8, .5))+
    stat_summary(fun.y=mean, size=2, color="red", geom="point")+
    facet_grid(group~participant) + 
    labs(title="mean log-transformed rt per subject split by symbol pair/group") -> p
  print(p)  
  d %>% add_predicted_samples(mod, n=100) %>%
    group_by(day, participant, group, symbl_pair, .iteration) %>%
    summarize(pred=mean(pred),reaction_time=mean(reaction_time)) %>%
    ggplot(aes(x=symbl_pair, y=reaction_time))+
    stat_pointinterval(aes(y = pred), .prob = c(.99, .95, .8, .5))+
    stat_summary(fun.y=mean, size=2, color="red", geom="point")+
    facet_grid(day~group~participant) + 
    labs(title="mean log-transformed rt per subject split by symbol pair/group/day") -> p
  print(p) 
  dev.off()
  
  return(mod)
}

#========================
## model definitions
#========================
models <- list(
  formula(reaction_time ~ (1|subj)),
  formula(reaction_time ~ ztrial * symbl_pair + day + day:ztrial + (1 + ztrial + symbl_pair|subj)),
  formula(reaction_time ~ ztrial * symbl_pair + day + day:ztrial + group + (1 + ztrial + symbl_pair|subj)),
  formula(reaction_time ~ ztrial * symbl_pair + day + day:ztrial + group + group:day + (1 + ztrial + symbl_pair|subj))

)
names(models) <- sprintf("mod%02i", 1:length(models))

#========================
## fit models
#========================
library(parallel)
models.fitted=mcmapply(fit_and_plot, names(models), models, mc.cores = 4, SIMPLIFY = FALSE)
# stop()
#========================
## model-selection
#========================

loos = if.cached.load("loos", {
  map(models.fitted, LOO, pointwise=F) 
}, base=bname)
loos=map2(names(models), loos, function(mname,mod){mod$model_name=mname; mod})

sink(log.filename("models.log", base=bname))

tibble(
  model_num=names(models),
  model=paste(format(models), sep="")
) %>% print

compare_ic(x=loos) %>% print

sink(NULL)

mod04<-load.cache.var("mod04", bname)
mcmc_intervals(as.matrix(mod04), regex_pars = "b_.*")+xlim(-.1,.2)
hypothesis(mod04, c("day2:groupp_cond<day2:groupp_cntrl",
                    "day2:groupp_cond<day2:groupn_cond"))

m4=as.matrix(mod04)
colnames(m4)
eff = m4[,"b_day2:groupp_cond"]-m4[,"b_day2:groupp_cntrl"]
hist(eff)
sum(eff<0)/sum(eff>=0)
sum(eff<0)/length(eff)

library(tidyverse)
library(HDInterval)
d4=as.data.frame(mod04)
hdi.lower=function(x) hdi(x)[1]
hdi.upper=function(x) hdi(x)[2]

d4 %>% mutate(
  diff_pcond_pctrl=b_day2.groupp_cond-b_day2.groupp_cntrl,
  diff_pcond_nctrl=b_day2.groupp_cond-b_day2.groupn_cntrl,
  diff_pcond_ncond=b_day2.groupp_cond-b_day2.groupn_cond
) %>% select(starts_with("diff")) %>%
  summarize_all(funs(mean,hdi.lower,hdi.upper))
names(d4)

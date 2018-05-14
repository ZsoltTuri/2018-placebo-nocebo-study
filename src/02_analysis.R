library(ProjectTemplate)
load.project()

library(brms)
library(bayesplot)
library(tidybayes)
library(forcats)

options(mc.cores=parallel::detectCores())
theme_set(theme_bw())

bname<-tools::file_path_sans_ext(basename(this.file.name()))
# stop()

#======================== 
## data preparation
#========================
learn <- rbind(day1.learn, day2.learn)
learn %>%
  mutate(group = fct_relevel(group, "nhg"),
         symbl_pair = factor(symbl_pair, levels = c("AB", "CD", "EF")),
         symbl_pair = fct_relevel(symbl_pair, "AB"),
         day = fct_relevel(day, "1")) %>%
  group_by(day, group, participant, symbl_pair) %>%
  mutate(ztrial = scale(trial)) %>% ungroup -> d
str(d)

#========================
## function applied to every model for fitting and plotting
#========================

fit_and_plot <- function(mod.name,frm){
  #mod.name = formula.name.fname(frm)
  
  if(!is.cached.var(mod.name, base=bname)){
    mod <- brm(frm, data = d, family = bernoulli())
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
    summarize(pred=mean(pred),accuracy=mean(accuracy)) %>%
    ggplot(aes(x=symbl_pair, y=accuracy))+
    stat_pointinterval(aes(y = pred), .prob = c(.99, .95, .8, .5))+
    #geom_point(size=2, color="red")+
    stat_summary(fun.y=mean, size=2, color="red", geom="point")+
    facet_wrap(~participant) -> p
  print(p)
  
  d %>% add_predicted_samples(mod, n=100) %>%
    group_by(day, participant, group, symbl_pair, .iteration) %>%
    summarize(pred=mean(pred),accuracy=mean(accuracy)) %>%
    ggplot(aes(x=symbl_pair, y=accuracy))+
    stat_pointinterval(aes(y = pred), .prob = c(.99, .95, .8, .5))+
    stat_summary(fun.y=mean, size=2, color="red", geom="point")+
    facet_grid(group~participant) + 
    labs(title="mean acc per subject split by symbol pair/group") -> p
  print(p)  
  d %>% add_predicted_samples(mod, n=100) %>%
    group_by(day, participant, group, symbl_pair, .iteration) %>%
    summarize(pred=mean(pred),accuracy=mean(accuracy)) %>%
    ggplot(aes(x=symbl_pair, y=accuracy))+
    stat_pointinterval(aes(y = pred), .prob = c(.99, .95, .8, .5))+
    stat_summary(fun.y=mean, size=2, color="red", geom="point")+
    facet_grid(day~group~participant) + 
    labs(title="mean acc per subject split by symbol pair/group/day") -> p
  print(p) 
  dev.off()
  
  return(mod)
}

#========================
## model definitions
#========================
models <- list(
  formula(accuracy ~ (1|subj)),
  formula(accuracy ~ ztrial + (1|subj)),
  formula(accuracy ~ symbl_pair + (1|subj)),
  formula(accuracy ~ ztrial + symbl_pair + (1|subj)),
  formula(accuracy ~ ztrial * symbl_pair + (1|subj)),
  formula(accuracy ~ ztrial * symbl_pair + day + (1|subj)),
  formula(accuracy ~ ztrial * symbl_pair + day + day:ztrial + (1|subj)),
  formula(accuracy ~ ztrial * symbl_pair + day + day:ztrial + (1 + ztrial|subj)),
  formula(accuracy ~ ztrial * symbl_pair + day + day:ztrial + (1 + ztrial + symbl_pair|subj)),
  formula(accuracy ~ ztrial * symbl_pair + day + day:ztrial + group + (1 + ztrial + symbl_pair|subj)),
  formula(accuracy ~ ztrial * symbl_pair + day + day:ztrial + group + group:day + (1 + ztrial + symbl_pair|subj)),
  formula(accuracy ~ ztrial * symbl_pair + day + day:ztrial + group + group:day + day:ztrial:group + (1 + ztrial + symbl_pair|subj))
  
  # formula(accuracy ~ ztrial * symbl_pair + group + (1|participant)),
  # formula(accuracy ~ ztrial * symbl_pair * group + (1|participant)),
  # formula(accuracy ~ ztrial * symbl_pair * group + day + (1|participant)),
  # formula(accuracy ~ ztrial * symbl_pair * group * day + (1|participant))
  
  # formula(accuracy ~ cue_type * ztrial + (1+cue_type|participant)),
  # formula(accuracy ~ cue_type * ztrial + (1+ztrial|participant)),
  # formula(accuracy ~ cue_type * ztrial + (1+cue_type+ztrial|participant)),
  # formula(accuracy ~ cue_type * ztrial + session + (1+cue_type+ztrial|participant)),
  # formula(accuracy ~ cue_type * ztrial * session + (1+cue_type+ztrial|participant)),
  # formula(accuracy ~ cue_type * ztrial * session + (1+cue_type+ztrial+session|participant))
)
names(models) <- sprintf("mod%02i", 1:length(models))

#========================
## fit models
#========================
library(parallel)
#models.fitted=map(models, fit_and_plot)
#models.fitted=mclapply(models, fit_and_plot, mc.cores = 4)
#models.fitted=map2(names(models), models, fit_and_plot)
models.fitted=mcmapply(fit_and_plot, names(models), models, mc.cores = 4, SIMPLIFY = FALSE)
# stop()
#========================
## model-selection
#========================

loos = if.cached.load("loos", {
  #mclapply(models.fitted, LOO, mc.cores = 4)
  map(models.fitted, LOO, pointwise=F) 
}, base=bname)
##loos = mclapply(models.fitted, LOO, mc.cores=6)
loos=map2(names(models), loos, function(mname,mod){mod$model_name=mname; mod})

sink(log.filename("modsel.log", base=bname))

tibble(
  model_num=names(models),
  model=paste(format(models), sep="")
) %>% print
compare_ic(x=loos)

sink(NULL)

# m9=as.matrix(mod09)
# mcmc_intervals(m9, regex_pars = c("b_"))
#========================
## exploratory
#========================
# mod11
mcmc_intervals(as.matrix(mod12), regex_pars = c("b_"))

# based on model-selection does adding session improve the model?

#d %>% group_by(participant,session) %>% summarise(perc_nogo=sum(is.na(reaction_time))/n()*100) %>%
#  ggplot(aes(x=participant, y=perc_nogo))+geom_bar(stat="identity")+facet_wrap(~session)

# 
# d %>% 
#   ggplot(aes(x = day,y = accuracy))+
#   stat_summary(fun.y = mean, geom = "bar")+
#   facet_wrap(~day)
# 
# d %>% filter(participant==7, session=="stim_7Hz") %>% View
# group_by(cue_type) %>% 
#   summarize(mean(accuracy), n(), sum(accuracy==1), sum(accuracy==0), sum(is.na(reaction_time)))# min(accuracy),max(accuracy))
# 
# 
# d %>% group_by(participant,session) %>% summarise(perc_nogo=sum(is.na(reaction_time))/n()*100) %>%
#   ggplot(aes(x=participant, y=perc_nogo))+geom_bar(stat="identity")+facet_wrap(~session)

library(tidyverse)
d11=as.data.frame(mod11)
hdi.lower=function(x) hdi(x)[1]
hdi.upper=function(x) hdi(x)[2]

d11 %>% mutate(
  diff_pcond_pctrl=b_day2.groupp_cond-b_day2.groupp_cntrl,
  diff_pcond_nctrl=b_day2.groupp_cond-b_day2.groupn_cntrl,
  diff_pcond_ncond=b_day2.groupp_cond-b_day2.groupn_cond
) %>% select(starts_with("diff")) %>%
  summarize_all(funs(mean,hdi.lower,hdi.upper))
names(d11)
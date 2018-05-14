library(ProjectTemplate)
load.project()

library(brms)
library(bayesplot)
library(tidybayes)
library(forcats)

options(mc.cores=parallel::detectCores())
theme_set(theme_bw())

bname<-tools::file_path_sans_ext(basename(this.file.name()))
#stop()

str(subj.outcomes)
subj.outcomes %>%
  gather(time, value, expected, perceived) %>%
  mutate(group = fct_relevel(group, "n_cntrl"),
         time = factor(time),
         value = factor(value, levels=c("decline", "neutral", "improve")),
         value = fct_relevel(value, "neutral"),
         effect = fct_relevel(effect, "placebo"),
         type = fct_relevel(type, "control")) -> d
str(d)

# library(brms)
options(mc.cores=parallel::detectCores())

if(F){
mod1 <- brm(value ~ effect + (1|subj), data = d, family = "categorical", prior=c(set_prior ("cauchy (0, 1)")))
mod2 <- brm(value ~ effect + type + (1|subj), data = d, family = "categorical", prior=c(set_prior ("cauchy (0, 1)")))
mod3 <- brm(value ~ effect * type + (1|subj), data = d, family = "categorical", prior=c(set_prior ("cauchy (0, 1)")))
mod4 <- brm(value ~ effect * type + time + (1|subj), data = d, family = "categorical", prior=c(set_prior ("cauchy (0, 1)")))
mod5 <- brm(value ~ effect * type * time + (1|subj), data = d, family = "categorical", prior=c(set_prior ("cauchy (0, 1)")))
}

mod6 <- if.cached.load("mod6", 
                       brm(value ~ effect * type + time + time:type+ (1|subj), data = d, family = "categorical", prior=c(set_prior ("cauchy (0, 1)"))),
                       base=bname)



# plot.postpred <- function(mod, nrep = 100){
#   m = posterior_predict(mod0, nsamples = nrep) 
#   d.rep = cbind(d %>% select(-value), data.frame(m%>%t) %>% gather(value)) %>% 
#     mutate(value = value - 1) %>% 
#     group_by(groups, time) %>% summarize(count=n()) 
#   
#   d %>%
#     ggplot(aes(x = value)) + geom_bar(stat = "bin") + facet_grid(time ~ groups)+
#     stat_summary(fun.data = mean_sdl, aes(x = time, y = count),fun.args = list(mult = 1),
#                  geom = "pointrange", color = "grey",
#                  data = d.rep)
# }

# plot(mod0)
# plot(mod1)

if(F){
mods <- list(mod1, mod2, mod3, mod4, mod5, mod6)
loos <- map(mods, LOO)

LOO(mod1, mod2, mod3, mod4, mod5, mod6)
m3=as.matrix(mod3)
mcmc_intervals(m3, regex_pars = c("b_"))
m4=as.matrix(mod4)
mcmc_intervals(m4, regex_pars = c("b_"))
LOO(mod1, mod2, mod3, mod4, mod5, mod6, reloo = TRUE)
}
# plot.postpred(mod0)
# ggsave(filename=plot.filename("ppred_mod0.pdf", base = bname))
# plot.postpred(mod1)
# ggsave(filename=plot.filename("ppred_mod1.pdf", base = bname))


softmax <- function(x){
  if(is.null(dim(x))){
    dim(x) <- c(1,length(x))
  }
  r<-exp(x)/rowSums(exp(x))
  if(dim(x)[1]==1){
    r<-as.vector(r)
  }
  r
}
## probs per group/day
# value ~ effect + type + time + time:type + effect:type
# effect: placebo/nocebo
# type: control/conditioning
# time: expected/perceived
posterior_samples(mod6, pars="b_.*") -> ms
intercpt   = with(ms, cbind(b_decline_Intercept,0,b_improve_Intercept))
effectncb  = with(ms, cbind(b_decline_effectnocebo,0,b_improve_effectnocebo))
typecond   = with(ms, cbind(b_decline_typeconditioning,0,b_improve_typeconditioning))
timeperc   = with(ms, cbind(b_decline_timeperceived,0,b_improve_timeperceived))
noc.x.cond = with(ms, cbind(`b_decline_effectnocebo:typeconditioning`,0,`b_improve_effectnocebo:typeconditioning`))
prc.x.cond = with(ms, cbind(`b_decline_typeconditioning:timeperceived`,0,`b_improve_typeconditioning:timeperceived`))

labs=c("decline","neutral","improve")
## expectations
p.exp.plc.ctrl = softmax(intercpt) %>% data.frame %>% setNames(labs) 
p.exp.noc.ctrl = softmax(intercpt + effectncb) %>% data.frame %>% setNames(labs) 
p.exp.plc.cond = softmax(intercpt + typecond) %>% data.frame %>% setNames(labs) 
p.exp.noc.cond = softmax(intercpt + typecond + effectncb + noc.x.cond) %>% data.frame %>% setNames(labs) 

## perceived
p.prc.plc.ctrl = softmax(intercpt + timeperc) %>% data.frame %>% setNames(labs) 
p.prc.noc.ctrl = softmax(intercpt + timeperc + effectncb + prc.x.cond) %>% data.frame %>% setNames(labs) 
p.prc.plc.cond = softmax(intercpt + timeperc + typecond + prc.x.cond) %>% data.frame %>% setNames(labs) 
p.prc.noc.cond = softmax(intercpt + timeperc + effectncb + typecond + prc.x.cond + noc.x.cond) %>% data.frame %>% setNames(labs) 

library(ggtern)
nsamp=dim(p.exp.plc.ctrl)[1]
pexp=cbind(
  rbind(p.exp.plc.ctrl, p.exp.noc.ctrl, p.exp.plc.cond, p.exp.noc.cond),
  placebo=rep(rep(c("placebo","nocebo"), each=nsamp), 2),
  conditioning=rep(c("control","experimental"), each=nsamp*2)
  )

pall=cbind(
  rbind(p.exp.plc.ctrl, p.exp.noc.ctrl, p.exp.plc.cond, p.exp.noc.cond,
        p.prc.plc.ctrl, p.prc.noc.ctrl, p.prc.plc.cond, p.prc.noc.cond),
  placebo=rep(rep(c("placebo","nocebo"), each=nsamp), 4),
  conditioning=rep(rep(c("control","experimental"), each=nsamp*2), 2),
  time=rep(c("expected", "perceived"), each=nsamp*4)
)

hdi.lower <- function(x) hdi(x)[1]
hdi.upper <- function(x) hdi(x)[2]

pall %>% group_by(placebo,conditioning,time) %>%
  summarize_all(funs(mean,hdi.lower,hdi.upper))

  
pall %>%
  ggtern(aes(decline,neutral,improve,color=placebo))+
  #geom_point(colour='grey', alpha=0.2)+
  geom_density_tern()+theme_bw()+
  facet_grid(time~conditioning)

cache.var("mod1", base = bname)
cache.var("mod2", base = bname)
cache.var("mod3", base = bname)
cache.var("mod4", base = bname)
cache.var("mod5", base = bname)
cache.var("mod6", base = bname)
load.cache.var("mod6", bname)

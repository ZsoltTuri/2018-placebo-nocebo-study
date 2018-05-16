# This script is part of 
#
# Turi, Z., Bjørkedal, E., Gunkel, L., Antal, A., Paulus, W. & Mittner, M. (2018).
# Evidence for Cognitive Placebo and Nocebo Effects in Healthy Individuals.
#
# Computational modelling.
#
#

library(ProjectTemplate)
load.project()
theme_set(theme_bw())
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(latex2exp)
library(patchwork)

bname<-tools::file_path_sans_ext(basename(this.file.name()))
modnum=str_split(bname, "_")[[1]][1]
# enabls --force for stallo etc.
options <- commandArgs(trailingOnly = TRUE)
if( "--force" %in% options)
  uncache.all()

#stop()

### MCMC parameters iterations
n.chains=8
n.cores=8
n.iter=2000
n.warmup=1000

## priors for mu_alphaG, mu_alphaL and mu_beta
n=1000000
mu_p=rnorm(n,-.5,.6)
sigma=extraDistr::rhcauchy(n, .01)
err=rnorm(n,0,1)
alpha=pnorm( mu_p + sigma*err)

mu_p2=rnorm(n,-1.5,.8)
#sigma2=extraDistr::rhcauchy(n, 1)
sigma2=extraDistr::rhnorm(n, .3)
err=rnorm(n,0,1)
beta=exp( mu_p2 + sigma2*err)

sigma_b=extraDistr::rhnorm(n,1)
b = rnorm(n,0,sigma_b)

tibble(alpha,beta,b) %>% filter(beta<2, b<3, b>-3) %>%
  gather(var,val) %>% mutate(var=fct_relevel(var,"alpha","beta","b")) %>% ggplot(aes(x=val))+stat_density(geom="line")+#geom_density()+
  facet_wrap(~var, scales="free") -> p
p + labs(x="parameter value",  y="density")
ggsave(plot.filename("prior.pdf", bname), width=9,height=2)


##================================================================================================
# model fitting
if(!is.cached.var("mod", base=bname)){
  #n.cores=1;  n.chains=1; n.iter=100; n.warmup=10
  mod = stan(file=sprintf("src/%s.stan", bname), data = stan.data, cores=n.cores, chains = n.chains, iter = n.iter, warmup = n.warmup)
  cache.var("mod", bname)
} else {
  mod <- load.cache.var("mod",bname)
}

#stop()
##================================================================================================
## diagnostics
pdf(file=plot.filename("diagnostics.pdf", base=bname), onefile=TRUE)
mcmc_rhat(rhat(mod)) %>% print
mcmc_neff(neff_ratio(mod)) %>% print


npages=30
modm = as.matrix(mod)
moda= as.array(mod)
vnames=colnames(modm)
nvar=dim(modm)[2]
for(p in split(1:nvar, ceiling(seq_along(1:nvar)/(nvar/npages)))){
  mcmc_trace(moda[,,p], np = nuts_params(mod)) %>% print
}
dev.off()


##================================================================================================
# parameters
params=c("alpha_g", "alpha_l", "beta")
mcmc_hist(modm, pars=paste("mu",params, sep="_"))
ggsave(filename = plot.filename("grouppars.pdf", base=bname))

mcmc_hist(modm, regex_pars = "sigma*")
#stan_hist(mod, pars=paste("mu",params, sep="_"))
ggsave(filename = plot.filename("sigmas.pdf", base=bname))

pl <- lapply(params, function(p) mcmc_intervals(as.matrix(mod), regex_pars = sprintf("%s\\[.*",p)))

pdf(file=plot.filename("indpars.pdf", base=bname), width = 20)
multiplot(plotlist=pl, cols = length(pl))
dev.off()



for( p in as.data.frame(combn(params, 2))){
  p=as.character(p)
  samples %>%
    ggplot(aes_string(x=p[1], y=p[2]))+geom_point(alpha=0.1)+facet_wrap(~subj, scales="free") -> pp
  
  ggsave(plot=pp, filename = plot.filename(sprintf("scatter_%s_%s.png", p[1], p[2]), base=bname), width=40, height=40)
}

## effects
mcmc_intervals(as.matrix(mod), regex_pars="b_.*")
ggsave(plot.filename("fixed_effects.pdf",bname))
##================================================================================================
# model selection
loo.tmp=if.cached.load(sprintf("loo%s",modnum), {
  loo(extract_log_lik(mod))
})
assign(sprintf("loo%s",modnum), loo.tmp)
##================================================================================================
# posterior predictive
source(sprintf("src/%s_ppred.R", bname))

nrep=100

ppred=if.cached.load("ppred",{
  get.ppred(mod,nrep,stan.data)
}, base = bname)

presp=if.cached.load("presp",{
  get.presp(mod,nrep,stan.data)
}, base = bname)


## ppred per group and day
# -------------------------
ppred %>% filter(nrep==0) %>% group_by(nrep,group,day,trial) %>%
  summarize(pgo=mean(accuracy)) -> dd

ppred %>% filter(nrep==00) %>% group_by(nrep,subj,day,pair) %>%
  mutate(ptrial=1:n()) %>% ungroup %>%
  group_by(nrep,group,day,pair,ptrial) %>%
  summarize(pcorr=mean(accuracy)) %>% ungroup %>%
  mutate(nrep=factor(nrep),day=factor(day),pair=factor(pair)) -> dd

ppred %>% filter(nrep>0) %>% group_by(nrep,subj,day,pair) %>%
  mutate(ptrial=1:n()) %>% ungroup %>%
  group_by(nrep,group,day,pair,ptrial) %>%
  summarize(pcorr=mean(accuracy)) %>% ungroup %>%
  mutate(nrep=factor(nrep),day=factor(day),pair=factor(pair))-> ppred.grp

ppred.grp %>%
  ggplot(aes(x=ptrial,y=pcorr))+
  #stat_summary(fun.data = mean_cl_normal, geom="ribbon")+
  geom_line(aes(group=nrep),color="grey",alpha=0.3)+
  geom_line(data=dd, color="red")+
  facet_grid(day~group+pair)

ggsave(plot.filename("ppred_pcorr.pdf", base=bname), width=20,height=10)


ppred.grp %>%
  ggplot(aes(x=ptrial,y=pcorr))+
  #geom_line(aes(group=nrep),color="grey",alpha=0.3)+
  stat_summary(fun.data = mean_hdi, fill="grey", geom="ribbon")+
  stat_summary(fun.y = mean, geom="line",color="blue")+
  geom_line(data=dd, color="red")+
  facet_grid(day~group+pair)

ggsave(plot.filename("ppred_pcorr_ribbon.pdf", base=bname), width=20,height=10)


## individual parameter estimates
##===========================================
mod %>%
  spread_samples(alpha_g[subj], alpha_l[subj], beta[subj]) -> samples

samples %>% gather(var,val,alpha_g, alpha_l,beta) %>%
  group_by(subj,var) %>%
  summarize(mean=mean(val), lower=hdi(val)[1], upper=hdi(val)[2]) %>%
  ggplot(aes(y=mean,x=subj,ymin=lower,ymax=upper,color=var))+
  geom_pointrange(position=position_dodge(width=0.9))+coord_flip()
ggsave(plot.filename("params_subj_overlaid.pdf", bname), height=15, width=7)


mod %>%
  spread_samples(alpha_g_day2[subj], alpha_l_day2[subj], beta_day2[subj]) -> samples_day2

samples_day2 %>% gather(var,val,alpha_g_day2, alpha_l_day2,beta_day2) %>%
  group_by(subj,var) %>%
  summarize(mean=mean(val), lower=hdi(val)[1], upper=hdi(val)[2]) %>% 
  gather(ix,val, mean,lower,upper,-subj) %>%
  unite(temp,var,ix) %>%
  spread(temp, val) -> df
attach(stan.data)
mat=cbind(group_p_cntrl,group_n_cntrl, group_p_cond, group_n_cond)
groupfac=factor(mat%*%(1:4), labels = c("nhg",colnames(mat)) )
detach("stan.data")

transp=0.2
cbind(df,group=groupfac) %>%
  ggplot(aes(x=alpha_g_day2_mean, xmin=alpha_g_day2_lower, xmax=alpha_g_day2_upper,
             y=alpha_l_day2_mean, ymin=alpha_l_day2_lower, ymax=alpha_l_day2_upper,
             color=group))+
  geom_point()+geom_errorbar(alpha=transp)+geom_errorbarh(alpha=transp)+xlim(0,1)+ylim(0,1)+
  theme(legend.position="none")-> p1

cbind(df,group=groupfac) %>%
  ggplot(aes(x=alpha_g_day2_mean, xmin=alpha_g_day2_lower, xmax=alpha_g_day2_upper,
             y=beta_day2_mean, ymin=beta_day2_lower, ymax=beta_day2_upper,
             color=group))+
  geom_point()+geom_errorbar(alpha=transp)+geom_errorbarh(alpha=transp)+xlim(0,1)+
  theme(legend.position="none")-> p2
cbind(df,group=groupfac) %>%
  ggplot(aes(x=alpha_l_day2_mean, xmin=alpha_l_day2_lower, xmax=alpha_l_day2_upper,
             y=beta_day2_mean, ymin=beta_day2_lower, ymax=beta_day2_upper,
             color=group))+
  geom_point()+geom_errorbar(alpha=transp)+geom_errorbarh(alpha=transp)+xlim(0,1) -> p3

p<-(p1+p2+p3)
ggsave(plot.filename("indpars_pairs_hdi.pdf", bname),width=14,height=4.5,plot=p)


## posterior mean scatter
#===========================

as.data.frame(mod) %>% 
  mutate(mu_alg_day2_nhg=pnorm(`mu_p[1]`+b_alg_day2),
         mu_alg_day2_pcond=pnorm(`mu_p[1]`+b_alg_day2+b_alg_pcond),
         mu_alg_day2_ncond=pnorm(`mu_p[1]`+b_alg_day2+b_alg_ncond),
         mu_alg_day2_pctrl=pnorm(`mu_p[1]`+b_alg_day2+b_alg_pctrl),
         mu_alg_day2_nctrl=pnorm(`mu_p[1]`+b_alg_day2+b_alg_nctrl),
         
         mu_all_day2_nhg=pnorm(`mu_p[2]`+b_all_day2),
         mu_all_day2_pcond=pnorm(`mu_p[2]`+b_all_day2+b_all_pcond),
         mu_all_day2_ncond=pnorm(`mu_p[2]`+b_all_day2+b_all_ncond),
         mu_all_day2_pctrl=pnorm(`mu_p[2]`+b_all_day2+b_all_pctrl),
         mu_all_day2_nctrl=pnorm(`mu_p[2]`+b_all_day2+b_all_nctrl),
         
         mu_be_day2_nhg=exp(`mu_p[3]`+b_be_day2),
         mu_be_day2_pcond=exp(`mu_p[3]`+b_be_day2+b_be_pcond),
         mu_be_day2_ncond=exp(`mu_p[3]`+b_be_day2+b_be_ncond),
         mu_be_day2_pctrl=exp(`mu_p[3]`+b_be_day2+b_be_pctrl),
         mu_be_day2_nctrl=exp(`mu_p[3]`+b_be_day2+b_be_nctrl)
         ) %>% select(starts_with("mu_"), -starts_with("mu_p")) %>% 
  rename(mu_alg_day1_all=mu_alpha_g, mu_all_day1_all=mu_alpha_l,mu_be_day1_all=mu_beta) %>% 
  gather(var,val) %>%
  separate(var, c("mu","parameter", "day", "group"), "_", remove=F) %>%
  group_by(parameter, day, group) %>%
  summarize(mean=mean(val), lower=hdi(val)[1], upper=hdi(val)[2]) %>%
  gather(ix,val, mean,lower,upper,-day,-group) %>%
  unite(temp,parameter,ix) %>%
  spread(temp, val) -> df

df %>%
  ggplot(aes(x=alg_mean,
             y=all_mean,
             color=group))+
  geom_point(aes(shape=day),size=3)+
  geom_errorbar(aes(ymin=all_lower,ymax=all_upper), width=0,alpha=transp)+
  geom_errorbarh(aes(xmin=alg_lower,xmax=alg_upper), height=0,alpha=transp)+
  labs(x="alphaG", y="alphaL")+
  theme(legend.position="none")-> p1
df %>%
  ggplot(aes(x=alg_mean,xmin=alg_lower,xmax=alg_upper,
             y=be_mean,ymin=be_lower,ymax=be_upper,
             color=group))+
  geom_point(aes(shape=day),size=3)+
  geom_errorbar(width=0,alpha=transp)+geom_errorbarh(height=0,alpha=transp)+
  labs(x="alphaG", y="beta")+
  theme(legend.position="none")-> p2
df %>%
  ggplot(aes(x=be_mean,xmin=be_lower,xmax=be_upper,
             y=all_mean,ymin=all_lower,ymax=all_upper,
             color=group))+
  geom_point(aes(shape=day),size=3)+
  geom_errorbar(width=0,alpha=transp)+geom_errorbarh(height=0,alpha=transp)+
  labs(x="beta", y="alphaL")-> p3

p<-(p1+p2+p3)
ggsave(plot.filename("grouppars_pairs_hdi.pdf", bname),width=10,height=3,plot=p)

### OVERLAY PLOT
#==================================
result=load.cache.var("result", base="optimal_pars")

#good.shapes=c(15,16,17,18,24,25)
good.shapes=c(0,1,2,6,17,25)
good.shapes=c(1,25,6,0,17,2)# ,1,2,6,17,25)
symbol.size=6
bg.transp=.8
transp=0.2

ggplot(result,aes(x=alphaG, y=alphaL, fill=result))+
  geom_raster(alpha=bg.transp)+  scale_fill_gradientn(colours = jet.colors(7))+
  geom_point(aes(x=alg_mean,y=all_mean,fill=NULL,color=NULL,shape=group),data=df,size=symbol.size,fill="black")+scale_shape_manual(values=good.shapes)+
  geom_errorbar(aes(x=alg_mean,y=all_mean,ymin=all_lower,ymax=all_upper,fill=NULL,color=NULL), width=0,alpha=transp,data=df)+
  geom_errorbarh(aes(x=alg_mean,y=all_mean,xmin=alg_lower,xmax=alg_upper,fill=NULL,color=NULL), height=0,alpha=transp,data=df)+
  labs(x=TeX("$\\alpha_G$"), y=TeX("$\\alpha_L$"))+
  theme(legend.position="none")+
  coord_cartesian(xlim = c(0.055, 0.245), c(0.025,0.175)) -> p1

ggplot(result,aes(x=alphaG, y=beta, fill=result))+
  geom_raster(alpha=bg.transp)+  scale_fill_gradientn(colours = jet.colors(7))+
  geom_point(aes(x=alg_mean,y=be_mean,fill=NULL,color=NULL,shape=group),data=df,size=symbol.size,fill="black")+scale_shape_manual(values=good.shapes)+
  geom_errorbar(aes(x=alg_mean,y=be_mean,ymin=be_lower,ymax=be_upper,fill=NULL,color=NULL), width=0,alpha=transp,data=df)+
  geom_errorbarh(aes(x=alg_mean,y=be_mean,xmin=alg_lower,xmax=alg_upper,fill=NULL,color=NULL), height=0,alpha=transp,data=df)+
  labs(x=TeX("$\\alpha_G$"), y=TeX("$\\beta$"))+
  theme(legend.position="none") +
  coord_cartesian(xlim = c(0.055, 0.245), c(0.108,0.342)) -> p2

ggplot(result,aes(y=beta, x=alphaL, fill=result))+
  geom_raster(alpha=bg.transp)+  scale_fill_gradientn(colours = jet.colors(7))+
  geom_point(aes(y=be_mean,x=all_mean,fill=NULL,color=NULL,shape=group),data=df,size=symbol.size,fill="black")+scale_shape_manual(values=good.shapes)+
  geom_errorbar(aes(y=be_mean,x=all_mean,ymin=be_lower,ymax=be_upper,fill=NULL,color=NULL), width=0,alpha=transp,data=df)+
  geom_errorbarh(aes(y=be_mean,x=all_mean,xmin=all_lower,xmax=all_upper,fill=NULL,color=NULL), height=0,alpha=transp,data=df)+
  labs(y=TeX("$\\beta$"), x=TeX("$\\alpha_L$"))+
  #theme(legend.position="none") +
  coord_cartesian(xlim = c(0.025, 0.175), ylim= c(0.108,0.342)) -> p3

p<-(p1+p2+p3) & theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
figdims=c(9,3)*1.3
ggsave(plot.filename("grouppars_optpars_overlay.pdf", bname), plot=p,width=figdims[1], height=figdims[2])


## table
vnames=as.data.frame(mod) %>% colnames
bnames=vnames[str_detect(vnames, "b_.*")]
tab <- summary(mod, pars = bnames, probs = c(0.05, 0.95))$summary
data.frame(tab) %>% select(-se_mean,-sd, -n_eff,-Rhat) %>% setNames(c("mean","lower","upper")) %>%
  mutate(variable=rownames(.)) %>% separate(variable, into=c("b","var","group")) %>%
  mutate(text=sprintf("%.2f [%.2f,%.2f]", mean, lower,upper)) %>% select(var,text,group) %>% spread(var,text)
  

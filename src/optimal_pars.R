library(ProjectTemplate)
load.project()

bname<-tools::file_path_sans_ext(basename(this.file.name()))
#stop()
# enabls --force for stallo etc.
options <- commandArgs(trailingOnly = TRUE)
if( "--force" %in% options ){
  uncache.all()
} 


library(Rcpp)
Rcpp::sourceCpp(file.path('lib', 'clikelihood.cpp'))
ntrials=80

alphaGs=seq(0.05,  0.25,by=0.01)
alphaLs=seq(0.02,  0.18,by=0.01)
betas  =seq(0.10,  0.35,by=0.01)

length(alphaGs)*length(alphaLs)*length(betas)

## parallel?
library(parallel)
ncpu=20
nrep=10000#500000

pars<-expand.grid(alphaG=alphaGs, alphaL=alphaLs, beta=betas)

if(is.cached.var('result', bname)){
  result=load.cache.var('result', bname)
} else {
  res=mclapply(data.frame(t(pars)), function(x){ mean(replicate(nrep, csimandeval_posneg(x[1], x[2], x[3], ntrials))) }, mc.cores=ncpu)
  result=cbind(pars, result=unlist(res))
  cache.var('result', bname)  
}

optpars <- unlist(with(result, pars[which.max(result),]))
show(optpars)

breaks<-c(c(0.7, 0.75, 0.8), seq(.8, .89, length.out = 6),
          seq(.89, .9, length.out=20))


ggplot(result)+geom_raster(aes(x=alphaG, y=alphaL, fill=result))+theme_bw()+
  scale_fill_gradientn(colours = jet.colors(7))
ggplot(result)+geom_raster(aes(x=alphaG, y=beta, fill=result))+theme_bw()+
  scale_fill_gradientn(colours = jet.colors(7))
ggplot(result)+geom_raster(aes(x=beta, y=alphaL, fill=result))+theme_bw()+
  scale_fill_gradientn(colours = jet.colors(7))


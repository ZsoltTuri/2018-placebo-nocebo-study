printf <- function(s, ...){
  cat(sprintf(s, ...))
}

log.filename <- function(log.name, base=NULL,sub="",overwrite=T){
  if(is.null(base))
    bname<-tools::file_path_sans_ext(basename(this.file.name()))
  else
    bname<-tools::file_path_sans_ext(basename(base))
  
  dir.create(file.path('logs',sub, bname), showWarnings = FALSE)
  fname<-file.path('logs', sub,bname, log.name)
  if( file_exists(fname) && overwrite){
    cat("> overwrite ", fname, "\n");
    file_delete(fname)
  }
  cat("> ", fname, "\n")
  fname  
}


get.mods <- function(mod.nums){
  fpath="cache/vars/"
  mods <- list()
  for(mod.num in mod.nums){
    mod.fname <- list.files(fpath, pattern=sprintf("^%s_.*_mod.RData", mod.num))
    bname=tools::file_path_sans_ext(mod.fname)
    bname=substr(bname, 1, (nchar(bname)-4))
    mods[[mod.num]] <- load.cache.var("mod",bname)
  }
  mods
}

get.loos <- function(mod.nums){
  fpath="cache/vars/"
  loos <- list()
  for(mod.num in mod.nums){
    mod.fname <- list.files(fpath, pattern=sprintf("^%s_.*_mod.RData", mod.num))
    bname=tools::file_path_sans_ext(mod.fname)
    bname=substr(bname, 1, (nchar(bname)-4))
    mod <- load.cache.var("mod",bname)
    loos[[mod.num]] <- loo(extract_log_lik(mod))
  }
  loos
}

delta.loo <- function(loos){
  mods.ordered<-rev(rownames(compare(x=loos)))
  nmods=length(loos) 
  delta=list()
  for(i in 1:(nmods-1)){
    delta[[sprintf("%s-%s", mods.ordered[i], mods.ordered[i+1])]]<-compare(loos[[ mods.ordered[i] ]], loos[[ mods.ordered[i+1] ]])
  }
  do.call(rbind, delta)
}

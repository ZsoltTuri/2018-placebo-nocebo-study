## caching single variables based on filename

## http://stackoverflow.com/a/15373917/972791
this.file.name <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

formula.name.fname <- function(frm){
  frm_str=paste(format(frm), sep="", collapse="")
  gsub("\\s+", " ", str_trim(frm_str))
  str_replace_all(
    str_replace_all(
      str_replace_all(
        str_replace_all(
          str_replace_all(frm_str, "\\|", "by"),
          "~", "eq"),
        "\\+", "plus"),
      "\\*", "x"), 
    "[^[:alnum:]]", "_") -> frm_str
  gsub("_+", "_", frm_str) -> frm_str
  gsub("_$", "", frm_str) -> frm_str
  frm_str
}

cache.get.fname <- function(varname, base){
  if(is.null(base))
    bname<-tools::file_path_sans_ext(basename(this.file.name()))
  else
    bname<-tools::file_path_sans_ext(basename(base))
  fname<-file.path('cache', 'vars', sprintf("%s_%s.RData",bname,varname))
  fname  
}

is.cached.var <- function(varname, base=NULL){
  fname<-cache.get.fname(varname, base)
  file.exists(fname)
}
cache.var <- function(varname, base=NULL){
  fname<-cache.get.fname(varname, base)
  printf("CACHE> saving %s to %s\n", varname, fname)
  dir.create(dirname(fname), showWarnings = FALSE)
  save(list = varname, envir = .GlobalEnv, file = fname)
}
uncache.var <- function(varname, base=NULL){
  fname<-cache.get.fname(varname, base)
  cat(sprintf("Deleting %s\n", fname))
  unlink(fname)
}

uncache.all <- function(base=NULL){
  if(is.null(base))
    bname<-tools::file_path_sans_ext(basename(this.file.name()))
  else
    bname<-tools::file_path_sans_ext(basename(base))
  fnames=list.files(file.path('cache', 'vars'), pattern = sprintf("%s*", bname), full.names=T)
  for(fname in fnames){
    cat(sprintf("Deleting %s\n", fname))
    unlink(fname)
  }
}

if.cached.load <- function(vname, expr,base=bname){
  if(!is.cached.var(vname, base=bname)){
    val=eval.parent(expr)
    assign(vname, val, envir=.GlobalEnv)
    cache.var(vname, bname)
  } else {
    val <- load.cache.var(vname,bname)
  }
  return(val)
}

load.cache.var <- function(varname, base=NULL){
  fname<-cache.get.fname(varname, base)
  printf("CACHE> loading %s from %s\n", varname, fname)
  load(fname)
  return(eval(parse(text=varname)))
}

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))



plot.filename <- function(plot.name, base=NULL,sub=""){
  if(is.null(base))
    bname<-tools::file_path_sans_ext(basename(this.file.name()))
  else
    bname<-tools::file_path_sans_ext(basename(base))
  
  
  
  dir.create(file.path('graphs',sub, bname), showWarnings = FALSE)
  fname<-file.path('graphs', sub,bname, plot.name)
  cat("> ", fname, "\n")
  fname  
}



# Adds a highest density ellipse to an existing plot
# xy: A matrix or data frame with two columns.
#     If you have to variables just cbind(x, y) them.
# coverage: The percentage of points the ellipse should cover
# border: The color of the border of the ellipse, NA = no border
# fill: The filling color of the ellipse, NA = no fill
# ... : Passed on to the polygon() function
#
#http://www.sumsar.net/blog/2014/11/how-to-summarize-a-2d-posterior-using-a-highest-density-ellipse/
add_hd_ellipse <- function(xy, coverage, border = "blue", fill = NA, ...) {
  fit <- MASS::cov.mve(xy, quantile.used = round(nrow(xy) * coverage), nsamp=1000)
  points_in_ellipse <- xy[fit$best, ]
  ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
  polygon(ellipse_boundary, border=border, col = fill, ...)
}

stan_trace_multipage <- function(fit, N, pars=NULL, ...){
  if(is.null(pars)){
    pars=fit@model_pars
  } 
  nplots=sum(unlist(lapply(fit@par_dims, function(x){ if(is.null(x)){1}else{prod(x)}})))
  npages=ceiling(nplots/N)
  nplots.per.page=ceiling(nplots/npages)
  for( i in 1:2){
    fromix=((i-1)*nplots.per.page)+1
    toix=fromix+nplots.per.page
    if(toix>nplots){
      toix=nplots
    }
    stan_trace(fit, pars=pars[fromix:toix], ...)
  }
}

textplot <- function(text){
  ggplot() + 
    annotate("text", x = 4, y = 25, size=8, label = text, hjust=0.5) + 
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          plot.background=element_blank())
}

dotplot.old <- function(ms, real=NULL, exclude_vars=NULL){
  pm=posterior.mode(ms)
  vars=factor(names(pm))
  p=data.frame(variable=vars, pmode.hdi=pm) %>% cbind(hdi(ms)) %>% 
    filter(!(variable %in% exclude_vars)) %>%
    ggplot(aes(x=pmode.hdi,y=variable,xmin=lower,xmax=upper))+
    geom_errorbarh(height=.1)+geom_point(size=3)
  if(!is.null(real)){
    p=p+geom_point(aes(x=real), color="red", size=4)
  }
  p
}

dotplot <- function(..., real=NULL, exclude_vars=NULL){
  nms <- setdiff(as.character(match.call(expand.dots=TRUE)), 
                 as.character(match.call(expand.dots=FALSE)))
  pm<-NULL
  i=1
  for(ms in list(...)){
    dset=nms[i]
    vars=factor(colnames(as.matrix(ms)))
    pm <- rbind(pm, cbind(data.frame(dset=dset, variable=vars, pmode.hdi=posterior.mode(ms)),
                          hdi(ms)))
    i=i+1
  }
  
  p=pm %>% mutate(i=factor(i)) %>%
    filter(!(variable %in% exclude_vars)) %>%
    ggplot(aes(x=pmode.hdi,y=variable,color=dset,xmin=lower,xmax=upper))+
    geom_errorbarh(height=.2, alpha=.6)+geom_point(size=3)+
    scale_colour_manual(values=c("black","blue", "green", 'magenta', 'orange'))
  if(!is.null(real)){
    p=p+geom_point(aes(x=real), color="red", size=4)
  }
  p
}


#
# Example for beta-distribution: 
#   ms=ms[,c("alpha", "beta")]
#   plot.postpred( ms, dbeta, c("alpha", "beta"), 100, xlim=c(0,1))
plot.postpred <- function(ms, fun, par.names, npostpred, xlim=c(0,1)){
  m=as.matrix(ms)
  mpars<-m[sample(1:dim(m)[1], npostpred, replace=T),]
  b.pars=posterior.mode(ms)
  args=as.list(b.pars)
  names(args) <- par.names
  p<-ggplot(data.frame(x = xlim), aes(x)) + 
    stat_function(fun = fun, args=args, size=2, colour='red')
  for( i in 1:npostpred){
    args=as.list(m[i,])
    names(args)<-par.names
    p<-p+stat_function(fun = fun, alpha=0.5, color='grey',args=args)
  }
  p
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## MATLAB jet colors
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


## http://stackoverflow.com/a/8197703
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

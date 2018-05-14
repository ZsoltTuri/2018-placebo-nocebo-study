rbind(day1.learn,day2.learn) %>%
  mutate(symbl_pair=fct_relevel(symbl_pair, "AB","CD","EF")) -> d

d %>% group_by(subj) %>% summarize(group=first(group)) -> gd

stan.data=list(
  N=length(unique(d$subj)),
  T=length(unique(d$trial))*2, ## both days one after the other
  group_n_cond=as.integer(gd$group=="n_cond"),
  group_p_cond=as.integer(gd$group=="p_cond"),
  group_n_cntrl=as.integer(gd$group=="n_cntrl"),
  group_p_cntrl=as.integer(gd$group=="p_cntrl"),
  
  pair=do.call(rbind,
               d %>% split(d$subj) %>% map(function(df){ as.integer(df$symbl_pair)})
  ),
  day2=do.call(rbind,
               d %>% split(d$subj) %>% map(function(df){ as.integer(df$day==2)})
  ),
  outcome=do.call(rbind,
                  d %>% split(d$subj) %>% map(function(df){ as.integer(df$reward)})
  ),
  accuracy=do.call(rbind,
                   d %>% split(d$subj) %>% map(function(df){ as.integer(df$accuracy)})
  )
)
 
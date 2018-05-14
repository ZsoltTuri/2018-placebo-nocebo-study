source("data/groups.R")

day2.learn = NULL
group.pathes = c('1_nocebo_cond','2_placebo_cond','3_nocebo_control','4_placebo_control','5_natural_history')

for(group in group.pathes){
  dpath <- file.path('data', 'raw', group, 'day2_learning')
  fnames <- list.files(dpath, pattern = "*.csv")
  
  learn <- NULL
  for(fname in fnames){
    d <- read.table(file.path(dpath, fname), header = T, sep = ",") %>% 
      mutate(buttonBox_2.rt = as.numeric(buttonBox_2.rt)) %>%
      select(participant = Participant, group = Group, symbl_pair = pair, symbl_position = position, fb_corr, fb_incorr, accuracy = buttonBox_2.corr, reaction_time = buttonBox_2.rt)  %>%
      subset(symbl_pair == "AB" |  symbl_pair == "CD" | symbl_pair == "EF") %>%
      mutate(trial = 1:n(), 
             symbl_pair = droplevels(symbl_pair), 
             fb_corr = as.numeric(str_detect(fb_corr, "posf")), 
             fb_incorr = as.numeric(str_detect(fb_incorr, "posf")),
             reward = ifelse(accuracy > 0, fb_corr, fb_incorr))  %>%
      select(-fb_corr, -fb_incorr)
    learn <- rbind(learn, d) 
  }
  learn$participant <- factor(rep((1:16), each = 240))
  day2.learn = rbind(day2.learn,learn)
}
day2.learn$group <- factor(rep(groups, each = 16 * 240))
day2.learn$day = factor(2)
rm(learn, d)

day2.learn %>%
  mutate(subj=factor(paste(group, participant, sep="_"))) -> day2.learn


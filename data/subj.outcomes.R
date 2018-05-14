source("data/groups.R")

subj.outcomes = NULL
dpath <- file.path('data', 'raw')
fnames <- head(list.files(dpath, pattern = "*.csv"), 4)

for(fname in fnames){
  d <- read.csv(file.path(dpath, fname), header = T, sep = ";") %>% 
    mutate(PID = factor(PID), Group = factor(Group)) %>% 
    select(participant = PID, group = Group, expected = ExpDirection, perceived =  PerceivedDirction)
  subj.outcomes = rbind(subj.outcomes, d)
}
subj.outcomes$expected <- recode(subj.outcomes$expected, "-1" = "decline", "0" = "neutral", "1" = "improve")
subj.outcomes$perceived <- recode(subj.outcomes$perceived, "-1" = "decline", "0" = "neutral", "1" = "improve")
subj.outcomes$effect <- factor(rep(c("nocebo","placebo"), each = 16, times = 2))
subj.outcomes$type <- factor(rep(c("conditioning","control"), each = 32))
rm(d)

subj.outcomes %>%
  mutate(group=fct_recode(group,
                          n_cond="nocebo_cond",
                          p_cond="placebo_cond",
                          n_cntrl="nocebo_control",
                          p_cntrl="placebo_control")) %>%
  mutate(subj=factor(paste(group, participant, sep="_"))) -> subj.outcomes




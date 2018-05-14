source("data/groups.R")

demographic = NULL
dpath <- file.path('data', 'raw')
fnames <- list.files(dpath, pattern = "*.csv")

for(fname in fnames){
  d <- read.csv(file.path(dpath, fname), header = T, sep = ";") %>% 
    mutate(PID = factor(PID), Group = factor(Group)) %>% 
    select(participant = PID, group = Group, age = Age)
  demographic = rbind(demographic, d)
}
rm(d)


demographic %>%
  mutate(group=fct_recode(group,
                          n_cond="nocebo_cond",
                          p_cond="placebo_cond",
                          n_cntrl="nocebo_control",
                          p_cntrl="placebo_control")) %>%
  mutate(subj=factor(paste(group, participant, sep="_"))) -> demographic



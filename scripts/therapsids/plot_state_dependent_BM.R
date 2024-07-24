## load log files

library(ggplot2)
library(tidyverse)
library(tidyr)

#df1 <- read.table("output/state_dependent_BM_character_1_run_1.log", header = TRUE)
dfs <- list()
for (i in 1:4){
  df1 <- read.table(paste0("output/state_dependent_BM_character_", i, ".log"), header = TRUE)
  
  df1 <- df1 %>%
    select(-starts_with("state_branch")) %>%
    select(-starts_with("num_changes")) %>%
    select(-starts_with("branch_rates"))
  
  df1$zeta_diff <- df1$zeta.1. - df1$zeta.2.
  df1$character <- factor(i)
  dfs[[i]] <- df1
}
df2 <- bind_rows(dfs)

zeta_plots <- list()
for (i in 1:4){
  p <- dfs[[i]] %>%
    ggplot(aes(x = zeta_diff)) +
    #geom_density(alpha = 0.5, fill = "gray") +
    geom_histogram() +
    theme_classic()
  zeta_plots[[i]] <- p
}






ggplot(df1, aes(x = zeta.1.)) +
  geom_histogram()

sum(df1$zeta.1. > df1$zeta.2.) / nrow(df1)

hist(df1$zeta.1. - df1$zeta.2.)


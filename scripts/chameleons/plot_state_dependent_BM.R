## load log files

library(ggplot2)
library(tidyverse)
library(tidyr)

#df1 <- read.table("output/state_dependent_BM_character_1_run_1.log", header = TRUE)
df1 <- read.table(paste0("output/chameleon_BM_state_dependent_run_1.log"), header = TRUE)
df2 <- read.table(paste0("output/chameleon_BM_state_dependent_run_2.log"), header = TRUE)

df <- bind_rows(df1, df2)

df <- df %>%
  select(-starts_with("state_branch")) %>%
  select(-starts_with("num_changes")) %>%
  select(-starts_with("branch_rates"))

df$zeta_diff <- df$zeta.1. - df$zeta.2.
df$character <- factor(i)


zeta_plot <- df %>%
  ggplot(aes(x = zeta_diff)) +
  #geom_density(alpha = 0.5, fill = "gray") +
  geom_histogram() +
  theme_classic()




ggplot(df, aes(x = zeta.1.)) +
  geom_histogram()

sum(df$zeta.1. > df$zeta.2.) / nrow(df)

hist(df$zeta.1. - df$zeta.2.)


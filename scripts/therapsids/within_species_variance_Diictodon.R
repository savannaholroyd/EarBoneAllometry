library(tibble)
library(dplyr)

df <- read.csv("data/DiictodonError.csv") %>% 
  as_tibble()


species_df <- df %>%
  group_by(Species) %>%
  summarize(
    "logJaw_sample_size" = length(logJaw),
    "logRL_sample_size" = length(logRL),
    "logRL_mean" = mean(logRL),
    "logJaw_mean" = mean(logRL),
    "logRL_var" = var(logRL),
    "logJaw_var" = var(logJaw),
  ) %>%
  mutate(
    "logJaw_varmean" = logJaw_var / logJaw_sample_size,
    "logRL_varmean" = logRL_var / logRL_sample_size,
         )

library(tibble)
library(dplyr)

df <- read.csv("data/SynapsidError.csv") %>% 
  as_tibble()

df$logRL <- (df$logRL1 + df$logRL2) / 2.0


species_df <- df %>%
  group_by(Species) %>%
  summarize(
    "logRL_mean" = mean(logRL),
    "logRL_var" = var(logRL),
    "logRL_sample_size" = length(logRL)
  ) %>%
  mutate(
    "logRL_varmean" = logRL_var / logRL_sample_size,
  )

species_df

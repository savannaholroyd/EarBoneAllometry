library(tibble)
library(dplyr)

df <- read.csv("data/chameleons/error measurement.csv") %>% 
  as_tibble()

df$logpter_specimen_mean <- (df$logpter1 + df$logpter2 + df$logpter3) / 3.0

species_df <- df %>%
  group_by(Species) %>%
  summarize(
    "logpter_mean" = mean(logpter_specimen_mean),
    "logpter_var" = var(logpter_specimen_mean),
    "sample_size" = length(logpter_specimen_mean)
  ) %>%
  mutate("logpter_varmean" = logpter_var / sample_size)

## read the previous data
df2 <- read.csv("data/chameleons/largest_chameleons.csv") %>%
  as_tibble()

## why are "logpter_mean" and "logpter" not the same across the two files?
## logpter is consistently different
dplyr::left_join(species_df, df2, by = "Species")






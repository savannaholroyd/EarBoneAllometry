library(tibble)
library(dplyr)

df <- read.csv("data/chameleons/masterAll.csv") %>% 
  as_tibble()

## assume that Thoehneli => Thoehnelii and Tjacksoni => Tjacksonii
df$Species <- gsub("Thoehneli", "Thoehnelii", df$Species)
df$Species <- gsub("Tjacksoni", "Tjacksonii", df$Species)

species_df <- df %>%
  group_by(Species) %>%
  summarize(
    "logpter_mean" = mean(logpter),
    "logpter_var" = var(logpter),
    "logpter_sample_size" = length(logpter),
    "logbsl_mean" = mean(logbsl),
    "logbsl_var" = var(logbsl),
    "logbsl_sample_size" = length(logbsl)
  ) %>%
  mutate(
    "logpter_varmean" = logpter_var / logpter_sample_size,
    "logbsl_varmean" = logbsl_var / logbsl_sample_size
  )


## calculate "overall variance", weighted by sample size
## assume it is this value for the species with just 1 specimen
logpter_var <- sum(species_df$logpter_var * species_df$logpter_sample_size) / sum(species_df$logpter_sample_size)
logbsl_var <- sum(species_df$logbsl_var * species_df$logbsl_sample_size) / sum(species_df$logbsl_sample_size)


## read the previous data
df2 <- read.csv("data/chameleons/largest_chameleons.csv") %>%
  as_tibble()

df3 <- df2 %>%
  mutate(
    "logpter_mean" = logpter,
    "logpter_sample_size" = 1,
    "logpter_varmean" = NaN,
    "logpter_var" = NaN,
    "logbsl_mean" = logbsl,
    "logbsl_sample_size" = 1,
    "logbsl_var" = NaN,
    "logbsl_varmean" = NaN,
  ) %>%
  select(-logbsl, -logpter, -Specimen, -Regime) %>%
  filter(!(Species %in% species_df$Species)) %>%
  mutate(
    "logpter_varmean" = !!logpter_var / 1.0, ## exclamation mark !! is to use the `logpter_var` variable from the global environment
    "logbsl_varmean" = !!logbsl_var / 1.0,
  )

df4 <- bind_rows(df3, species_df)

regimes_df <- df2 %>%
  select(Species, Regime)


df5 <- left_join(regimes_df, df4, by = "Species")


## save the species level data
write.csv(df5, "data/chameleons/chameleons_speciesdata.csv")


## some plotting
library(ggplot2)
library(ggrepel)

p <- ggplot(df5, aes(x = logbsl_mean, y = logpter_mean, color = factor(Regime))) +
  geom_point() +
  theme_classic() +
  geom_errorbar(aes(ymax = logpter_mean + logpter_varmean, ymin = logpter_mean - logpter_varmean)) +
  geom_errorbarh(aes(xmax = logbsl_mean + logbsl_varmean, xmin = logbsl_mean - logbsl_varmean)) +
  geom_text_repel(aes(label = Species)) +
  labs(
    y = "Pterygoid area (log (mm^2))",
    x = "Basal skull length (log (mm))",
  )

ggsave("figures/chameleons_pter_bsl.pdf", p, width = 200, height = 150, units = "mm")

























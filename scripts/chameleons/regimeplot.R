library(ape)
library(tibble)
library(tidytree)
library(ggtree)
library(ggplot2)
library(tidyverse)
library(tidyr)

phy <- read.nexus("data/chameleons/Tolley2013trimmedL.trees")

df <- as_tibble(read.csv("data/chameleons/chameleons_speciesdata.csv"))


#ans <- ace(df$Regime_label, phy, type = "discrete")
ans <- ace(df$Regime_label, phy, type = "discrete", marginal = TRUE)


internal_node_states <- colnames(ans2$lik.anc)[apply(ans$lik.anc, 1, which.max)]
tip_states <- df$Regime_label

node_states <- c(tip_states, internal_node_states)
edge_states <- node_states[phy$edge[,2]]

phy_df <- as_tibble(phy)

phy_df$regime <- node_states

td@data <- tibble(
  "node" = 1:39,
  #"edge" = 1:39,
  "regime" = node_states,
  "label" = c(phy$tip.label, rep("", 19))
)

td <- as.treedata(phy_df)

p1 <- ggtree(td, aes(color = regime)) + 
  geom_tree(linewidth = 1) +
  geom_tiplab(fontface = 3) +
  scale_color_manual(values = c("#ea078b", "black")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.1, 0.8)) +
  xlim(c(-0, 73))

#ggsave("figures/chameleon_ace_regimes.pdf", p, width = 200, height = 130, units = "mm")

## allometry figure
library(slouch)

m0 <- slouch.fit(
  phy,
  df$Species,
  response = df$logpter_mean,
  mv.response = df$logpter_varmean,
  direct.cov = df$logbsl_mean,
  mv.direct.cov = df$logbsl_varmean,
)
a <- m0$beta_primary$coefficients[1,1]
b <- m0$beta_primary$coefficients[2,1]

p2 <- ggplot(df, aes(x = logbsl_mean, y = logpter_mean, color = Regime_label)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c("#ea078b", "black")) +
  geom_abline(intercept = a, slope = b) +
  geom_errorbar(aes(ymin = logpter_mean - sqrt(logpter_varmean), ymax = logpter_mean + sqrt(logpter_varmean)), width = 0.02) +
  geom_errorbarh(aes(xmin = logbsl_mean - sqrt(logbsl_varmean), xmax = logbsl_mean + sqrt(logbsl_varmean)), height = 0.02) +
  labs(
    x = "log basal skull length (cm)",
    y = "log pterygoid area (cm^2)",
  )


## read the results from the RevBayes analysis
df1 <- read.table(paste0("output/chameleon_BM_state_dependent_run_1.log"), header = TRUE)
df2 <- read.table(paste0("output/chameleon_BM_state_dependent_run_2.log"), header = TRUE)
df3 <- read.table(paste0("output/chameleon_BM_state_dependent_run_3.log"), header = TRUE)
df4 <- read.table(paste0("output/chameleon_BM_state_dependent_run_4.log"), header = TRUE)

rb_logs <- bind_rows(df1, df2, df3, df4) %>%
  select(-starts_with("state_branch")) %>%
  select(-starts_with("num_changes")) %>%
  select(-starts_with("branch_rates"))

rb_logs$zeta_diff <- rb_logs$zeta.1. - rb_logs$zeta.2.
rb_logs$above <- rb_logs$zeta.1. > rb_logs$zeta.2.


zeta_plot <- rb_logs %>%
  ggplot(aes(x = zeta_diff, fill = above)) +
  geom_histogram() +
  theme_classic() +
  scale_fill_manual(values = c("black", "#ea078b")) +
  geom_vline(xintercept = 0) +
  labs(y = "posterior density",
       x = "zeta no ear - zeta ear")

p <- (p1 | p2 | zeta_plot) &
  theme(legend.position = "none")


ggsave("figures/chameleon_figure.pdf", p, width = 220, height = 80, units = "mm")




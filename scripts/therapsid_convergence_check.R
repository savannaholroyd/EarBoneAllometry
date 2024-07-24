library(tibble)
library(dplyr)

therapsid_log <- read.table("output/state_dependent_BM_character_1.log", sep = "\t", header = TRUE) |> 
  as_tibble() |>
  select(!starts_with("state_branch_rate")) |>
  select(!starts_with("branch_rates")) |>
  select(!starts_with("num_changes")) |>
  select(!starts_with("Posterior")) |>
  select(!starts_with("Likelihood")) |>
  select(!starts_with("Prior")) |>
  select(!starts_with("Iteration"))

therapsid_long <- pivot_longer(therapsid_log, !starts_with("Replicate_ID"))


therapsid_log0 <- dplyr::filter(therapsid_log, Replicate_ID == 0)
therapsid_log1 <- dplyr::filter(therapsid_log, Replicate_ID == 1)


cnames <- tail(names(therapsid_log0), n = -1)

D <- list()
pvalue <- list()
for (cname in cnames){
  ks <- ks.test(therapsid_log0[[cname]], therapsid_log1[[cname]])
  D[[cname]] <- ks$statistic[[1]]
  pvalue[[cname]] <- ks$p.value
}


df_convergence <- tibble(
  "parameter_name" = cnames,
  "D" = unname(unlist(D)),
  "p_value" = unname(unlist(pvalue))
)

## D critical = 0.09
## all parameters less than 0.09
print(df_convergence)



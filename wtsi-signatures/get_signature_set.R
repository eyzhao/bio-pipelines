' get_signature_set.R

Usage: get_signature_set.R -s SIGNATURES -m METRICS -o OUTPUT

Options:
    -s --signatures SIGNATURES      Path to cohort signatures
    -m --metrics METRICS            Path to cohort metrics for model selection
    -o --output OUTPUT              Path to output reference mutation set for cohort
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

signatures <- read_tsv(args[['signatures']])
metrics <- read_tsv(args[['metrics']])

auto_chosen_model = metrics %>%
  spread(metric, value) %>%
  mutate(
    stability_proportion = (stability - min(stability)) / (max(stability) - min(stability)),
    reconstruction_proportion = (reconstructionError - min(reconstructionError)) / (max(reconstructionError) - min(reconstructionError)),
    combined = stability_proportion - reconstruction_proportion
  ) %>%
  filter(combined == max(combined)) %>%
  filter(row_number() == n()) %>%
  .$n_signatures

signatures %>%
  filter(n_signatures == auto_chosen_model) %>%
  mutate(
    mutation_type = gsub('(.>.) (.)\\.(.)', '\\2[\\1]\\3', class),
    signature = factor(signature, levels = paste0('V', 1:length(unique(signature))))
  ) %>%
  select(
    mutation_type, signature, proportion
  ) %>%
  spread(signature, proportion) %>%
  write_tsv(args[['output']])

' mutation_timing_linear_fit.R

Usage: mutation_timing_linear_fit.R -i INPUT -o OUTPUT [ -r REF -s SIGNIT ]

Options:
    -i --input INPUT        Path to input signature timing summary
    -o --output OUTPUT      Path to output - summarized signature timing linear model fits
    -r --ref REF            Path to reference signatures. By default uses COSMIC 30 signatures
    -s --signit SIGNIT      Path to SignIT package, if it is not installed
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(parallel)
library(doParallel)

registerDoParallel(detectCores()-1)

if (is.null(args[['signit']])) {
    library(signit)
} else {
    devtools::load_all(args[['signit']])
}

if (is.null(args[['ref']])) {
    ref <- get_reference_signatures()
} else {
    ref <- read_tsv(args[['ref']])
}

signature_names <- ref %>% select(-mutation_type) %>% colnames

pog_signature_timing_raw <- read_tsv(args[['input']]) %>%
  mutate(
    patient = gsub('.*?patient\\/(.*?)\\/.*', '\\1', paths),
    sample = gsub('.*?sample\\/(.*?)\\/.*', '\\1', paths)
  ) %>%
  select(-paths)

pog_signature_timing <- pog_signature_timing_raw %>%
  group_by(
    patient,
    sample
  ) %>%
  filter(
    waic == min(waic)
  ) %>%
  distinct(
    patient,
    sample,
    population,
    population_proportion,
    prevalence_mean,
    n_mutations,
    n_populations,
    waic
  ) %>%
  tidyr::crossing(
    signature = signature_names
  ) %>%
  left_join(
    pog_signature_timing_raw %>% select(-bic),
    by = c(
      "population",
      "population_proportion",
      "prevalence_mean", 
      "n_mutations", 
      "n_populations",
      "waic",
      "patient",
      "sample",
      "signature"
    )
  ) %>%
  replace_na(list(
    mean = 0,
    median = 0,
    mode = 0,
    sd = 0
  ))

pog_signature_timing %>%
  filter(n_populations > 1) %>%
  group_by(
    patient, sample, signature
  ) %>%
  plyr::ddply(c('patient', 'sample', 'signature'), function(z) {
    message(sprintf('Patient: %s | Sample: %s | Signature: %s', z$patient[1], z$sample[1], z$signature[1]))
    model <- lm(
      mean ~ prevalence_mean,
      weights = population_proportion,
      data = z
    )
    tibble(
      slope = model$coefficients[['prevalence_mean']],
      intercept = model$coefficients[['(Intercept)']],
      min_prevalence = min(z$prevalence_mean),
      max_prevalence = max(z$prevalence_mean)
    ) %>%
      mutate(
        early_exposure = slope * max_prevalence + intercept,
        late_exposure = slope * min_prevalence + intercept
      )
  }, .parallel = TRUE) %>%
  write_tsv(args[['output']])

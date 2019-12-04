' compare_to_gtex.R

Use --tissuetypes to see a list of available tissue types.

Usage:
  compare_to_gtex.R -i INPUT -t TISSUE -s SAMPLE -r GTEX -o OUTPUT
  compare_to_gtex.R -s SAMPLE --tissuetypes

Options:
  -i --input INPUT      Path to input data (TPM values)
  -t --tissue TISSUE    Tissue comparator to use
  -s --sample SAMPLE    Path to GTEX sample attributes data.
  -r --gtex GTEX        Path to GTEX expression data (Expression TPM files)
  -o --output OUTPUT    Path to output comparison data
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(pander)


print_tissue_types <- function(sample_path) {
  print(sample_path)
  gtex_sample_data <- read_tsv(sample_path, guess_max=10000)
  gtex_sample_data %>%
    distinct(SMTS) %>%
    rename(`Available Tissue Types` = SMTS) %>%
    pander()
}

get_relevant_samples <- function(sample_path, tissue_type) {
  gtex_sample_data <- read_tsv(sample_path, guess_max=10000)
  tissue_data <- gtex_sample_data %>%
    filter(SMTS == tissue_type)

  if (nrow(tissue_data) == 0) {
    stop(sprintf('No samples found for tissue type: %s', tissue_type))
  }

  return(tissue_data$SAMPID)
}

get_expression_data <- function(expression_path, samples) {
  return(read_tsv(expression_path, skip=2, guess_max=10)[c('Name', 'Description', samples)])
}

get_expression_for_genes <- function(expression_data, genes) {
  expression_for_genes <- expression_data %>%
    filter(Description %in% genes) %>%
    gather(sample, tpm, -Name, -Description)

  print(expression_for_genes %>% as.data.frame())
}

main <- function(args) {
  if (args[['--tissuetypes']]) {
    suppressWarnings(print_tissue_types(args[['--sample']]))
    message("Tip: Use quotes around your tissue type when specifying it as an argument, as in 'Fallopian Tube'\n")
  } else {
    samples <- get_relevant_samples(args[['--sample']], args[['--tissue']])
    expression_data <- get_expression_data(args[['--gtex']], samples)
    get_expression_for_genes(expression_data, c('DDX11L1', 'WASH7P', 'FAM138A'))
  }
}

main(args)

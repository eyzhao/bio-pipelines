' run_signatureestimation.R - Wrapper for deconstructsigs to compute mutation signatures

Usage: deconstructsigs_script.R -i INPUT -o OUTPUT [ options ]

Options:
    -i --input INPUT    Path to mutation catalog (TSV file with columns mutation_type and count)
    -o --output OUTPUT  Output path to mutation signature exposures

    -x --exome          Flag to indicate the variants were called from exome, which will
                            trigger the script to apply trinucleotide frequency correction
                            and adjust for mutation

    -s --sigest SIGEST  Path to SignatureEstimation package, if not installed

    -n --signit SIGNIT  Path to SignIT package, if not installed

' -> doc

library(docopt)
args <- docopt(doc)
print(args)

library(devtools)
library(tidyverse)
library(deconstructSigs)
library(stringr)

if (! is.null(args[['sigest']])) {
    load_all(args[['sigest']])
} else {
    library(SignatureEstimation)
}

if (! is.null(args[['signit']])) {
    load_all(args[['signit']])
} else {
    library(signit)
}

catalog <- read_tsv(args[['input']])

if (args[['exome']]) {
    reference_signatures <- get_reference_signatures('cosmic_30_exome')
} else {
    reference_signatures <- get_reference_signatures('cosmic_30')
}

signature_names <- reference_signatures %>% select(-mutation_type) %>% colnames
n_mutations <- catalog$count %>% sum

exposures_output <- findSigExposures(
    M = catalog$count %>% as.matrix(),
    P = reference_signatures %>% reference_signatures_as_matrix(catalog),
    decomposition.method = decomposeSA
) %>%
    .$exposures

exposures_output %>%
    as.data.frame %>%
    rownames_to_column('signature') %>%
    rename(count = V1) %>%
    mutate(count = count * n_mutations) %>%
    write_tsv(args[['output']])

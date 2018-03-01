' run_mapsigs.R

Usage: run_mapsigs.R -i INPUT -o OUTPUT [ options ]

Options:
    -i --input INPUT            Table of SNVs, with columns chr, pos, ref, alt
    -o --output OUTPUT          Path to output RDS file
    -c --chrom CHROM            Chromosome name to filter for.
    -r --ref REF                Reference signatures to use (default: COSMIC 30 signature set)
    -m --mapsigs MAPSIGS        Path to MapSigs package (if not installed)
    -s --signit SIGNIT          Path to SignIT package (if not installed)

    --half-life HALFLIFE        Half-life of beta concentration decay factor (default: 200000).
                                    Higher values will result in more influence between neighbouring
                                    mutations, which assumes that mutational processes change less
                                    across the genome.

    --kappa-start KAPPA         The beta concentration between neighbouring mutations (default: 500).
                                    The higher the kappa-start, the less difference in mutation signatures
                                    between closely neighbouring mutations.

    --chains CHAINS             Number of chains to run (default: 3)
    --iter ITER                 Number of sampling iterations per chain (default: 100)
    --warmup WARMUP             Number of warmup iterations per chain (default: 1900)

    --recompile                 If present, recompiles the MapSigs code
' -> doc

library(docopt)
args <- docopt(doc)
print(args)

library(foreach)
library(tidyverse)
library(signit)
library(rstan)
library(loo)

if (is.null(args[['mapsigs']])) {
    library(mapsigs)
} else {
    devtools::load_all(args[['mapsigs']], recompile = args[['recompile']])
}

if (is.null(args[['signit']])) {
    library(signit)
} else {
    devtools::load_all(args[['signit']])
}

if (is.null(args[['ref']])) {
    reference_signatures = get_reference_signatures() %>% signit:::normalize_reference_signatures()
} else {
    reference_signatures = read_tsv(args[['ref']]) %>% signit:::normalize_reference_signatures()
}

mutations <- read_tsv(
        args[['input']],
        col_types = cols(
            chr = col_character(),
            pos = col_integer(),
            ref = col_character(),
            alt = col_character()
        )
    ) %>%
    filter(gsub('chr', '', chr) %in% as.character(1:22)) %>%
    mutate(mutation_type = get_snv_mutation_type(chr, pos, ref, alt))
print('mutations')
print(mutations)

start_time <- Sys.time()

mapsigs_output = MapSigs(
    mutations,
    chromosome = args[['chrom']],
    reference_signatures = reference_signatures,
    n_chains = ifelse(is.null(args[['chains']]), 3, as.numeric(args[['chains']])),
    n_iter = ifelse(is.null(args[['iter']]), 100, as.numeric(args[['iter']])),
    n_warmup = ifelse(is.null(args[['warmup']]), 1900, as.numeric(args[['warmup']])),
    half_life = ifelse(is.null(args[['half-life']]), 200000, as.numeric(args[['half-life']])),
    kappa_start = ifelse(is.null(args[['kappa-start']]), 500, as.numeric(args[['kappa-start']]))
)

end_time <- Sys.time()

mapsigs_output['run_time'] = as.double(end_time - start_time, units = 'secs')

mapsigs_output %>%
    saveRDS(args[['output']])

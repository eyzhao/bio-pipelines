' test_on_simulated_genomes.R

Usage: test_on_simulated_genomes.R -o OUTPUT -n NMUT -s NSIG -p PERTURB [ --signit SIGNIT --msimr MSIMR --sigest SIGEST ]

Options:
    -o --output OUTPUT              Path to output file (.tsv table)

    -n --nmut NMUT                  Number of mutations in simulated genome

    -s --nsig NSIG                  Number of active signatures to be randomly included in the model.

    -p --perturb PERTURB            Percentage perturbation of signatures. For example, entering
                                        10 means signatures will be perturbed randomly according to
                                        a normal probability distribution with an SD of 10% of the mean
                                        signature probability for each mutation type.

    --signit SIGNIT                 Path to SignIT package, if it is not installed.

    --msimr MSIMR                   Path to mSimR package, if it is not installed

    --sigest SIGEST                 Path to SignatureEstimation package, if it is not installed
' -> doc

library(docopt)
args <- docopt(doc)

library(devtools)
library(tidyverse)
library(rjags)
library(rstan)
library(nnls)
library(Rtsne)
library(doParallel)
library(deconstructSigs)
library(GenSA)

if (is.null(args[['signit']])) {
    library(signit)
} else {
    load_all(args[['signit']])
}

if (is.null(args[['msimr']])) {
    library(msimr)
} else {
    load_all(args[['msimr']])
}

if (is.null(args[['sigest']])) {
    library(SignatureEstimation)
} else {
    load_all(args[['sigest']])
}

results <- compare_signature_methods(
    c(
      SignIT_runner,
      deconstructSigs_runner,
      nnls_runner,
      SA_runner
    ), 
    simulation_args = list(
        perturbation_percent_deviation = args[['perturb']] %>% as.numeric,
        n_active_signatures = args[['nsig']] %>% as.numeric, 
        n_mutations = args[['nmut']] %>% as.numeric
    ), 
    full_output = FALSE,
    run_in_parallel = FALSE
)

print('Simulations complete.')

write_tsv(results, args[['output']])

print(paste0('Results saved to ', args[['output']]))

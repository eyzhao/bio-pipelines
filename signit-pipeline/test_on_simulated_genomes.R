' test_on_simulated_genomes.R

Usage: test_on_simulated_genomes.R -o OUTPUT -m METHOD -n NMUT -s NSIG -p PERTURB [ --nsim NSIM --signit SIGNIT --msimr MSIMR ]

Options:
    -o --output OUTPUT              Path to output file (.tsv table)

    -m --method METHOD              Can be one of the following:
                                        deconstructSigs
                                        SignIT
                                        nnls
    
    -n --nmut NMUT                  Number of mutations in simulated genome

    -s --nsig NSIG                  Number of active signatures to be randomly included in the model.

    -p --perturb PERTURB            Percentage perturbation of signatures. For example, entering
                                        10 means signatures will be perturbed randomly according to
                                        a normal probability distribution with an SD of 10% of the mean
                                        signature probability for each mutation type.

    --nsim NSIM                     Number of independent simulations to run (default is 1000)

    --signit SIGNIT                 Path to SignIT package, if it is not installed.

    --msimr MSIMR                   Path to mSimR package, if it is not installed
' -> doc

library(docopt)
args <- docopt(doc)

library(devtools)
library(tidyverse)
library(rjags)
library(nnls)
library(Rtsne)
library(doParallel)
library(deconstructSigs)

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

if (is.null(args[['nsim']])) {
    n_sim = 1000
} else {
    n_sim = args[['nsim']] %>% as.numeric
}

if (args[['method']] == 'SignIT') {
    run_function = SignIT_runner
} else if (args[['method']] == 'deconstructSigs') {
    run_function = deconstructSigs_runner
} else if (args[['method']] == 'nnls') {
    run_function = nnls_runner
} else {
    stop('Please provide a valid method.')
}

results <- test_signature_method(
    run_function, 
    n_simulations = n_sim, 
    simulation_args = list(
        perturbation_percent_deviation = args[['perturb']] %>% as.numeric,
        n_active_signatures = args[['nsig']] %>% as.numeric, 
        n_mutations = args[['nmut']] %>% as.numeric
    ), 
    full_output = FALSE,
    run_in_parallel = TRUE
)

print('Simulations complete.')

write_tsv(results, args[['output']])

print(paste0('Results saved to ', args[['output']]))

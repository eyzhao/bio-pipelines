' run_mutsignatures.R

Usage: run_mutsignatures.R -i INPUT -o OUTPUT -n NSIG [ -c NCORE -b NITER ]

Options:
    -i --input INPUT            Path to mutSignatures input file (.tsv with mutation types as rownames)

    -o --output OUTPUT          Path to mutSignatures output file (.Rds serialized object)

    -n --nsig NSIG              The number of signatures to run the model with

    -c --ncore NCORE            The number of cores to run with. If not specified, then will run with
                                the 3/4 of the maximum number of cores available by default.

    -b --niter NITER            The minimum target number of iterations to run with. The default value
                                is 1000. The number of iterations per core will be set to be
                                NITER / NCORE rounded up to the nearest integer.
' -> doc

library(docopt)
library(parallel)
args <- docopt(doc)

library(tidyverse)
library(mutSignatures)

input <- read.table(args[['input']], sep = '\t')
input_data <- setMutCountObject(input)

if (is.null(args[['ncore']])) {
    n_cores <- detectCores() - 1
} else {
    n_cores <- as.numeric(args[['ncore']])
}

if (is.null(args[['niter']])) {
    n_iter <- 1000
} else {
    n_iter <- as.numeric(args[['niter']])
}

n_iter_per_core <- ceiling(n_iter / n_cores)

message(
    sprintf(
        'Running analysis with a total of %s iterations across %s cores (%s iterations per core)', 
        n_cores * n_iter_per_core, 
        n_cores, 
        n_iter_per_core
    )
)

param <- setMutClusterParams(
    num.processes.toextract = as.numeric(args[['nsig']]),
    tot.iterations = n_iter_per_core,
    tot.cores = n_cores
)

output <- decipherMutationalProcesses(input_data, param)
saveRDS(output, args[['output']])

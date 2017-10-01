'
Usage: run_signeR -i INPUT -n NSIG -o OUTPUT

Options:
    -i --input INPUT            Path to input matrix. Input must be structured as described in
                                    https://bioconductor.org/packages/devel/bioc/vignettes/signeR/inst/doc/signeR-vignette.html.
                                    An example can be seen at 
                                    https://github.com/rvalieris/signeR/blob/devel/inst/extdata/21_breast_cancers.mutations.txt.
                                    The R object resulting from loading this can be seen in the interactive R terminal by running
                                    mut <- read.table(system.file("extdata","21_breast_cancers.mutations.txt",
                                        package="signeR"), header=TRUE, check.names=FALSE)

    -n --nsig NSIG              Number of signatures to evaluate

    -o --output OUTPUT          Path to output RDS file containing the signeR mutation signature object
' -> doc

library(docopt)
args <- docopt(doc)

library(signeR)

mut <- read.table(args[['input']], header=TRUE, check.names=FALSE)
signatures <- signeR(mut, 
                     nsig=as.numeric(args[['nsig']]))

saveRDS(signatures, args[['output']])


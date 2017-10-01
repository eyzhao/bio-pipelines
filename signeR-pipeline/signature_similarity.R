' signature_similarity.R

Determines the similarity of de novo mutation signatures to a defined
set of reference signatures.

Usage: signature_similarity.R -s SIGNATURES -r REFERENCE -o OUTPUT

Options:
    -s --signatures SIGNATURES  Path to mutation signature table. Tidy formatting with at minimum four columns named as follows.
                                    n_sigs: number of signatures used in the model
                                    signature_index: the name or numbering of each signature deciphered
                                    mutation_type: mutation types in the same format as reference signatures
                                    signature_proportion: the proportion of mutations of mutation_type in the specified signature
                                Example:
                                           n_sigs   mutation_type   signature_index signature_proportion
                                           <dbl>    <chr>           <chr>           <dbl>
                                    1      20       A[C>A]A         Signature 1     0.0159118426
                                    2      20       A[C>A]C         Signature 1     0.0043557273
                                    3      20       A[C>A]G         Signature 1     0.0033075889
                                    4      20       A[C>A]T         Signature 1     0.0067257024
                                    5      20       A[C>A]A         Signature 1     0.0062728748 

    -r --reference REFERENCE    Path to reference signature table. Download from:
                                http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt

    -o --output OUTPUT          Path to output table.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
signatures <- read_tsv(args[['signatures']])
reference <- read_tsv(args[['reference']])

reference_tidy <- reference %>%
    gather(
        reference_signature,
        reference_proportion,
        -`Substitution Type`, 
        -Trinucleotide, 
        -`Somatic Mutation Type`
    ) 

cosine_similarity <- function(x, y) {
    similarity = x %*% y / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
    return(as.numeric(similarity))
}

signatures %>%
    inner_join(
        reference_tidy,
        by = c('mutation_type' = 'Somatic Mutation Type')
    ) %>%
    group_by(
        n_sigs,
        signature_index,
        reference_signature
    ) %>%
    summarise(
        cosine_similarity = cosine_similarity(signature_proportion, reference_proportion)
    ) %>%
    ungroup() %>%
    write_tsv(args[['output']])

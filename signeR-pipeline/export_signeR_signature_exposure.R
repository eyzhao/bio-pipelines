'
Usage: export_signeR_signature_exposures.R -i INPUT -s SIGOUT -e EXPOUT

Options:
    -i --input INPUT            RDS file containing output object from the SigneR algorithm
    -s --signatures SIGOUT      Path to table with "tidy" mutation signatures
    -e --exposures EXPOUT       Path to table with "tidy" mutation signature exposures
' -> doc

library(docopt)
args <- docopt(doc)

library(signeR)
library(tidyverse)

print('Loading data')
signature_object <- readRDS(args[['input']])
mutations <- signature_object$SignExposures@mutations
samples <- signature_object$SignExposures@samples

print('Processing signatures')
signatures <- Average_sign(signature_object$SignExposures) %>% 
    as_tibble %>% 
    `colnames<-`(paste('Signature', 1:signature_object$Nsign)) %>%
    mutate(n_sigs = signature_object$Nsign,
           mutation_type = gsub('(.)>(.):(.)(.)(.)', '\\3[\\1>\\2]\\5', mutations)
    ) %>%
    gather(signature_index, signature_proportion, -n_sigs, -mutation_type)

print('Processing exposures')
exposures <- Average_exp(signature_object$SignExposures) %>%
    as_tibble %>%
    `colnames<-`(samples) %>%
    mutate(n_sigs = signature_object$Nsign,
           signature_index = paste('Signature', 1:signature_object$Nsign)
    ) %>%
    gather(sample, exposure, -n_sigs, -signature_index)

print('Writing signatures and exposures')
signatures %>% write_tsv(args[['signatures']])
exposures %>% write_tsv(args[['exposures']])

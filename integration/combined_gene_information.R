' Usage: combined_gene_information.R -d DRUGTARGET -g GERMLINE -s SNV -i INDEL -o OUTPUT
' -> doc

library(docopt)
args <- docopt(doc)

library(readr)

drug_target <- unique(read_tsv('meta/drugtarget_table_clean.txt', col_types = cols(
    seqnames = col_character(), copy_change_list = col_character())))
germline <- unique(read_tsv('meta/germline_subset.txt', col_types = cols(
    seqnames = col_character())))
snv <- unique(readRDS('rdata/somatic_mutations_filtered.RData'))
indel <- unique(readRDS('rdata/indels_filtered.RData'))


metadata <- read_tsv('meta/case_metadata.txt')

library(dplyr)

# compute germline zygosities

germline$patient <- sapply(strsplit(germline$tumour_id, '_'), function(z) {z[[1]]})
germline_drug_target <- merge(germline, drug_target, by.x = c('patient', 'gene'), by.y = c('patient', 'hugo'))

ERROR_RATE <- 0.01
germline_likelihoods <- cbind(dbinom(x = germline_drug_target$normal_var_depth, size = germline_drug_target$normal_depth, prob = ERROR_RATE),
                              dbinom(x = germline_drug_target$normal_var_depth, size = germline_drug_target$normal_depth, prob = 0.5),
                              dbinom(x = germline_drug_target$normal_var_depth, size = germline_drug_target$normal_depth, prob = 1 - ERROR_RATE))
germline_drug_target$germline_zygosity <- apply(germline_likelihoods, 1, function(z) {c(0, 0.5, 1)[which(z == max(z))]})

# bayesian model of somatic zygosities

n_t <- 

' Usage: subset_germline_calls.R -i INPUT -o OUTPUT -l GENESET
' -> doc

library(readr)
library(docopt)
args <- docopt(doc)

library(dplyr)

all_chr <- c(1:22, 'X', 'Y')

input_table <- read_tsv(args[['INPUT']])
input_table <- rename(input_table, start=Pos); input_table$end = input_table$start

gene_set <- read.table(args[['GENESET']],
                     header=F,
                     sep='\t',
                     col.names = c('chr', 'start', 'end', 'gene', 'ensembl'))

library(GenomicRanges)
input_gr <- GRanges(input_table); seqlevels(input_gr) <- all_chr
gene_set_gr <- GRanges(gene_set); seqlevels(gene_set_gr) <- all_chr

overlaps <- as.data.frame(findOverlaps(input_gr, gene_set_gr))
subset <- as.data.frame(input_gr[overlaps[, 1], ])

write_tsv(subset, args[['OUTPUT']])

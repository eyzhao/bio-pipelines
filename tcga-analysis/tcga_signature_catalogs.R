' tcga_signature_catalogs.R - Computes mutation signature catalogs for many TCGA MAFs

Usage: tcga_signature_catalogs.R -p PATHSFILE -o OUTPUT [ -c NCORES ]

Options:
    -p PATHSFILE    Path to file which contains the paths (one per line) of merged TCGA MAFs created by tcga_maf_merge.R
    -o OUTPUT       Path to output file, into which mutation catalogs will be deposited
    -c NCORES       Number of cores to use for parallel processing

The script outputs a large table. Each row represents a single sample within a single tumour type.
The first two columns are cancer_type, sample, which denote the tumour type and TCGA tumour ID.
The remaining 96 columns are the mutation counts for each mutation type.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(deconstructSigs)
library(doParallel)
library(BSgenome.Hsapiens.NCBI.GRCh38)


### Register cores for multicore processing

if (!is.null(args[['c']])) {
    registerDoParallel(cores = as.numeric(args[['c']]))
    run_in_parallel = TRUE
}


### Parse metadata for project file paths

metadata <- tibble(path = read_lines(args[['p']])) %>%
    separate(path,
             c('project',
               'cancer_type'),
             '\\.',
             remove = FALSE,
             extra = 'drop') %>%
    select(path, cancer_type)


### Import MAF files and compute catalogs
    
catalog <- plyr::ddply(metadata, 'cancer_type', function(cancertype_table) {
    cancertype <- cancertype_table[['cancer_type']][1]
    cat(sprintf('Computing catalogs for cancertype: %s\n', cancertype))
    suppressMessages(merged_maf <- read_tsv(cancertype_table[['path']][1], progress = FALSE))
    merged_maf <- merged_maf %>%
        mutate(var_caller_count = as.numeric(var_caller_count),
               Chromosome = gsub('chr', '', Chromosome)) %>%
        filter(var_caller_count >= 2) %>%
        as.data.frame()

    mutation_catalog <- mut.to.sigs.input(mut.ref = merged_maf,
                                          sample.id = 'Tumor_Sample_Barcode',
                                          chr = 'Chromosome',
                                          pos = 'Start_Position',
                                          ref = 'Reference_Allele',
                                          alt = 'Allele',
                                          bsg = BSgenome.Hsapiens.NCBI.GRCh38) %>%
        as.data.frame() %>%
        rownames_to_column('sample')

    return(mutation_catalog)
}, .parallel = run_in_parallel)

catalog %>% write_tsv(args[['o']])
print(sprintf('Output catalog to %s', args[['o']]))

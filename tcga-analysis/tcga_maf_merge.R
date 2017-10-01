' tcga_maf_merge.R - Merges TCGA MAF files and keeps track of variant caller hits

Usage: tcga_maf_merge.R -p PATHSFILE -o OUTPUTDIR [ -c NCORES ]

Options:
    -p PATHSFILE    Path to file which contains the paths of TCGA MAFs, one per line
    -o OUTPUTDIR    Path to output directory, in which merged MAFs will be deposited
    -c NCORES       Number of cores to use for parallel processing
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(doParallel)

### Register cores for multicore processing

if (!is.null(args[['c']])) {
    registerDoParallel(cores = as.numeric(args[['c']]))
    run_in_parallel = TRUE
}


### Parse metadata for project file paths

metadata <- tibble(path = read_lines(args[['p']])) %>%
    separate(path,
             c('project',
               'cancer_type',
               'variant_caller',
               'uuid',
               'mutation_type',
               'version',
               'version2',
               'filetype'),
             '\\.',
             remove = FALSE) %>%
    select(-version, -version2)


### Import and merge MAF files specified in the metadata
    
merged_mafs <- plyr::ddply(metadata, 'cancer_type', function(cancertype_table) {
    cancertype = cancertype_table[['cancer_type']][1]

    ### Note that TCGA uses multiple variant callers (current 4 caller types).
    ### For each cancer type, create list of MAFs, with one entry per variant caller type
    ### Create an extra column for each named after the variant caller, and set all values to 1

    maf_list <- plyr::dlply(cancertype_table, 'variant_caller', function(z) {
        variant_caller_name <- z[['variant_caller']][1]
        print(sprintf('Reading MAF: %s', z[['path']][1]))
        suppressMessages(maf <- read_tsv(z[['path']][1], comment = '#', progress = F) %>%
            select(Chromosome, Start_Position, End_Position, Strand, Reference_Allele, Allele, Tumor_Sample_Barcode))
        maf[[variant_caller_name]] = as.integer(1)
        return(maf)
    })

    ### Join all MAFs from different variant callers which have the same sample, chr, pos, alleles, and strand.
    ### Easiest way is to set merged_snv to the first entry of maf_list, then join the rest with a for loop.

    merged_snv <- maf_list[[1]]
    for (i in 2:length(maf_list)) {
        merged_snv <- merged_snv %>%
            full_join(maf_list[[i]],
                      by = c('Tumor_Sample_Barcode',
                             'Chromosome',
                             'Start_Position',
                             'End_Position', 'Reference_Allele', 
                             'Allele',
                             'Strand')) 
    }

    ### The columns named after each variant caller have values of either 1 or NA.
    ### Gather these columns and compute rowSums to get a vector which indicates how many
    ### variant callers found each mutation. Return the resulting table.

    merged_snv['var_caller_count'] = rowSums(merged_snv[, cancertype_table[['variant_caller']]], na.rm = T)
    merged_snv <- merged_snv %>%
        arrange(Chromosome, Start_Position, End_Position)

    write_tsv(merged_snv, sprintf('%s/TCGA.%s.somatic.merged.maf', args[['o']], cancertype))

    return(1)
}, .parallel = run_in_parallel)

'
Usage: tcga_maf_separate_samples.R -i INPUT -d OUTPUTDIR -s SAMPLECOLNAME [ -f FILTER -c CALLERCOUNTCOL -t CANCERTYPE ]

Options:
    -i INPUT            MAF file
    -d OUTPUTDIR        Output directory
    -s SAMPLECOLNAME    Column name of the sample ID, on which to separate into individual sample files

    -f FILTER           A number (minimum 1) indicating how many variant callers must agree to call the variant.
                            If this argument is included, CALLERCOUNTCOL must also be supplied.

    -c CALLERCOUNTCOL   The name of the column specifying how many variant callers agree on the variant.

    -t CANCERTYPE       The cancertype. If provided, this will be appended to the end of the filename
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

sample_column = args[['s']]
output_dir = args[['d']]

suppressMessages(maf <- read_tsv(args[['i']], comment = '#'))

if (!is.null(args[['f']])) {
    if (is.null(args[['c']])) {
        stop('If you provide -f, you must also provide -c')
    }
    min_caller_count = args[['f']] %>% as.numeric
    caller_count_col = args[['c']] %>% as.character

    maf[[caller_count_col]] <- as.numeric(maf[[caller_count_col]])
    maf <- maf[maf[[caller_count_col]] >= min_caller_count, ]
}

filtered <- maf %>% plyr::ddply(sample_column, function(sample_maf) {
    sample_name = sample_maf[[sample_column]]
    if (is.null(args[['t']])) {
        sample_maf %>% write_tsv(sprintf('%s/%s.maf', output_dir, sample_name[1]))
    } else {
        sample_maf %>% write_tsv(sprintf('%s/%s_%s.maf', output_dir, sample_name[1], args[['t']]))
    }
})

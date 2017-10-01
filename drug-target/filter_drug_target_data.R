' Usage: filter_drug_target_data.R -i PATHSFILE -o OUTPUT -l GENESET [ -c NCORES ]
' -> doc

library('docopt')
library('plyr')
library('dplyr')
library('tibble')
library('readr')
library('doParallel')
args <- docopt(doc)

paths <- data.frame(path = readLines(args[['PATHSFILE']]))

library('GenomicRanges')

if (!is.null(args[['GENESET']])) {
    geneset <- read_tsv(args[['GENESET']], col_names = c('chr', 'start', 'end', 'gene', 'ensg'))
    geneset_gr <- unique(GRanges(geneset))
    seqlevels(geneset_gr) <- c(1:22, 'X', 'Y')
}

if (!is.null(args[['NCORES']])) {
    print(sprintf('Opening %s cores...', args[['NCORES']]))
    registerDoParallel(cores = as.numeric(args[['NCORES']]))
    run_in_parallel = TRUE
}

output <- dlply(paths, 'path', function(z) {
    p <- as.character(z$path[1])
    print(p)

    print(sprintf('Reading file %s', p))
    data <- as_tibble(read_tsv(p, comment = '#', col_types = cols(
            chr = col_character(),
            copy_change_list = col_character()
        )
    ))
    print('Finished reading file')

    gr <- GRanges(data[!is.na(data$chr) & !is.na(data$start) & !is.na(data$end), ])

    print('Filtering using set of regions')
    overlaps <- findOverlaps(gr, geneset_gr)
    out <- as_tibble(gr[queryHits(overlaps), ])

    out$gene_name <- as.data.frame(geneset_gr)[subjectHits(overlaps), 'gene']

    fc_colname <- names(data)[grepl('FC_breast_Bodymap', names(data))][1]
    tcga_colname <- names(data)[grepl('BRCA_percentile', names(data))][1]
    norm_colname <- names(data)[grepl('BRCA_norm_percentile', names(data))][1]

    if (any(c(is.na(fc_colname), is.na(tcga_colname), is.na(norm_colname)))) {
        selected_cols = c('seqnames', 'start', 'end', 'hugo', 'copy_category', 'loh', 'copy_change_list', 'breakpoint',
                          'ploidy_corrected_copy_change', 'ploidy_corrected_copy_number','sv_genome', 'sv_transcriptome',
                          'gene_name')
        out <- out[, selected_cols]
        out$foldchange <- NA
        out$percentile <- NA
        out$normal_percentile <- NA
    } else {
        selected_cols = c('seqnames', 'start', 'end', 'hugo', 'copy_category', 'loh', 'copy_change_list', 'breakpoint',
                          'ploidy_corrected_copy_change', 'ploidy_corrected_copy_number', 'sv_genome', 'sv_transcriptome',
                          'gene_name', fc_colname, tcga_colname, norm_colname)

        out <- out[, selected_cols]
        out <- rename_(out, foldchange = fc_colname, percentile = tcga_colname, normal_percentile = norm_colname)
    }

    out$patient <- gsub('.*?(POG\\d\\d\\d).*', '\\1', p)
    out$prefix <- gsub('.*?([archbiop]+\\d).*', '\\1', p)
    out$file_path <- p

    return(out)
}, .parallel = run_in_parallel)

print('Combining tables and outputting')

output_table <- do.call('rbind', output[!vapply(output, is.null, logical(1))])
write_tsv(output_table, args[['OUTPUT']], col_names = T)

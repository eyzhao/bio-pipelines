' get_purdom_timing_input.R - Creates event timing input data for Purdom method

In order to use this code, the following R packages must be installed:
    docopt, dplyr, readr, tidyr, doParallel, devtools, stringr, cancerTiming, ggplot2

Also, the following must be installed from Bioconductor:
    VariantAnnotation, GenomicRanges

Usage: get_purdom_timing_input.R -l LOHPATH -s SNVPATH -m METADATA --hrdtools HRDTOOLS -o OUTPUT -d OUTPUTDIR [ -c NCORES ]

Options:
    -l LOHPATH          Path to TSV file, where column 1 is sample ID and column 2 is paths to LOH files
    -s SNVPATH          Path to SNV file, where column 1 is sample ID and column 2 is paths to VCF files
    -m METADATA         Path to aggregated metadata file, 
                            created using https://svn01.bcgsc.ca/svn/personal/ezhao/projects/signatures/trunk/pipelines/pog-files/aggregate_metadata.R
    --hrdtools HRD      Path to copy of the HRDtools code,
                            available at https://svn.bcgsc.ca/svn/personal/ezhao/projects/signatures/trunk/hrdtools
    -o OUTPUT           Path where complete output (containing all samples data) will be written to
    -d OUTPUTDIR        Path where timing input and timing output files for individual samples can be dumped to
    -c NCORES           Number of cores to run in parallel (leave out if no parallel processing desired)
' -> doc

#library(docopt)
#args <- docopt(doc)

library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(doParallel)
library(ggplot2)

library(devtools)
load_all(args[['--hrdtools']])

metadata <- read_tsv(args[['m']])
contamination_table <- metadata %>%
    dplyr::filter(!is.na(biofx_tc)) %>%
    dplyr::mutate(pog_id = str_extract(path, 'POG\\d\\d\\d')) %>%
    unite(sample_id, pog_id, sample_prefix) %>%
    dplyr::mutate(contamination = (100 - as.numeric(biofx_tc)) / 100) %>%
    dplyr::select(sample_id, contamination) %>%
    distinct()

paths <- read_tsv(args[['l']]) %>%
    dplyr::inner_join(read_tsv(args[['s']]), by = 'sample_id') %>%
    dplyr::inner_join(contamination_table, by = 'sample_id')

#' Carries out the event timing function from cancerTiming
#' 
#' @param inputTable Table input to cancerTiming
#' @param contamination Float between 0 and 1 representing the normal cell fraction
#' @importFrom cancerTiming eventTimingOverList
#' @export

event_timing <- function(inputTable, contamination, iterations = 10) {
    inputList <- list(inputTable)

    eventArgs = list()
    eventArgs['returnAssignments'] = TRUE

    if (iterations >= 3) {
        eventArgs['bootstrapCI'] = 'nonparametric'
        eventArgs['B'] = iterations
    }

    print('Running mutation timing algorithm...')
    timingObject <- eventTimingOverList(inputList, contamination, eventArgs)
    return(timingObject)
}

#' Assembles data for the CancerTiming method
#' 
#' @param loh_path Path to input LOH file, a TSV with at least columns for chr, start, end, copy number, and loh state
#' @param vcf_path Path to input VCF file
#' @param loh_colname String indicating the column name of the LOH type classification in file at loh_path
#' @param copynumber_colname String indicating name of the copy number column in file at loh_path
#' @param genomeVersion String to log the version of the genome
#' @importFrom hrdtools import_ranges
#' @importFrom hrdtools import_vcf
#' @importFrom GRanges seqlevels seqnames ranges findOverlaps
#' @importFrom plyr ddply
#' @export

get_timing_input <- function(loh_path, vcf_path, loh_colname, copynumber_colname, file_colnames = NULL, genomeVersion='hg19') {
    col_names <- c(loh_colname, copynumber_colname)

    print('Importing data...')
    print(loh_path)

    if (is.null(file_colnames)) {
        file_colnames <- c('chr', 'start', 'end', 'width', 'vars', 'copy_number', 'lohtype', 'allele1', 'allele2')
    }

    # Import data and sort
    loh_ranges <- sort(import_ranges(loh_path, col.names=file_colnames))
    snv_ranges <- sort(import_vcf(vcf_path, genomeVersion=genomeVersion))
    print(loh_ranges)

    loh_ranges <- contig(loh_ranges, loh_colname, copynumber_colname)    # Fill in seg gaps
    print(loh_ranges)
    loh_ranges <- filter_ranges(loh_ranges, loh_colname, copynumber_colname)          # Filter out very short segs
    loh_ranges$segId <- 1:length(loh_ranges)                    # Create a segId column

    merged <- get_merged_data(loh_ranges, snv_ranges, loh_colname, copynumber_colname)

    return(list('loh_ranges'=loh_ranges,
                'snv_ranges'=snv_ranges,
                'merged'=merged
                )
    )
    
}

get_merged_data <- function(loh_ranges, snv_ranges, loh_colname, copynumber_colname) {
    print('Computing overlaps...')

    overlaps <- as.data.frame(findOverlaps(loh_ranges, snv_ranges))

    print('Merging data...')

    loh_df = as.data.frame(loh_ranges)  # Convert to data frames
    snv_df = as.data.frame(snv_ranges)
    loh_idx <- overlaps$queryHits       # Retrieve overlapping indices
    snv_idx <- overlaps$subjectHits

    merged <- data.frame(
        mutationId = 1:length(loh_idx), 
        chr = as.character(snv_df[snv_idx, 'seqnames']), 
        pos = as.numeric(snv_df[snv_idx, 'start']),
        ref = snv_df[snv_idx, 'ref'],
        var = snv_df[snv_idx, 'alt'],
        snvType = paste(snv_df[snv_idx, 'ref'], snv_df[snv_idx, 'alt'], sep='>'),
        segId = loh_df[loh_idx, 'segId'],
        type = get_type_value(loh_df, loh_colname, copynumber_colname)[loh_idx],
        nMutAllele = snv_df[snv_idx, 'altDepth'],
        nReads = snv_df[snv_idx, 'totalDepth'],
        fraction = snv_df[snv_idx, 'altDepth'] / snv_df[snv_idx, 'totalDepth']
    )

    return(merged)
}

combine_timing_inputs <- function(timing_input_1, timing_input_2, loh_colname, copynumber_colname) {
    gr_1 <- timing_input_1$loh_ranges
    gr_2 <- timing_input_2$loh_ranges
    col_names <- c(loh_colname, copynumber_colname)

    unique_values <- as.data.frame(unique(values(c(timing_input_1$loh_ranges, timing_input_2$loh_ranges))[, col_names]))
    merged_loh_ranges <- apply(unique_values, 1, function(z) {
          subset_1 <- gr_1[values(gr_1)[, loh_colname] == z[1] & values(gr_1)[, copynumber_colname] == z[2]]
          subset_2 <- gr_2[values(gr_2)[, loh_colname] == z[1] & values(gr_2)[, copynumber_colname] == z[2]]
          int <- intersect(subset_1, subset_2)
          int <- attach_values(int, c(loh_colname, copynumber_colname), z)
          return(int)
    })
    merged_loh_ranges <- sort(do.call('c', merged_loh_ranges))
    merged_loh_ranges$segId <- 1:length(merged_loh_ranges)

    merged_1 <- get_merged_data(merged_loh_ranges, timing_input_1$snv_ranges, loh_colname, copynumber_colname)
    merged_2 <- get_merged_data(merged_loh_ranges, timing_input_2$snv_ranges, loh_colname, copynumber_colname)

    timing_input_1$loh_ranges <- merged_loh_ranges; timing_input_2$loh_ranges <- merged_loh_ranges
    timing_input_1$merged <- merged_1; timing_input_2$merged <- merged_2

    return(list('timing_input_1'=timing_input_1, 'timing_input_2'=timing_input_2))
}

get_type_value <- function(loh_df, loh_colname, copynumber_colname) {
    regionStart = loh_df[, 'start']
    regionEnd = loh_df[, 'end']
    lohtype = loh_df[, loh_colname]
    is_loh = lohtype %in% c('ALOH', 'DLOH', 'NLOH')
    copynumber = loh_df[, copynumber_colname]

    d <- data.frame(id=1:length(is_loh), is_loh, copynumber, lohtype)

    type <- ddply(d, 'id', function(z) {
        if (z$lohtype == 'NLOH') {
            t = "CNLOH"
        } else if (!z$is_loh && z$copynumber == 2) {
            t = "Diploid"
        } else if (!z$is_loh && z$copynumber == 3) {
            t = "SingleGain"
        } else if (!z$is_loh && z$copynumber == 4) {
            t = "DoubleGain"
        } else {
            t = "Other"
        }
        return(t) 
    }); colnames(type) <- c('id', 'type')
    return(type$type)
}

main <- function(loh_path, vcf_path, loh_colname, copynumber_colname, loh_file_col_names, contamination,
                 timing_input_destination='timing_input.RData', timing_output_destination='timing_output.RData') {
    library('VariantAnnotation')
    library('GenomicRanges')
    library('plyr')
    library('ggplot2')
    library('cancerTiming')

    col_names <- c(loh_colname, copynumber_colname)
    input <- get_timing_input(loh_path, vcf_path, loh_colname, copynumber_colname, loh_file_col_names)
    saveRDS(input, file = timing_input_destination)

    timing_object <- event_timing(input$merged, contamination)
    saveRDS(timing_object, file = timing_output_destination)

    return(timing_object)
}

if (! is.null(args[['c']])) {
    run_in_parallel = T
    print('Running in parallel')
    print(paste('Opening', args[['c']]))
    registerDoParallel(cores = as.numeric(args[['c']]))
} else {
    run_in_parallel = F
}

library(cancerTiming)
library(GenomicRanges)
library(VariantAnnotation)

complete_output <- plyr::dlply(paths, 'sample_id', function(z) {
    print(paste('Running for sample:', z$sample_id))
    z <- z[1, ]
    
    main(loh_path = z$loh_path, 
         vcf_path = z$snv_path,
         loh_colname = 'lohtype',
         copynumber_colname = 'copy_number',
         loh_file_col_names = c('chr', 'start', 'end', 'width', 'var_count', 'copy_number', 'lohtype', 'allele1', 'allele2'),
         contamination = z$contamination,
         timing_input_destination=paste0(args[['d']], '/', z$sample_id, '_metastasis_timing_input.RData'),
         timing_output_destination=paste0(args[['d']], '/', z$sample_id, '_metastasis_timing_output.RData')
    )

    print(paste('Completed sample:', z$sample_id))

}, .parallel = run_in_parallel)

saveRDS(complete_output, args[['o']])
print(paste('Complete output written to', args[['o']]))

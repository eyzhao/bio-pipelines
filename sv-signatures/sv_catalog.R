' Computes SV mutation catalogs, with counts for each SV type and subtype

Usage: calculate_sv_abyss_merged_catalogs.R -i INPUT -p POGID -s SAMPLEPREFIX -o OUTPUT

Options:
    -i INPUT            Input path (merged SV, .tsv) with at minimum the following columns indicating breakpoint coordinates:
                        chr1, pos1, chr2, pos2, sv_type
    
    -p POGID            String indicating the ID of the POG case being processed (i.e. POG303)

    -s SAMPLEPREFIX     Sample prefix (i.e. biop1)

    -o OUTPUT           Path to output file (.tsv) which will contain the resulting mutation catalog
' -> doc

label_clustered <- function(sv_table) {
    sv_table <- cbind(id=1:dim(sv_table)[1], sv_table)
    id_table <- cbind(sv_table[, c('id', 'chr1', 'pos1', 'chr2', 'pos2')])  # Use trans-ABySS breakpoint coordinates
    top <- id_table[, c(1,2,3)]; bottom <- id_table[, c(1,4,5)] 
    colnames(top) <- c('id', 'chr', 'pos'); colnames(bottom) <- c('id', 'chr', 'pos')
    breakpoints <- rbind(top, bottom)
    breakpoints$chr <- factor(as.character(breakpoints$chr), levels=c(as.character(1:22), 'X', 'Y'))
    breakpoints$pos <- as.numeric(breakpoints$pos)
    breakpoints <- breakpoints[order(as.numeric(breakpoints$chr), as.numeric(breakpoints$pos)), ]
    breakpoints <- plyr::ddply(breakpoints, 'chr', function(z) { 
        prev_pos <- c(0, z$pos[1:length(z$pos)-1])
        data.frame(id = z$id, chr = z$chr, pos = z$pos, distance = z$pos - prev_pos)
    })

    fit <- pcf(data = breakpoints[, c('chr', 'pos', 'distance')], kmin=10, gamma=25)
    mean_distance <- mean(breakpoints$distance)
    clustered_regions <- fit[fit$mean < mean_distance/10, c('chrom', 'start.pos', 'end.pos')]

    if (dim(clustered_regions)[1] == 0) {
        breakpoints$clustered <- FALSE
    } else {
        breakpoints$clustered <- apply(apply(clustered_regions, 1, function(z) {
            return(as.character(breakpoints$chr) == as.character(z['chrom'])
                   & as.numeric(breakpoints$pos) > as.numeric(z['start.pos'])
                   & as.numeric(breakpoints$pos) < as.numeric(z['end.pos'])) }),
        1, any)
    }

    breakpoints <- plyr::ddply(breakpoints[, c('id', 'clustered')], 'id', function(z) {data.frame(clustered=any(z$clustered))})

    merged <- merge(sv_table, breakpoints, by='id')

    return(merged)
}

calculate_sv_catalog <- function(sv_table, id, sample_prefix) {
    SV_LENGTH_BREAKS <- c(0, 1000, 10000, 100000, 1000000, 10000000, Inf)

    if (dim(sv_table)[1] > 0) {
        sv_table <- label_clustered(sv_table)

        if (! is.null(args[['MERGEOUT']])) {
            merge_path <- paste0(args[['MERGEOUT']], '/', id, '_', sample_prefix, '.sv-abyss-merge.txt')
            merge_path_clustered <- paste0(args[['MERGEOUT']], '/', id, '_', sample_prefix, '.sv-abyss-merge-clustered.txt')
            write_tsv(sv_table, merge_path)
            write_tsv(sv_table[sv_table$clustered, ], merge_path_clustered)
            print(paste0('Output dumped to ', merge_path))
        }

        catalog <- do.call('rbind', lapply(c(TRUE, FALSE), function(clustered_bool) {
            if (any(sv_table$clustered == clustered_bool)) {
                cluster_subtable <- sv_table[sv_table$clustered == clustered_bool, ]
                
                subcatalog <- do.call('rbind', lapply(c('deletion', 'duplication', 'inversion'), function(type) {
                    type_table <- cluster_subtable[cluster_subtable$sv_type == type, ]
                    if (dim(type_table)[1] == 0) {
                        counts <- rep(0, 6)
                    } else {
                        sv_length <- abs(as.numeric(type_table$pos2) - as.numeric(type_table$pos1))
                        counts <- hist(as.numeric(sv_length), SV_LENGTH_BREAKS, plot=FALSE)$counts
                    }

                    out <- data.frame(subtype = SV_LENGTH_BREAKS[1:6], value = counts)
                    out$sv_type <- type
                    return(out)
                }))
                subcatalog <- rbind(subcatalog, data.frame(sv_type = 'translocation', subtype = 0, value = c(sum(cluster_subtable$sv_type == 'translocation'))))
            } else {
                subcatalog <- data.frame(sv_type = c(rep('deletion', 6), rep('duplication', 6), rep('inversion', 6), 'translocation'),
                                            subtype = c(rep(SV_LENGTH_BREAKS[1:6], 3), 0),
                                            value = rep(0, 6+6+6+1))
            }
            subcatalog$clustered <- clustered_bool

            return(subcatalog)
        }))
        return(catalog) 
    } else {
        subtype = rep(c(rep(SV_LENGTH_BREAKS[1:6], 3), 0), 2)
        value = rep(0, 38)
        sv_type = rep(c(rep('deletion', 6), rep('duplication', 6), rep('inversion', 6), 'translocation'), 2)
        clustered = c(rep(TRUE, 19), rep(FALSE, 19))
        return(data.frame(subtype, value, sv_type, clustered))
    }
}

main <- function(input, pog_id, sample_prefix) {
    match <- read_tsv(input,
                      col_types = cols(chr1 = col_character(),
                                       chr2 = col_character(),
                                       pos1 = col_number(),
                                       pos2 = col_number())) %>%
        dplyr::select(chr1, pos1, chr2, pos2, sv_type)

    print('Computing SV catalog')

    sv_catalog <- calculate_sv_catalog(match, pog_id, sample_prefix) %>%
        filter(sv_type == 'translocation' | subtype > 0)

    return(sv_catalog)
}

load_libraries <- function() {
    library('tidyverse')
    library('copynumber')
    library('VariantAnnotation')
}

if (! interactive()) {
    library('docopt')
    args <- docopt(doc)

    load_libraries()

    input = args[['i']]
    pog_id = args[['p']]
    sample_prefix = args[['s']]

    sv_catalog <- main(input, pog_id, sample_prefix)

    write_tsv(sv_catalog, args[['o']])
}

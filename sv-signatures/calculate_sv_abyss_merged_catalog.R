' Usage: calculate_sv_abyss_merged_catalogs.R -p PATHSFILE -o OUTPUT [ --ncores NCORES --mergeout MERGEOUT ]
' -> doc

library('docopt')
args <- docopt(doc)

library('tidyr')
library('readr')
library('plyr')
library('copynumber')
library('VariantAnnotation')
library('doParallel')

label_clustered <- function(sv_table) {
    sv_table <- cbind(id=1:dim(sv_table)[1], sv_table)
    id_table <- cbind(sv_table[, c('id', 'chr1', 'pos1', 'chr2', 'pos2')])  # Use trans-ABySS breakpoint coordinates
    top <- id_table[, c(1,2,3)]; bottom <- id_table[, c(1,4,5)] 
    colnames(top) <- c('id', 'chr', 'pos'); colnames(bottom) <- c('id', 'chr', 'pos')
    breakpoints <- rbind(top, bottom)
    breakpoints$chr <- factor(as.character(breakpoints$chr), levels=c(as.character(1:22), 'X', 'Y'))
    breakpoints$pos <- as.numeric(breakpoints$pos)
    breakpoints <- breakpoints[order(as.numeric(breakpoints$chr), as.numeric(breakpoints$pos)), ]
    breakpoints <- ddply(breakpoints, 'chr', function(z) { 
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

    breakpoints <- ddply(breakpoints[, c('id', 'clustered')], 'id', function(z) {data.frame(clustered=any(z$clustered))})

    merged <- merge(sv_table, breakpoints, by='id')

    return(merged)
}

calculate_sv_catalog <- function(sv_table, p) {
    SV_LENGTH_BREAKS <- c(0, 1000, 10000, 100000, 1000000, 10000000, Inf)

    if (dim(sv_table)[1] > 0) {
        sv_table <- label_clustered(sv_table)

        if (! is.null(args[['MERGEOUT']])) {
            merge_path <- paste0(args[['MERGEOUT']], '/', p$id, '_', p$sample_prefix, '.sv-abyss-merge.txt')
            merge_path_clustered <- paste0(args[['MERGEOUT']], '/', p$id, '_', p$sample_prefix, '.sv-abyss-merge-clustered.txt')
            write_tsv(sv_table, merge_path)
            write_tsv(sv_table[sv_table$clustered, ], merge_path_clustered)
            print(paste0('Output dumped to ', merge_path))
        }

        catalog <- do.call('rbind', lapply(c(TRUE, FALSE), function(clustered_bool) {
            if (any(sv_table$clustered == clustered_bool)) {
                cluster_subtable <- sv_table[sv_table$clustered == clustered_bool, ]
                
                subcatalog <- do.call('rbind', lapply(c('DEL', 'DUP', 'INV'), function(type) {
                    type_table <- cluster_subtable[cluster_subtable$SVTYPE == type, ]
                    if (dim(type_table)[1] == 0) {
                        counts <- rep(0, 6)
                    } else {
                        sv_length <- abs(as.numeric(type_table$pos2) - as.numeric(type_table$pos1))
                        counts <- hist(as.numeric(sv_length), SV_LENGTH_BREAKS, plot=FALSE)$counts
                    }

                    out <- data.frame(subtype = SV_LENGTH_BREAKS[1:6], value = counts)
                    out$SVTYPE <- type
                    return(out)
                }))
                subcatalog <- rbind(subcatalog, data.frame(SVTYPE = 'TRA', subtype = 0, value = c(sum(cluster_subtable$SVTYPE == 'TRA'))))
            } else {
                subcatalog <- data.frame(SVTYPE = c(rep('DEL', 6), rep('DUP', 6), rep('INV', 6), 'TRA'),
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
        SVTYPE = rep(c(rep('DEL', 6), rep('DUP', 6), rep('INV', 6), 'TRA'), 2)
        clustered = c(rep(TRUE, 19), rep(FALSE, 19))
        return(data.frame(subtype, value, SVTYPE, clustered))
    }
}

path_to_catalog <- function(p) {
    delly_path <- p$delly_path[1]
    abyss_path <- p$abyss_path[1]
    
    print(delly_path)

    vcf <- readVcf(delly_path, 'hg19')
    table <- cbind(as.data.frame(rowRanges(vcf)), as.data.frame(info(vcf))) 
    somatic <- table[table$Somatic, ]
    
    abyss <- read_tsv(abyss_path)[, c('breakpoint')]
    abyss_tidy <- separate(data = abyss, col = breakpoint, into = c('chr1', 'pos1', 'chr2', 'pos2'), sep = '[\\|\\:]')

    match <- ddply(somatic, c('seqnames', 'start'), function(z) {
        z <- z[1, ]
        z$start <- as.numeric(z$start); z$end <- as.numeric(z$end);
        abyss_tidy$pos1 <- as.numeric(abyss_tidy$pos1); abyss_tidy$pos2 <- as.numeric(abyss_tidy$pos2);

        m <- abyss_tidy[((abyss_tidy$chr1 == z$seqnames & abyss_tidy$chr2 == z$CHR2 &
                          abs(z$start - abyss_tidy$pos1) < 20 & abs(z$END - abyss_tidy$pos2) < 20)
                        | (abyss_tidy$chr2 == z$seqnames & abyss_tidy$chr1 == z$CHR2 &
                          abs(z$END - abyss_tidy$pos1) < 20 & abs(z$start - abyss_tidy$pos2) < 20)), ]
        if(dim(m)[1] > 0) {
            return(cbind(z, m[1, ]))
        }
    })

    return(calculate_sv_catalog(match, p))
}

paths <- read_tsv(args[['PATHSFILE']])

if (! is.null(args[['NCORES']])) {
    run_in_parallel = T
    print('Running in parallel')
    print(paste('Opening', args[['NCORES']]))
    registerDoParallel(cores = as.numeric(args[['NCORES']]))
} else {
    run_in_parallel = F
}

print('Performing Calculations')

output_list <- dlply(paths, c('id', 'sample_prefix'), .parallel = run_in_parallel, .fun = path_to_catalog)

saveRDS(output_list, file = args[['OUTPUT']])

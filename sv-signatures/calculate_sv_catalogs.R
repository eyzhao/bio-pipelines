' Usage: calculate_sv_catalogs.R -p PATHSFILE -o OUTPUT
' -> doc

library('docopt')
args <- docopt(doc)

paths <- readLines(args[['PATHSFILE']])

library('plyr')
library('copynumber')
library('VariantAnnotation')

label_clustered <- function(sv_table) {
    sv_table <- cbind(id=1:dim(sv_table)[1], sv_table)
    id_table <- cbind(sv_table[, c('id', 'seqnames', 'start', 'CHR2', 'END')])
    top <- id_table[, c(1,2,3)]; bottom <- id_table[, c(1,4,5)] 
    colnames(top) <- c('id', 'chr', 'pos'); colnames(bottom) <- c('id', 'chr', 'pos')
    breakpoints <- rbind(top, bottom)
    breakpoints$chr <- factor(as.character(breakpoints$chr), levels=c(as.character(1:22), 'X', 'Y'))
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

calculate_sv_catalog <- function(sv_table) {
    SV_LENGTH_BREAKS <- c(0, 1000, 10000, 100000, 1000000, 10000000, Inf)
    sv_table <- label_clustered(sv_table)

    catalog <- do.call('rbind', lapply(c(TRUE, FALSE), function(clustered_bool) {
        if (any(sv_table$clustered == clustered_bool)) {
            cluster_subtable <- sv_table[sv_table$clustered == clustered_bool, ]
            
            subcatalog <- do.call('rbind', lapply(c('DEL', 'DUP', 'INV'), function(type) {
                type_table <- cluster_subtable[cluster_subtable$SVTYPE == type, ]
                if (dim(type_table)[1] == 0) {
                    counts <- rep(0, 6)
                } else {
                    counts <- hist(as.numeric(type_table$SVLEN), SV_LENGTH_BREAKS, plot=FALSE)$counts
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

        print(subcatalog)
        return(subcatalog)
    }))

    catalog <- catalog[catalog$subtype != 0 | catalog$SVTYPE == 'TRA', ]
    rownames(catalog) <- 1:(dim(catalog)[1])
    return(catalog)
}

output_list <- lapply(paths, function(p) {
    vcf <- readVcf(p, 'hg19')    
    table <- cbind(as.data.frame(rowRanges(vcf)), as.data.frame(info(vcf))) 
    somatic <- table[table$Somatic, ]
    calculate_sv_catalog(somatic)
})

names(output_list) <- paths

saveRDS(output_list, file = args[['OUTPUT']])

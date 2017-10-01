' Usage: count_cnv_overlaps -i CNVPATHSFILE -o OUTPUT
' -> doc

library('docopt')
args <- docopt(doc)

library('GenomicRanges')
library('plyr')
library('dplyr')
library('tidyr')
library('tibble')
library('readr')
library('reshape2')

paths <- readLines(args[['CNVPATHSFILE']])

complete_input <- lapply(paths, function(p) {
    table <- read_tsv(p, col_names = c('chr', 'start', 'end', 'state'))
    return(table[table$state != 2, ])
})
full_table <- do.call('rbind', complete_input)

recurrence_table_by_state <- ddply(full_table, 'state', function(state_table) {
    breakpoints <- as_tibble(melt(state_table[, 1:3], id.vars='chr'))
    breakpoints <- unique(breakpoints[, c('chr', 'value')])
    breakpoints <- arrange(breakpoints, chr, value)

    all_breakpoints <- GRanges(ddply(breakpoints, 'chr', function(z) {
        index <- 1:(dim(z)[1] - 1)
        start = z$value[index]
        end = z$value[index+1] - 1
        return(data.frame(start, end))
    }))

    gr <- GRanges(state_table)

    overlaps <- findOverlaps(gr, all_breakpoints)
    sh <- factor(subjectHits(overlaps), levels=1:length(all_breakpoints))
    all_breakpoints$count <- table(sh)

    return(as_tibble(all_breakpoints))
})

write.table(recurrence_table_by_state, file=args[['OUTPUT']], sep='\t', col.names=T, row.names=F, quote=F)

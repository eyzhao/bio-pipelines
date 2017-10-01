' Usage: plot_cnv_overlaps.R -i CNVTABLE -o PLOTPATH [ --wide ]
' -> doc

library('docopt')
args <- docopt(doc)

library('readr')
library('ggplot2')
library('dplyr')

print('Reading data')
cnv <- read_tsv(args[['CNVTABLE']])
cnv$seqnames <- factor(cnv$seqnames)
levels(cnv$seqnames) <- gsub('23', 'X', as.character(levels(cnv$seqnames)))
cnv$state <- factor(cnv$state); levels(cnv$state) <- c('HOMD', 'LOSS', 'GAIN', 'DOUBLEGAIN', 'AMPLIFICATION')

print('Generating plot object')
p <- ggplot(cnv, aes(xmin = start, xmax = end, ymin = 0, ymax = count.Freq)) + 
    geom_rect() +
    facet_grid(state ~ seqnames, space = 'free_x', scales = 'free_x')

if (args[['--wide']]) {
    print('Making WIDE version')
    p <- p + theme(panel.spacing.x = unit(0, 'lines'),
                   axis.text.x=element_text(angle = 45, hjust = 1))
    w = 30
} else {
    p <- p + theme(panel.spacing.x = unit(0, 'lines'),
                   axis.text.x=element_blank())
    w = 12
}

print('Outputting plot to PDF')
pdf(args[['PLOTPATH']], height = 6, width = w)
print(p)
dev.off()

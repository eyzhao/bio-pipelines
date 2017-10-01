' Usage: plot_sv_signatures.R -i INPUT -p PDFOUT
' -> doc

library(docopt)
args <- docopt(doc)

library(ggplot2)
library(reshape2)
library(readr)

signatures <- read_tsv(args[['INPUT']])
signatures$clustered <- c('Non-Clustered', 'Clustered')[as.numeric(grepl('clustered', signatures$type)) + 1]
signatures$sv_type <- sapply(strsplit(signatures$type, '_'), function(z) z[[length(z)]])
signatures <- signatures[, c('clustered', 'sv_type', colnames(signatures)[grepl('V\\d+', colnames(signatures))])]

signatures_melt <- melt(signatures, id.vars = c('clustered', 'sv_type'),
                        variable.name = 'Signature', value.name = 'Probability')

p <- ggplot(signatures_melt) +
  geom_bar(aes(x = sv_type, y = Probability), stat = 'identity') +
  facet_grid(Signature ~ clustered)

pdf(args[['PDFOUT']], height = 5, width = 5)
print(p)
dev.off()

' export_wtsi_metrics.R

Usage: export_wtsi_metrics.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT        Path to WTSI code output
    -o --output OUTPUT      Path to exported metrics
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(R.matlab)

wtsi <- readMat(args[['input']])

metrics <- data.frame(
    metric = c(
        'stability',
        'reconstructionError'
    ),
    value = c(
        mean(wtsi$processStabAvg),
        norm(wtsi$input[, , 1]$originalGenomes - wtsi$processes %*% wtsi$exposures, type = 'F')
    )
)

metrics %>%
    mutate(n_signatures = ncol(wtsi$processes)) %>%
    write_tsv(args[['output']])

' summarise_sciclone_output.R

Usage: summarise_sciclone_output.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT            SciClone output data object (RDS).
    -o --output OUTPUT          Output SciClone summary table (TSV).
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(sciClone)

sc <- readRDS(args[['input']])

if (is.null(sc)) {
    summary <- tribble(
        ~cluster, ~mean_vaf, ~proportion
    )
} else {
    summary <- tibble(
        cluster = 1:length(sc@clust$cluster.means),
        mean_vaf = as.numeric(sc@clust$cluster.means),
        proportion = sc@clust$pi
    )
}

write_tsv(summary, args[['output']])

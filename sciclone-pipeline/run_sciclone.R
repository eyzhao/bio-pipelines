' run_sciclone.R

Usage: run_sciclone.R -n NAME -v VAF -c CNV -e EXCLUDE -o OUTPUT

Options:
    -n --name NAME          Name of sample.
    -v --vaf VAF            Variant allele frequency input.
    -c --cnv CNV            Copy number variant input.
    -e --exclude EXCLUDE    Excluded regions (e.g. LOH) input.
    -o --output OUTPUT      Output location for the results of sciclone.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(sciClone)
library(tryCatchLog)

vaf <- read_tsv(
    args[['vaf']],
    col_types=cols(chr=col_character())
)

cnv <- read_tsv(
    args[['cnv']],
    col_types=cols(
        chr=col_character(),
        start=col_double(),
        stop=col_double()
    )
)

exclude <- read_tsv(
    args[['exclude']],
    col_types=cols(chr=col_character())
)

sc = NULL
tryCatch({
    sc = sciClone(
        vafs=vaf %>% as.data.frame,
        copyNumberCalls=cnv %>% as.data.frame,
        regionsToExclude=exclude %>% as.data.frame,
        sampleNames=args[['name']]
    )
}, error = function(e) {
    sc = NULL
})

saveRDS(sc, file=args[['output']])

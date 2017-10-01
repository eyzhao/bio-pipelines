' Usage: get_sv_paths.R -p GLOBSTRING -m METADATA -o OUTPUT
' -> doc

library('docopt')
library('readr')
library('tidyr')

args <- docopt(doc)

paths <- Sys.glob(args[['GLOBSTRING']])
paths_table <- data.frame(paths = basename(paths)) %>%
    separate(paths, into = c('ids'), sep='\\.') %>%
    separate(ids, c('pog', 'library1', 'library2'), sep='_')
paths_table <- cbind(data.frame(paths), paths_table)

metadata <- read_tsv(args[['METADATA']])
metadata <- metadata[grepl('biop', metadata$sample_prefix) | grepl('arch', metadata$sample_prefix), ]

merged1 <- merge(metadata, paths_table, by.y = 'library1', by.x = 'library_name')
output1 <- as.character(merged1[grepl('biop', merged1$sample_prefix) | grepl('arch', merged1$sample_prefix), 'paths'])
merged2 <- merge(metadata, paths_table, by.y = 'library2', by.x = 'library_name')
output2 <- as.character(merged2[grepl('biop', merged2$sample_prefix) | grepl('arch', merged2$sample_prefix), 'paths'])

output <- unique(sort(c(output1, output2)))

writeLines(output, con = args[['OUTPUT']])

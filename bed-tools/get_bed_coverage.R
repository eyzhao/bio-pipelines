' get_bed_coverage.R

Usage: get_bed_coverage.R -i INPUT -o OUTPUT

Options:
    -i INPUT        Globstring for all input files
    -o OUTPUT       Output for BED coverage table
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(doParallel)

registerDoParallel(cores = 30)

input_paths <- data.frame(path = Sys.glob(args[['i']]))

output <- plyr::ddply(input_paths, 'path', function(row) {
    print(row$path)
    input <- read_tsv(as.character(row$path), 
                      col_names = c('chr', 'start', 'end'),
                      col_types = cols(chr = col_character(),
                                       start = col_number(),
                                       end = col_number())) %>%
        filter(chr %in% c(as.character(1:22), 'X', 'Y')) %>%
        mutate(width = end - start)
    return(data.frame(path = row$path,
                      coverage = sum(input$width)))
}, .parallel = TRUE)

output %>% write_tsv(args[['o']])

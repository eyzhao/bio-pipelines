' get_purdom_summary.R

Usage: get_purdom_summary.R -t TIMINGOBJECT -s SEGFILE -o OUTPUT

Options:
    -t TIMINGOBJECT         Serialized R data object (.Rds file) output 
                            from cancerTiming analysis across multiple
                            events for a single sample

    -s SEGFILE              File containing table of segs, the events
                            which were timed. One event per row.
                            Must at minimum contain the columns: 
                                event_id (unique ID, must match TIMINGOBJECT seg IDs)
                                chr
                                start
                                end

    -o OUTPUT               Summary to be generated (.tsv file) for
                            the events contained in TIMINGOBJECT
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(cancerTiming)

timing_output <- readRDS(args[['t']])
segs <- read_tsv(args[['s']], col_types = cols(chr = col_character()))
summary <- getPi0Summary(timing_output) %>% 
    mutate(segId = segId %>% as.character %>% as.numeric) %>%
    inner_join(segs, by = c('segId' = 'event_id'))

summary %>% write_tsv(args[['o']])

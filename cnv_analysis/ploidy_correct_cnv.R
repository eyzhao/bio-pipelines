' ploidy_correct_cnv.R

Usage: ploidy_correct_cnv.R -i INPUT -o OUTPUT -p PLOIDY -t TC [-c COLS]

Options:
    -i --input INPUT        Path to input CNV file with columns chr, start, end, ratio, state
    -o --output OUTPUT      Path to output CNV segments with ploidy corrected copy numbers
    -p --ploidy PLOIDY      Ploidy as a number, i.e. diploid is 2, triploid is 3
    -t --tc TC              Tumour content as a percentage, i.e. 95 is a tumour content of 0.95

    -c --cols COLS          Comma separated column name values, if the CNV input file has unnamed
                                columns or differently named columns. Ensure that column names
                                chr, start, end, ratio, and state are among your custom colnames.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

if (! is.null(args[['cols']])) {
    data_colnames = strsplit(args[['cols']], ',')[[1]]
    data <- read_tsv(args[['input']], col_names = data_colnames)
} else {
    data <- read_tsv(args[['input']])
}

stopifnot(
    all(c('chr', 'start', 'end', 'ratio', 'state') %in% colnames(data))
)

tumour_content = as.numeric(args[['tc']]) / 100
ploidy = as.integer(args[['ploidy']])

data %>%
  filter(!is.na(chr)) %>%
  mutate(
    segment = cumsum(as.integer(state != c(0, state[1:(n()-1)])))
  ) %>%
  filter(! is.na(segment)) %>%
  group_by(chr, segment, state) %>%
  mutate(
    estimated_cn = ((1 + ratio) * ((tumour_content * ploidy) + ((1-tumour_content) * 2)) - (1-tumour_content) * 2) / tumour_content
  ) %>%
  summarise(
    mean_cn = sum((end - start) * estimated_cn) / sum(end - start),
    copy_number = round(mean_cn),
    start = start[1],
    end = end[n()]
  ) %>%
  ungroup() %>%
  mutate(
    state = factor(state, levels = state %>% unique %>% sort),
    copy_number = factor(copy_number, levels = copy_number %>% unique %>% sort)
  ) %>%
  write_tsv(args[['output']])

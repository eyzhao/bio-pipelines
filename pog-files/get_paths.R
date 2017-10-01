'
Usage: get_paths.R -G GLOB -r IDREGEX -n GROUPS -p PATHS -m MERGE

Options:
    -G --glob GLOB              Globstring which matches the files of interest

    -r --regex IDREGEX          Regex string which matches a unique ID per sample on which to merge.
                                    Groups will be captured. Multiple groups will be concatenated
                                    using an underscore by default or --concat-sep

    -n --group-names GROUPS     Names of groups, comma separated. I.e. patient,sample,run

    -p --path-name PATHS        Name of the column holding the paths. I.e. snv_paths

    -m --merge-col MERGE        Name of column to merge with STDIN. Must match one of the names
                                    in GROUPS.
' -> doc

suppressMessages(library(docopt))
suppressMessages(args <- docopt(doc))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))

group_names <- strsplit(args[['group-names']], ',')[[1]]
merge_cols <- strsplit(args[['merge-col']], ',')[[1]]
path_column_name <- args[['path-name']]

message(paste('Retrieving files for', path_column_name))

paths <- Sys.glob(args[['glob']])

message(paste(length(paths), 'files for', path_column_name))

path_table <- str_match(paths, args[['regex']]) %>% 
    as_tibble %>%
    select(-1) %>% 
    mutate(path = paths) %>%
    `colnames<-`(c(group_names, path_column_name)) %>%
    distinct()

path_table <- path_table[complete.cases(path_table), ]

message(paste(dim(path_table)[1], 'distinct files with valid IDs for', path_column_name))

if (! isatty(stdin())) {
    message(paste('Merging standard input with', path_column_name))
    input_table <- read_tsv(file('stdin'))
    path_table <- input_table %>% 
        inner_join(path_table, by = merge_cols)
}

cat(format_tsv(path_table))

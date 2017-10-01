' soft_link_delly.R

Usage: soft_link_delly.R -p DELLYPATHS -d OUTDIR -o OUTPUT

Options:
    -p DELLYPATHS       Points to file listing Delly paths, one per line
    -d OUTDIR           Prefix directory to output into
    -o OUTPUT           Where to dump the resulting bash commands
' -> doc

library(docopt)
args <- docopt(doc)
prefix <- args[['d']]

library(tidyverse)

delly <- read_tsv(args[['p']], col_names = 'delly_path')

output <- delly %>% plyr::ddply('delly_path', function(row) {
    path = row$delly_path[1]
    pog_id = gsub('.*?(POG\\d\\d\\d).*', '\\1', path)
    libraries = gsub('.*?delly\\/delly.*?\\/(.*?)\\/.*', '\\1', path)
    library_1 = strsplit(libraries, '_')[[1]][1]
    library_2 = strsplit(libraries, '_')[[1]][2]

    suppressMessages(
        flatfile <- read_tsv(sprintf('/projects/POG/POG_data/%s/flatfile/%s.tab', pog_id, pog_id)) %>%
            filter(library_name == library_1 | library_name == library_2,
                   ! startsWith(sample_prefix, 'blood'))
    )
    metadata <- flatfile[1, ]
    
    t_lib = metadata[['library_name']]
    n_lib = if_else(library_1 == t_lib, library_2, library_1)
    sample_prefix = metadata[['sample_prefix']]
    comparison = paste(sample_prefix, 't', t_lib, 'blood1', 'n', n_lib, sep='_')

    return(data.frame(pog_id, sample_prefix, t_lib, n_lib, comparison))
}) %>%
    filter(!is.na(pog_id) & !is.na(sample_prefix) & !is.na(t_lib) & !is.na(n_lib)) %>%
    mutate(command = sprintf('ln -sf %s %s/%s.%s.delly.vcf', delly_path, prefix, pog_id, comparison))

write_lines(output[['command']], args[['o']])

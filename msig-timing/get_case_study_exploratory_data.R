' get_case_study_exploratory_data.R

Usage: get_case_study_exploratory_data.R -p POGID -s SAMPLEPREFIX -o OUTPUT
' -> doc

library(docopt)

args <- docopt(doc)

library(dplyr)
library(readr)
library(tidyr)

flatfile_path <- sprintf('/projects/POG/POG_data/%s/flatfile/%s.tab', args[['POGID']], args[['POGID']])
flatfile_df <- read_tsv(flatfile_path)
libraries <- flatfile_df %>% 
    select(library_name, sample_prefix) %>%
    filter(sample_prefix == args[['SAMPLEPREFIX']]) %>%
    unique()

abyss_path <- Sys.glob(sprintf('/projects/POG/POG_data/%s/wgs/%s_*/trans-ABySS/fusions/LSR.tsv', args[['POGID']], args[['SAMPLEPREFIX']]))[1]
loh_path <- Sys.glob(sprintf('/projects/POG/POG_data/%s/wgs/%s_*_blood1_*/reviewed/loh/results/apolloh_out_segs.txt', args[['POGID']], args[['SAMPLEPREFIX']]))[1]

delly_paths <- sapply(libraries[['library_name']], function(z) {
    Sys.glob(sprintf('/projects/POG/POG_data/%s/sv/delly/delly-0.6.1/*%s*/Somatic_Germline_Quality_tagged.vcf',
             args[['POGID']], z))
})

delly_path <- as.character(delly_paths[file.exists(as.character(delly_paths))])[1]

abyss_df <- read_tsv(abyss_path)
loh_df <- read_tsv(loh_path, col_names = c('chr', 'start', 'end', 'width', 'nvar', 'cn', 'type', 'major_allele', 'minor_allele'))

library(VariantAnnotation)

delly_vr <- readVcf(delly_path, genome = 'hg19')

saveRDS(list(abyss = abyss_df, loh = loh_df, delly = delly_vr), args[['OUTPUT']])

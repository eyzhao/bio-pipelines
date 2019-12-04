' bionlp_annotation.R

Annotates OICR sample mutations with data from CancerMine or CIViCMine

Usage: bionlp_annotation.R -m MUTATIONS -c CLINICAL -n CNV -s SAMPLE -o OUTPUT -d DATA

Options:
  -m --mutations MUTATIONS        Path to cbioportal mutation data file
  -n --cnv CNV                    Path to cbioportal CNV data file
  -c --clinical CLINICAL          Path to cbioportal clinical samples data file
  -s --sample SAMPLE              ID of the sample for annotation
  -o --output OUTPUT              Path to output file
  -d --data DATA                  Path to cancermine/civicmine resource file, downloadable from 
                                    http://zenodo.org/record/3525385/files/cancermine_sentences.tsv
                                    http://zenodo.org/record/3529793/files/civicmine_sentences.tsv
' -> doc

library(docopt)
args <- docopt(doc)
print(args)

library(tidyverse)

data_clinical_samples <- read_tsv(
  args[['clinical']],
  comment='#'
)

cancertype <- data_clinical_samples %>%
  filter(SAMPLE_ID == args[['sample']]) %>%
  .$CANCER_TYPE_DETAILED %>%
  unique()

if (length(cancertype) > 1) {
  stop(sprintf('Nonunique cancer type for sample: %s', args[['sample']]))
}

data_mutations <- read_tsv(
  args[['mutations']],
  col_types = cols(
    Entrez_Gene_Id = col_character(),
    LEVEL_1 = col_character(),
    LEVEL_2A = col_character(), 
    LEVEL_2B = col_character(), 
    LEVEL_3A = col_character(), 
    LEVEL_3B = col_character(), 
    LEVEL_4 = col_character(), 
    LEVEL_R1 = col_character(),
    Highest_level = col_character()
  )
)

data_cnv <- read_tsv(
  args[['cnv']]
)
sample_cnv <- data_cnv[, c('Hugo_Symbol', args[['sample']])]
names(sample_cnv) <- c('Hugo_Symbol', 'cnv')
reported_cnv <- sample_cnv %>%
  filter(
    cnv >= 2 | cnv <= -2
  )

cancermine <- read_tsv(
  args[['data']],
  col_types = cols(
    gene_entrez_id = col_character(),
    title = col_character(),
    pmid = col_character(),
    year = col_integer(),
    month = col_character(),
    day = col_character()
  )
)

doid_mapping <- read_tsv('manual/octane_cancertype_mapping.tsv')
doid_values = doid_mapping %>%
  filter(CANCER_TYPE_DETAILED == cancertype) %>%
  .$Disease_Ontology %>%
  unique() %>%
  strsplit('\\|') %>%
  .[[1]]

mutation_matches <- data_mutations %>% 
  mutate(
    mutation_report = sprintf('%s: %s', BIOTYPE, HGVSp_Short)
  ) %>%
  filter(Tumor_Sample_Barcode == args[['sample']]) %>% 
  group_by(Hugo_Symbol, Entrez_Gene_Id) %>%
  summarise(mutations = paste(mutation_report, collapse=' \\| ')) %>%
  ungroup() %>%
  crossing(doid = doid_values) %>%
  left_join(
    cancermine,
    by = c(
      'Entrez_Gene_Id' = 'gene_entrez_id',
      'doid' = 'cancer_id'
    )
  ) %>%
    filter(! is.na(matching_id)) %>%
  mutate(
    url = if_else(
      !is.na(pmid),
      paste0('www.ncbi.nlm.nih.gov/pubmed/', pmid),
      'No annotation available'
    )
  ) %>%
  select(-gene_normalized)

message(sprintf('%s reported CNVs', nrow(reported_cnv)))

cnv_matches <- reported_cnv %>%
  mutate(mutation_report = sprintf('CNV: %s', cnv)) %>%
  group_by(Hugo_Symbol) %>%
  summarise(mutations = paste(mutation_report, collapse=' \\| ')) %>%
  ungroup() %>%
  crossing(doid = doid_values) %>%
  left_join(
    cancermine,
    by = c(
      'Hugo_Symbol' = 'gene_normalized',
      'doid' = 'cancer_id'
    )
  ) %>%
  filter(! is.na(matching_id)) %>%
  mutate(
    url = if_else(
      !is.na(pmid),
      paste0('www.ncbi.nlm.nih.gov/pubmed/', pmid),
      'No annotation available'
    )
  ) %>%
  rename(Entrez_Gene_Id = gene_entrez_id)

message(sprintf('%s matched CNVs', nrow(cnv_matches)))

matches <- rbind(mutation_matches, cnv_matches)

if (grepl('civicmine', args[['data']])) {
  matches %>%
    select(
      Hugo_Symbol,
      mutations,
      cancer_normalized,
      evidencetype,
      url,
      year,
      month,
      day,
      everything()
    ) %>%
    select(-sentence) %>%
    write_tsv(args[['output']])
} else if (grepl('cancermine', args[['data']])) {
  matches %>%
    select(
      Hugo_Symbol,
      mutations,
      cancer_normalized,
      role,
      url,
      year,
      month,
      day,
      everything()
    ) %>%
    select(-sentence) %>%
  write_tsv(args[['output']])
}

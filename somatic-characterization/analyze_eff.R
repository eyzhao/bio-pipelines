' analyze_eff -v VROBJ -o OUTPUT
' -> doc

library('docopt')
args <- docopt(doc)

library('tibble')
library('tidyr')
library('plyr')
library('dplyr')

vr <- readRDS(args[['VROBJ']])

select_eff <- function(eff) {
    eff_df = data.frame(eff);
    eff_df <- separate(eff_df, eff, c('Effect', 'Annotation'), sep='[\\(\\)]', extra='drop');
    eff_df <- separate(eff_df, Annotation, sep='\\|', c('effect_impact', 'functional_class', 'codon_change', 'amino_acid_change', 'amino_acid_length', 'gene_name', 'transcript_biotype', 'gene_coding', 'transcript_id', 'exon_rank', 'genotype_number'), extra='warn');
    return(eff_df)
}

vr$type <- sapply(vr$eff, function(z) {
    eff <- select_eff(z)
    if (any(eff$functional_class == 'NONSENSE')) {
        return('NONSENSE')
    } else if (any(eff$functional_class == 'MISSENSE')) {
        return('MISSENSE')
    } else if (any(eff$Effect == 'SPLICE_SITE_DONOR' | eff$Effect == 'SPLICE_SITE_ACCEPTOR')) {
        return('SPLICE')
    } else {
        return(NA)
    } 
})

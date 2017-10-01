' deconstructsigs_script.R - Wrapper for deconstructsigs to compute mutation signatures

Usage: deconstructsigs_script.R -p PATHSFILE -o OUTPUT

Options:
    -p PATHSFILE        File containing paths to SNV calls as VCF files, one path per line
    -o OUTPUT           Path to output, which will contain the deconstructsigs matrix

' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(deconstructSigs)
library(VariantAnnotation)
library(doParallel)
library(stringr)

paths <- data.frame(snv = read_lines(args[['p']]))

registerDoParallel(cores = 30)

exposures <- paths %>%
    plyr::ddply('snv', function(p) {
          path <- as.character(p$snv[1])
          pog_id <- str_extract(path, 'POG\\d\\d\\d')
          sample_prefix <- str_extract(path, '[bioparch]+\\d')

          cat(paste('Current sample:', pog_id, sample_prefix, '\n'))

          cat(paste('Reading VCF for', pog_id, sample_prefix, '\n'))
          if (file.exists(path)) {
              vr <- readVcf(path)
          } else {
              cat('Warning: FILE DOES NOT EXIST', paste(as.character(path)), '\n')
          }

          mutations <- rowRanges(vr) %>%
            as.data.frame() %>%
            dplyr::filter(elementNROWS(CharacterList(ALT)) == 1) %>%
            dplyr::mutate(sample_name = str_extract(path, 'POG\\d\\d\\d'),
                          ALT = unlist(CharacterList(ALT)))

          cat(paste('Computing Catalogs', pog_id, sample_prefix, '\n'))
          catalog <- mut.to.sigs.input(mut.ref = mutations,
                                        sample.id = 'sample_name',
                                        chr = 'seqnames',
                                        pos = 'start',
                                        ref = 'REF',
                                        alt = 'ALT')

          cat(paste('Computing exposures', pog_id, sample_prefix, '\n'))
          e_vector <- whichSignatures(catalog,
                                     contexts.needed = TRUE,
                                     signatures.ref = signatures.cosmic)
          return(e_vector$weights)
    }, .parallel = TRUE)

write_tsv(exposures, args[['o']])
print(paste('Wrote output to', args[['o']]))

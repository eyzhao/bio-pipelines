' Usage: pair_sv_abyss_paths.R -s SVPATHFILE -a ABYSSPATHFILE -m METADATAFILE -o OUTPUT
' -> doc

library(docopt)
args <- docopt(doc)

library(readr)

sv_paths <- readLines(args[['SVPATHFILE']])
abyss_paths <- readLines(args[['ABYSSPATHFILE']])
metadata <- read_tsv(args[['METADATAFILE']])

sv <- data.frame(id = gsub('.*?(POG\\d\\d\\d).*', '\\1', sv_paths),
                 library1 = gsub('.*?POG\\d+\\_(.*?)\\_.*', '\\1', sv_paths),
                 library2 = gsub('.*?POG\\d+\\_.*?\\_(.*?)\\..*', '\\1', sv_paths),
                 delly_path = sv_paths)
print(head(sv))

abyss <- data.frame(id = gsub('.*?(POG\\d\\d\\d).*', '\\1', abyss_paths),
                    sample_prefix = gsub('.*?POG\\d+\\_(.*?)\\..*', '\\1', abyss_paths),
                    abyss_path = abyss_paths)
print(head(abyss))

sv_metadata_1 <- merge(sv[, c('id', 'library1', 'delly_path')],
                       metadata[, c('library_name', 'sample_prefix')],
                       by.x = 'library1', by.y='library_name')

sv_metadata_2 <- merge(sv[, c('id', 'library2', 'delly_path')],
                       metadata[, c('library_name', 'sample_prefix')], 
                       by.x = 'library2', by.y='library_name')

colnames(sv_metadata_1)[1] <- 'library'
colnames(sv_metadata_2)[1] <- 'library'
sv_metadata <- rbind(sv_metadata_1, sv_metadata_2)

sv_metadata_abyss <- merge(sv_metadata, abyss, by = c('id', 'sample_prefix'))

write_tsv(sv_metadata_abyss, args[['OUTPUT']])

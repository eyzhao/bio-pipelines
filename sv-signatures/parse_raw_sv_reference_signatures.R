' Usage: parse_raw_sv_reference_signatures.R -i REFERENCE -o OUTPUT
' -> doc

library(docopt)
args <- docopt(doc)

library(readr)

reference_raw <- read_tsv(args[['REFERENCE']])

is_clustered <- reference_raw$Type == 'clustered'
split <- strsplit(reference_raw$Size, ':')
type_raw <- toupper(sapply(split, function(z) { return(z[[1]]) }))
type_sv <- gsub('TRANS', 'TRA', gsub('TDS', 'DUP', type_raw))
clustered_status <- rep('', length(is_clustered)); clustered_status[is_clustered] = 'clustered_'

length_string <- sapply(split, function(z) { if (length(z) == 2) {return(z[[2]])} else {return('tra')} })
length_string[length_string == '1-10kb'] <- '1 kb'
length_string[length_string == '10-100kb'] <- '10 kb'
length_string[length_string == '100kb-1Mb'] <- '100 kb'
length_string[length_string == '1Mb-10Mb'] <- '1 Mb'
length_string[length_string == '>10Mb'] <- '10Mb'

type <- paste0(clustered_status, type_sv)
subtype <- length_string
class <- paste(type, subtype)
values <- apply(reference_raw[, c(3:dim(reference_raw)[2])], 2, function(z) { num <- as.numeric(gsub('\\%', '', z)); return(num/sum(num)) })

type_df <- data.frame(type, subtype, class)
out_df <- cbind(type_df, values)

write_tsv(out_df, path = args[['OUTPUT']])

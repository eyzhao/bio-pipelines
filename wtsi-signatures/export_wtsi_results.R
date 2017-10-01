' Usage: export_wtsi_results.R -c CONFIG -e EXPOSURES -s SIGNATURES
' -> doc

library('R.matlab')
library('docopt')

args <- docopt(doc)

file_path <- readLines(args[['CONFIG']], n=1)

output <- readMat(file_path)

sample_names <- unlist(output$input[, , 1]$sampleNames)
exposures <- as.data.frame(t(output$exposures))
exposures <- cbind(sample_names, exposures)

write.table(exposures, file=args[['EXPOSURES']], col.names=T, row.names=F, sep='\t', quote=F)

type <- unlist(output$input[, , 1]$types)
subtype <- unlist(output$input[, , 1]$subtypes)
class <- paste(type, subtype)

signatures <- cbind(data.frame(type, subtype, class), as.data.frame(output$processes))

write.table(signatures, file=args[['SIGNATURES']], col.names=T, row.names=F, sep='\t', quote=F)

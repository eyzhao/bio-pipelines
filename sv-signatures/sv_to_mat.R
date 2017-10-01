' Usage: sv_to_mat.R -i INPUT -o OUTPUT -c CANCERTYPE
' -> doc

library(docopt)
library(plyr)
library('R.matlab')

args <- docopt(doc)

input_path <- args[['INPUT']]
output_path <- args[['OUTPUT']]

cancerType <- args[['CANCERTYPE']]

sv <- readRDS(input_path)

sv <- llply(sv, function(z) {
    z[z$SVTYPE == 'TRA' | z$subtype != 0, ]
})

sampleNames <- names(sv)

example_sv <- sv[[1]]
svtype <- example_sv$SVTYPE
clustered <- c('', 'clustered_')[as.numeric(example_sv$clustered) + 1]
types <- paste0(clustered, svtype)

convert_subtype <- data.frame(numeric = c(0, 1000, 10000, 100000, 1000000, 10000000),
                              string = c('tra', '1 kb', '10 kb', '100 kb', '1 Mb', '10Mb'))
subtypes <- as.character(sapply(example_sv$subtype, function(z) {
    convert_subtype[which(convert_subtype$numeric == z), 2]
}))

originalGenomes <- sapply(sv, function(z) { 
    as.double(z[order(-z$clustered, z$SVTYPE, z$subtype), 'value'])
})
originalGenomes <- as.matrix(originalGenomes)
rownames(originalGenomes) <- NULL

print(sampleNames)
print(types)
print(subtypes)
print(dim(originalGenomes))

print('writing output to matlab file')

writeMat(con=output_path,
         cancerType = cancerType,
         originalGenomes = originalGenomes,
         sampleNames = sampleNames,
         subtypes = subtypes,
         types = types)

#!/usr/bin/env Rscript
library('Cairo')
library('RSvgDevice')
library('RColorBrewer')
suppressWarnings(library(optparse))
options(echo=FALSE) # if you want see commands in output file

prepare_optparser <- function(){
 #Prepare optparser object. New options will be added in this function first.

 option_list <- list(
     make_option(c("-i","--input"), type="character", help="Mutation signature exposures file", dest="input_file"),
     make_option(c("-d", "--data"), dest="data_file", type="character", help="Path to file containing exposure data for all POGs", default="/projects/ezhao_prj/signatures/data/POG/figures/nnls_exposures.txt"),
     make_option(c("-o", "--output"), dest="output_file", type="character", help="Output destination to the plot", default=FALSE),
     make_option(c("-a", "--add"), dest="add_to_data", action="store_true", help="Add this case to the data store", default=FALSE),
     make_option(c("--allow-duplicate-row"), dest="allow_duplicate_row", action="store_true", help="Allow duplicate signature data into the data store", default=FALSE)
     )

 optparser <- OptionParser(option_list=option_list,
                           usage="usage: %prog -i <input_file> -d <reference data file> -o <output_file> [ -a ]",
                           description="Generate mutation signature comparison boxplot.")

 return(optparser)
}

optparser <- prepare_optparser()
args <- parse_args(optparser)
print(args$input_file)

# Read in the reference data file
exposureTable = read.table(args$data_file, sep='\t', header=TRUE, stringsAsFactors=FALSE)

headerNames = names(exposureTable)
signames = headerNames[2 : length(headerNames)]
#signames = c("1A", "1B", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "R1", "R2", "R3", "U1", "U2")

# Process the reference data file, turning exposure values into proportions
if (dim(exposureTable)[1] > 0) {
    data <- as.matrix(exposureTable[1:dim(exposureTable)[1], 2:dim(exposureTable)[2]])
    normalizedData = data
    for (i in 1:dim(normalizedData)[1]) {
        normalizedData[i, ] = normalizedData[i, ] / sum(normalizedData[i, ])
    }
}

input = read.table(args$input_file, header=TRUE, sep='\t', stringsAsFactors=FALSE)
queryVector = t(as.matrix(input[, 2]))

# double check that the data signature names match with the ones in the patient's data
sampleSignatureNames = make.names(input[, 1])
if (identical(signames, sampleSignatureNames)) {
    print('Data signature names match sample signature names.')
} else {
    print(signames)
    print(sampleSignatureNames)
    stop('Signature names do not match between data file and input file.')
}

# make figure
if (! args$output_file == FALSE) {
    joint = rbind(normalizedData, queryVector)
    percColPalette <- colorRampPalette(brewer.pal(9, 'YlOrRd'))(10)
    percVector = c()
    percColorVector = c()
    for (i in 1:dim(normalizedData)[2]) {
        perc = (trunc(rank(joint[, i])) / length(joint[, i])) [dim(joint)[1]]
        percVector = rbind(percVector, perc)
        percColor = percColPalette[ceiling(perc * 10)]
        percColorVector = rbind(percColorVector, percColor)
    }

    outputFileType = substr(args$output_file, sapply(gregexpr("\\.", args$output_file), tail, 1), nchar(args$output_file))
    print(outputFileType)

    if (outputFileType == '.pdf') {
        CairoPDF(file = args$output_file, width = 8, height = 8) # Output
    } else if (outputFileType == '.png') {
        CairoPNG(file = args$output_file, width = 600, height = 600, onefile = TRUE, bg = 'white')
    } else if (outputFileType == '.svg') {
        print('exporting SVG')
        devSVG(file = args$output_file, width = 5, height = 5,       
     onefile = TRUE)
    } else {
        print('Output file types supported: PDF, SVG, PNG')
        stop()
    }

    boxplot(normalizedData, names=signames, las=2, border = 'White', col='Grey', whisklty='solid', whiskcol='Grey', outcol=rgb(.5, .5, .5, .3), pch=20, xlab='Signature', ylab='Mutation fraction contributed')
    rect(-.25+(1:dim(normalizedData)[2]), queryVector-.009, 0.25 + (1:dim(normalizedData)[2]), queryVector+.009, col=percColorVector, border=NA)
    legend(24, 1.02, legend = paste(rev(seq(0, 90, 10)), rev(seq(10, 100, 10)), sep=' - '), col=rev(percColPalette), pch=15, title='Percentile', bty='n')

    dev.off()
}

# add patient's data to the data file

if (args$add_to_data) {
    duplicateExists = FALSE
    if (dim(exposureTable)[1] > 0) {
        for (i in 1:dim(data)[1]) {
            if (identical(queryVector, data[i, ])) {
                print(queryVector)
                print(data[i, ])
                print(paste('Duplicate row with name:', input[i, 1]))
                duplicateExists = TRUE
            }
        }
    }

    if (basename(args$input_file) %in% row.names(exposureTable)) {
        stop(paste('Sample already exists:', basename(args$input_file)))
    } else if (duplicateExists) {
        stop(paste('Duplicate row exists in data under a different name. Use cmd line argument --allow-duplicate-row to override.'))
    } else {
        writeVector = c(basename(args$input_file), input[, 2])
        write(writeVector, file=args$data_file, ncolumns=length(writeVector), append=TRUE, sep='\t')
    }
}

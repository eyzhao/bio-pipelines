#!/usr/bin/env Rscript

' obtain_timing_input.R

Usage: obtain_timing_input.R -p PATHFILE

Options:
    -p --paths PATHFILE     File containing paths of data
' -> doc

library('docopt')

args = docopt(doc)

library('hrdtimR')
library('VariantAnnotation')
library('plyr')
library('parallel')

open_cores <- function(numberOfIterations) {
    maxCores <- detectCores() - 1
    if (numberOfIterations < maxCores) {
        nCores = numberOfIterations
    } else {
        nCores = maxCores
    }

    print(sprintf('Opening %s cores...', nCores))
    cl <- makeCluster(nCores)
    #clusterExport(cl, c('get_timing_input'))
    clusterCall(cl, function() {library('hrdtimR'); library('VariantAnnotation'); library('plyr')})
    return(cl)
}


paths = read.table(args$paths, header=TRUE, sep='\t', stringsAsFactors=FALSE)
paths = paths[! file.exists(paths[, 4]), ]

paths = ddply(paths, 'timing_input', function(z) {
      return(z[dim(z)[1], ])
})

print(paths)

cl = open_cores(dim(paths)[1])

parApply(cl, paths, 1, function(row) {
    output_path = row[4]
    progress_file = paste0(row[4], '.in-progress.txt')
    if (!file.exists(output_path)) {
        file.create(progress_file)

        pog_id = strsplit(row[1], '_')[[1]][1]

        logString = sprintf('Sending %s to %s', pog_id, output_path)
        write(logString, file=progress_file, append=TRUE)

        timing_input = get_timing_input(row[2], row[3])
        save(timing_input, file=output_path)

        unlink(progress_file)
    }
})


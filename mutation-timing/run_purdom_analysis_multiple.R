#!/usr/bin/env Rscript

' run_purdom_analysis.R

Usage: run_purdom_analysis.R -p PATHS

Options:
    -p --paths PATHS    Path to a file containing path information
' -> doc

library('docopt')
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
    clusterCall(cl, function() {library('hrdtimR'); library('VariantAnnotation'); library('plyr')})
    return(cl)
}

args = docopt(doc)

paths = read.table(args$paths, header=TRUE, sep='\t', stringsAsFactors=FALSE)

cl = open_cores(dim(paths)[1])

parApply(cl, paths, 1, function(row) {
    if (file.exists(row[4])) {
        pog_id = strsplit(row[1], '_')[[1]][1]
        output_path = paste0(row[4], '.output.Rdata')
        print(sprintf('Sending %s to %s', pog_id, output_path))
        load(row[4]) # Loads variable timing_input
        timing_output = event_timing(timing_input$merged, as.numeric(row[5]))
        save(timing_output, file=output_path)
    }
})


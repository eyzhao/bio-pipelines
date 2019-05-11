' output_table.R

Usage: output_table.R -i INPUT -o OUTPUT [ -p PLOT ]

Options:
    -i --input INPUT        Path to output object from sciClone (.RDS)
    -o --output OUTPUT      Path to write the SciClone table to (.tsv)
    -p --plot PLOT          Path to write the SciClone figure to (.pdf) - optional
' -> doc

library(docopt)
args <- docopt(doc)

library(sciClone)

input <- readRDS(args[['input']])
if (is.null(input)) {
    file.create(args[['output']])
    if (!is.null(args[['plot']])) {
        file.create(args[['plot']])
    }
} else {
    writeClusterTable(input, args[['output']])
    if (!is.null(args[['plot']])) {
        sc.plot1d(input, args[['plot']])
    }
}

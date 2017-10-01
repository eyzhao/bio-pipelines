' Usage: nnls_monte_carlo.R -v VCFPATH [ -o OUTPUTPATH -p PLOTPATH -f ]

Options:
    -v VCFPATH          Path to table of VCF paths (two columns: sample_id and snv_path)
    -o OUTPUTPATH       Path to output file, where NNLS results will be written to as a table
    -f                  Fraction. If present, will normalize mutation signatures to add to 1.
' -> doc

library('docopt')
args <- docopt(doc)

library('devtools'); load_all('/home/ezhao/Projects/signatures/svn/hrdtools')
library('ggplot2')
library('plyr')
library('nnls')
library('readr')
library('IRanges')
library('GenomicRanges')
library('VariantAnnotation')
library('BSgenome')
library('SomaticSignatures')

paths <- read_tsv(args[['v']])

exposures <- ddply(paths, 'sample_id', function(z) {
    z <- z[1, ]
    output <- run_snv(vcf_file = z$snv_path, fractions = args[['f']], iterations = 1000)
    return(data.frame(names = output$names,
                      fit = output$x,
                      means = output$means,
                      lCI = output$lCI,
                      uCI = output$uCI))
})

if (!is.null(args[['o']])) {
    write_tsv(as.data.frame(exposures), args[['o']])
}

' Usage: nnls_and_plot.R -v VCFPATH [ -o OUTPUTPATH -p PLOTPATH -f ]

Options:
    -v VCFPATH          Path to input VCF file of somatic SNVs
    -o OUTPUTPATH       Path to output file, where NNLS results will be written to as a table
    -p PLOTPATH         Path to output file, where a PDF output plot will be generated
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

exposures <- run_snv(vcf_file = args[['v']], fractions = args[['f']])
plot <- plot_nnls_exposures(exposures)

if (!is.null(args[['o']])) {
    write_tsv(as.data.frame(exposures), args[['o']])
}

if (!is.null(args[['p']])) {
    pdf(args[['p']], height=6, width=2.5,)
    print(plot)
    dev.off()
}

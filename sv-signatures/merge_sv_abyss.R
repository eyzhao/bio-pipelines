' Usage: merge_sv_abyss.R -p PATHSFILE -o OUTPUT
' -> doc

library('docopt')
args <- docopt(doc)

library(tidyr)
library(plyr)
library(VariantAnnotation)

paths <- read_tsv(args[['PATHSFILE']])

merged <- ddply(paths, c('id', 'sample_prefix'), function(p) {
    print(paste('Running on files:', p$delly_path, 'and', p$abyss_path))
    vcf <- readVcf(p$delly_path, 'hg19')
    table <- cbind(as.data.frame(rowRanges(vcf)), as.data.frame(info(vcf)))
    somatic <- table[table$Somatic, ]

    abyss <- read_tsv(p$abyss_path)[, c('breakpoint')]
    abyss_tidy <- separate(data = abyss, col = breakpoint, into = c('chr1', 'pos1', 'chr2', 'pos2'), sep = '[\\|\\:]')

    match <- ddply(somatic, c('seqnames', 'start'), function(z) {
        z <- z[1, ]
        m <- abyss_tidy[((abyss_tidy$chr1 == z$seqnames & abyss_tidy$chr2 == z$CHR2)
                        | (abyss_tidy$chr2 == z$seqnames & abyss_tidy$chr1 == z$CHR2))
                        & (abs(as.numeric(z$start) - as.numeric(abyss_tidy$pos1)) < 20 
                           | abs(as.numeric(z$start) - as.numeric(abyss_tidy$pos2))  < 20 )
                        & (abs(as.numeric(z$END) - as.numeric(abyss_tidy$pos2)) < 20 
                           | abs(as.numeric(z$END) - as.numeric(abyss_tidy$pos1)) < 20), ]
        if(dim(m)[1] > 0) {
            return(cbind(z, m[1, ]))
        }
    })

    return(match[, c('chr1', 'pos1', 'chr2', 'pos2', 'SVLEN', 'SVTYPE')])
})

write_tsv(merged, path = args[['OUTPUT']])

' Usage: microhomology_score.R ( -i INDELGLOB | -p PATHSFILE ) -o OUTPUT
' -> doc

library(docopt)
args <- docopt(doc)

library(plyr)
library(readr)
library(tibble)
library(BSgenome)
library(VariantAnnotation)

rev_char <- function(char) {
    sapply(char, function(a) {
           paste(rev(substring(a,1:nchar(a),1:nchar(a))),collapse="")
    })
}

get_match_length <- function(x, y) { 
    min_length <- min(nchar(x), nchar(y))
    if (min_length == 0) {
        return(0) 
    } else { 
        return(sum(sapply(1:min_length, function(i) {
            substr(x, 1, i) == substr(y, 1, i) 
        })))
    }
}

get_overlap_length <- function(seq1, seq2) {
  overlap_matrix <- sapply(1:nchar(seq1), function(i) {
    subseq1 <- substring(seq1, i, nchar(seq1))
    sapply(1:nchar(seq2), function(j) {
      subseq2 <- substring(seq2, j, nchar(seq2))
      get_match_length(subseq1, subseq2)
    })
  })
  
  return(max(overlap_matrix))
}

parse_indel <- function(indel_vcf_path) {
    vcf <- readVcf(indel_vcf_path)
    vr <- rowRanges(vcf)
    rowsToKeep <- elementNROWS(vr$ALT) == 1
    vr$ALT <- unlist(CharacterList(vr$ALT))
    all_indels <- as_tibble(vr[rowsToKeep, ])

    indels <- all_indels[paste0('chr', all_indels$seqnames) %in% seqlevels(hg19)
                         & nchar(all_indels$REF) > nchar(all_indels$ALT), ]
    return(indels)
}

get_flanking_sequence <- function(indel_vcf_path) {
    print(indel_vcf_path)
    indels <- parse_indel(as.character(indel_vcf_path))
    print(paste('Number of indels:', dim(indels)[1]))
    print('retrieving flanking sequences')
    five_prime_flank <- getSeq(hg19, paste0('chr', indels$seqnames), indels$start - 25, indels$start)
    three_prime_flank <- getSeq(hg19, paste0('chr', indels$seqnames), indels$end + 1, indels$end + 26)
    deleted <- as.character(sapply(indels$REF, function(z) {substring(z, 2, nchar(z))}))
    return(cbind(indels, tibble(deleted = as.character(deleted),
                                five_prime_flank = as.character(five_prime_flank), 
                                three_prime_flank = as.character(three_prime_flank))))
}


hg19 <- getBSgenome('BSgenome.Hsapiens.UCSC.hg19')

if (args[['i']]) {
    indel_paths <- Sys.glob(args[['INDELGLOB']])
} else if (args[['p']]) {
    indel_paths <- read_lines(args[['PATHSFILE']])
} else {
    stop('ERROR: Must provide either -p or -i flag')
}

indel_paths_df <- data.frame(id = gsub('.*?(POG\\d\\d\\d).*', '\\1', indel_paths),
                       sample_prefix = gsub('.*?([bioparch]+\\d).*', '\\1', indel_paths),
                       indel_path = indel_paths)

indel_flanked <- ddply(indel_paths_df, c('id', 'sample_prefix'), function(z) {
    p <- z$indel_path[length(z$indel_path)]
    flanked <- get_flanking_sequence(as.character(p))
})

write_tsv(indel_flanked, path = args[['OUTPUT']])
print(paste('Wrote INDEL flanking data to', args[['OUTPUT']]))

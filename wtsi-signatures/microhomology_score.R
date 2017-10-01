' Usage: microhomology_score.R -i INDELGLOB -d DELLYPATHFILE -m METADATAFILE -o OUTPUT
' -> doc

library(docopt)
args <- docopt(doc)

library(readr)
library(tibble)
library(BSgenome)
library(VariantAnnotation)

rev_char <- function(char) {
    sapply(char, function(a) {
           paste(rev(substring(a,1:nchar(a),1:nchar(a))),collapse="")
    })
}

match_length <- function(x, y) { 
    min_length <- min(nchar(x), nchar(y))
    if (min_length == 0) {
        return(0) 
    } else { 
        return(sum(sapply(1:min_length, function(i) {
            substr(x, 1, i) == substr(y, 1, i) 
        })))
    }
}

any_match_length <- function(seq1, seq2) {
    match <- TRUE
    index <- 0
    while (match) {
        index <- index + 1
        if (index > nchar(seq1) || index > nchar(seq2)) {
            match <- FALSE
        } else {
            match <- substring(seq1, 1, index) == substring(seq2, 1, index)
        }
    }
    return(index - 1)
}

parse_indel <- function(indel_vcf_path) {
    vcf <- readVcf(indel_vcf_path)
    vr <- rowRanges(vcf)
    rowsToKeep <- elementNROWS(vr$ALT) == 1
    vr$ALT <- unlist(CharacterList(vr$ALT))
    all_indels <- as_tibble(vr[rowsToKeep, ])

    indels <- all_indels[all_indels$width >= 3
                         & paste0('chr', all_indels$seqnames) %in% seqlevels(hg19)
                         & nchar(all_indels$REF) > nchar(all_indels$ALT), ]
    return(indels)
}

parse_sv <- function(sv_vcf_path) {
    vcf <- readVcf(sv_vcf_path)
    vr <- rowRanges(vcf)
    info <- info(vcf)

    rowsToKeep <- elementNROWS(vr$ALT) == 1 & info$PRECISE & info$Somatic
    vr <- vr[rowsToKeep, ]
    info <- info[rowsToKeep, ]
    vr$ALT <- unlist(CharacterList(vr$ALT))

    del_rows <- vr$ALT == '<DEL>'
    del_vr <- vr[del_rows, ]
    del_info <- info[del_rows, ]

    return(as_tibble(GRanges(data.frame(chr = seqnames(del_vr), start = start(del_vr), end = del_info$END))))
}

get_microhomology_score <- function(indel_vcf_path, sv_vcf_path = NULL) {
    # This function can accept SV path to account for SV deletion calls from DELLY
    # However, currently, only INDELS are being used in this analysis

    print(indel_vcf_path)
    indels <- parse_indel(as.character(indel_vcf_path))

    if (!is.null(sv_vcf_path)) {
        print(paste('Calculating Microhomology Score for', indel_vcf_path, 'and', sv_vcf_path))
        sv <- parse_sv(as.character(sv_vcf_path))
        print(paste('Number of precise somatic sv deletions:', dim(sv)[1]))

        deletions <- rbind(indels[, c('seqnames', 'start', 'end')], sv[, c('seqnames', 'start', 'end')])
        deletions <- unique(deletions[order(deletions$seqnames, deletions$start, deletions$end), ])
    } else {
        print(paste('Calculating Microhomology Score for', indel_vcf_path))
        deletions <- unique(indels[, c('seqnames', 'start', 'end')])
    }

    print(paste('Number of indels:', dim(indels)[1]))

    print(paste('Total number of deletions:', dim(deletions)[1]))

    print('retrieving flanking sequences')

    right_flank_1 <- getSeq(hg19, paste0('chr', deletions$seqnames), deletions$start + 1, deletions$start + 101)
    right_flank_2 <- getSeq(hg19, paste0('chr', deletions$seqnames), deletions$end + 1, deletions$end + 101)

    left_flank_1 <- getSeq(hg19, paste0('chr', deletions$seqnames), deletions$start - 100, deletions$start)
    left_flank_2 <- getSeq(hg19, paste0('chr', deletions$seqnames), deletions$end - 100, deletions$end)

    print('computing match statistics')

    right_match <- apply(data.frame(seq1 = right_flank_1, seq2 = right_flank_2), 1, function(z) {
        match_length(z['seq1'], z['seq2'])
    })

    left_match <- apply(data.frame(seq1 = rev_char(as.character(left_flank_1)),
                                   seq2 = rev_char(as.character(left_flank_2))), 1, function(z) {
        match_length(z['seq1'], z['seq2'])
    })

    max_match <- apply(data.frame(left_match, right_match), 1, max)
    print(paste('Number of microhomology events:', sum(max_match > 4)))

    microhomology_score <- sum(max_match > 4) / length(max_match)
    print(paste('Microhomology proportion:', microhomology_score))

    return(data.frame(count = sum(max_match > 4), total = length(max_match), proportion = microhomology_score))
}

hg19 <- getBSgenome('BSgenome.Hsapiens.UCSC.hg19')
indel_paths <- Sys.glob(args[['INDELGLOB']])
delly_paths <- readLines(args[['DELLYPATHFILE']])
metadata <- read_tsv(args[['METADATAFILE']])

indel_paths_df <- data.frame(id = gsub('.*?(POG\\d\\d\\d).*', '\\1', indel_paths),
                       sample_prefix = gsub('.*?([bioparch]+\\d).*', '\\1', indel_paths),
                       indel_path = indel_paths)

delly_paths_df <- data.frame(id = gsub('.*?(POG\\d\\d\\d).*', '\\1', delly_paths),
                             library1 = gsub('.*?POG\\d\\d\\d\\_(.*?)\\_.*', '\\1', delly_paths),
                             library2 = gsub('.*?POG\\d\\d\\d\\_.*?\\_(.*?).delly.vcf', '\\1', delly_paths),
                             delly_path = delly_paths)

delly_metadata_1 <- merge(delly_paths_df, metadata[, c('library_name', 'sample_prefix')], by.x = 'library1', by.y = 'library_name')[, c('id', 'sample_prefix', 'delly_path')]
delly_metadata_2 <- merge(delly_paths_df, metadata[, c('library_name', 'sample_prefix')], by.x = 'library2', by.y = 'library_name')[, c('id', 'sample_prefix', 'delly_path')]
delly_metadata <- rbind(delly_metadata_1, delly_metadata_2)
delly_metadata_indel <- merge(delly_metadata, indel_paths_df, by = c('id', 'sample_prefix'))

print(head(delly_metadata_indel))
microhomology <- do.call('rbind', apply(delly_metadata_indel, 1, function(z) { get_microhomology_score(z['indel_path'], z['delly_path']) }))

out <- cbind(id = delly_metadata_indel$id, prefix = delly_metadata_indel$sample_prefix, microhomology)
write_tsv(out, path = args[['OUTPUT']])
print(paste0('Wrote output to ', args[['OUTPUT']]))

' Usage: microhomology_score.R -i INDELSUMMARY -o OUTPUT
' -> doc

library(docopt)
args <- docopt(doc)

library(readr)
library(tibble)
library(plyr)

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

microhomology_summary <- function(indel_table) {
    print('computing match statistics')

    three_prime <- apply(indel_table, 1, function(z) {
        get_match_length(z['deleted'], z['three_prime_flank'])
    })

    five_prime <- apply(indel_table, 1, function(z) {
        get_match_length(rev_char(z['deleted']), rev_char(z['five_prime_flank']))
    })

    max_match <- apply(cbind(three_prime, five_prime), 1, max)

    indel_table$microhomology_length <- max_match
    indel_table$is_microhomology <- max_match > 2

    return(indel_table)
}

indels <- read_tsv(args[['INDELSUMMARY']])
indels <- indels[!is.na(indels$deleted)
                 & !is.na(indels$three_prime_flank)
                 & !is.na(indels$five_prime_flank), ]

indel_microhomology <- microhomology_summary(indels)

mh_scores <- ddply(indel_microhomology, c('id', 'sample_prefix'), function(z) {
    print(paste(z[['id']][1], z[['sample_prefix']][1]))
    mh_count <- sum(z$is_microhomology)
    indel_count <- length(z$is_microhomology)
    mh_score <- mh_count / indel_count
    data.frame(mh_count, indel_count, mh_score)
})

write_tsv(mh_scores, args[['OUTPUT']])

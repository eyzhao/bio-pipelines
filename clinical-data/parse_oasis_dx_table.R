' Usage: parse_oasis_dx_table.R -i INPUT -b BIRTHDATES -d DXDATES
' -> doc

library(docopt)
args <- docopt(doc)

library(readr)
library(plyr)

table <- read_tsv(args[['INPUT']])
table <- table[ !is.na(table$`Gsc Pog Id`), ]
table <- table[grepl('\\d\\d\\d', table$`Gsc Pog Id`), ]

table$gsc_pog_id <- paste0('POG', gsub('.*?(\\d\\d\\d).*', '\\1', table$`Gsc Pog Id`))

for (date_col in colnames(table)[grepl('Date', colnames(table))]) {
    table[[date_col]] <- parse_date(table[[date_col]], '%b %d %Y')
}

birthdate_table <- unique(table[, c('gsc_pog_id', 'Birth Date')])
colnames(birthdate_table) <- c('gsc_pog_id', 'birthdate')
birthdate_table <- birthdate_table[!is.na(birthdate_table$birthdate), ]

write_tsv(birthdate_table, args[['BIRTHDATES']])

all_dx <- table[ !is.na(table$`Diagnosis Date 1`), ]

dx_table <- all_dx[, c('gsc_pog_id', paste(rep(c("Diagnosis Date", "Total Age At Diagnosis", "Tumour Group"), 6), sort(rep(1:6, 3))))]

dx_list <- list(dx_table[, c(1, 2:4)], dx_table[, c(1, 5:7)], dx_table[, c(1, 8:10)], dx_table[, c(1, 11:13)])
dx_list <- lapply(dx_list, function(z) {
    colnames(z) <- c('gsc_pog_id', 'dx_date', 'dx_age', 'tumour_group'); return(z)
})

dx_df <- do.call('rbind', dx_list)
dx_df <- dx_df[!is.na(dx_df$dx_date), ]

dx_df <- ddply(dx_df, 'gsc_pog_id', function(z) { z <- z[order(z$dx_date), ]; z$dx_id <- 1:dim(z)[1]; return(z) })

write_tsv(dx_df, args[['DXDATES']])

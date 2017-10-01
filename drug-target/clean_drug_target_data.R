' Usage: clean_drug_target_data.R -i INPUT -o OUTPUT
' -> doc

library('docopt')
args <- docopt(doc)

table <- read.table(args[['INPUT']], header=T, sep='\t', stringsAsFactors=F)

clean <- unique(table[table$hugo == table$gene_name, ])
clean <- clean[!is.na(clean$seqnames), ]
clean <- clean[, !(names(clean) == 'gene_name')]

write.table(clean, args[['OUTPUT']], col.names = T, row.names = F, sep='\t', quote=F)

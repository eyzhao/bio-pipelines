' Usage: get_clinical_data.R -i INPUT -o OUTPUT -c CASELIST -b BCCA
' -> doc

library('docopt')
args <- docopt(doc)

clinical_table <- read.table(args[['INPUT']], sep='\t', header=T)
caselist <- read.table(args[['CASELIST']], header=F, col.names=c('gsc_pog_id'))
bcca_table <- read.table(args[['BCCA']], sep='\t', header=T, quote='"')

merged <- merge(caselist, clinical_table, by='gsc_pog_id')
merged <- merge(bcca_table[, c('id', 'bcca')], merged, by.x = 'id', by.y = 'gsc_pog_id')

write.table(merged, file = args[['OUTPUT']], sep='\t', col.names = T, row.names = F, quote=F)

' Usage: make_germline_bam_query_file.R -i GERMLINEPATH -b BAMDATA -o OUTPUT
' -> doc

library('docopt')
library('dplyr')

args <- docopt(doc)

germline_table <- read.table(args[['GERMLINEPATH']], header=T, sep='\t', quote="")
bam_table <- read.table(args[['BAMDATA']], header=T, sep='\t', quote="")

bam_table$normal_id <- paste(bam_table$pog_id, bam_table$normal_library, sep='_')
bam_table$tumour_id <- paste(bam_table$pog_id, bam_table$tumour_library, sep='_')
germline_df <- merge(germline_table, bam_table[, c('normal_id', 'tumour_id', 'normal_bam')], by.x='Sample.Name', by.y='normal_id')
somatic_df <- merge(germline_table, bam_table[, c('normal_id', 'tumour_id', 'tumour_bam')], by.x='Sample.Name', by.y='normal_id')
germline_df <- rename(germline_df, normal_id = Sample.Name); somatic_df <- rename(somatic_df, normal_id = Sample.Name)

germline_df$bam_type = 'normal'
somatic_df$bam_type = 'tumour'

ncol <- dim(germline_df)[2]
germline_df <- germline_df[, c(2:5, ncol-1, ncol-2, 1, ncol, 6:(ncol-4))]; names(germline_df)[5] <- 'bam_path'
somatic_df <- somatic_df[, c(2:5, ncol-1, ncol-2, 1, ncol, 6:(ncol-4))]; names(somatic_df)[5] <- 'bam_path'

complete_table <- rbind(germline_df, somatic_df)
complete_table <- complete_table[order(complete_table$normal_id, complete_table$Chr, complete_table$Pos), ]                                           

write.table(complete_table, file=args[['OUTPUT']], col.names=T, row.names=F, quote=F, sep='\t')

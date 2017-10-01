d <- readRDS('rdata/POG_BR_nnls_exposures.RData')

summary <- ldply(d, function(z) {data.frame(exposure=z$x[3], mean=z$means[3], lCI=z$lCI[3], uCI=z$uCI[3], signature.name=z$names[3])},
                 .id = 'gsc_pog_id')

summary <- summary[!grepl('arch', summary$gsc_pog_id), ]
summary$gsc_pog_id <- gsub('.*?(POG\\d\\d\\d).*', '\\1', summary$gsc_pog_id)

tx <- read.table('../../pog/clinical/20170123_pog_treatment_table.txt', header=TRUE, sep='\t', stringsAsFactors=F)

tx$course_begin_on <- as.Date(tx$course_begin_on); tx$course_end_on <- as.Date(tx$course_end_on)
candidates <- tx[grepl('PLATIN', tx$drug_list) & tx$pog_tumour_group == 'BRC' & tx$course_end_on != tx$course_begin_on, ]

ttf <- ddply(candidates, 'gsc_pog_id', function(z) {z <- z[z$course_begin_on > as.Date('2012-01-01'), ]; data.frame(ttf = sum(z$course_end_on - z$course_begin_on))}) # Change this to filter based on diagnosis date, and later potentially biopsy date too

m <- merge(summary, ttf, by='gsc_pog_id')

pdf('test.pdf', width=4, height=4); plot(m$mean, as.numeric(m$ttf)); dev.off();
cox <- coxph(Surv(as.numeric(m$ttf)) ~ m$exposure)
print(summary(cox))

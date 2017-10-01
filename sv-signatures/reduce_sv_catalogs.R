' Usage: reduce_sv_catalogs.R -i INPUT -r REDUCTION -o OUTPUT
' -> doc

library(docopt)

args <- docopt(doc)

library(plyr)

catalogs <- readRDS(args[['INPUT']])
reduction <- readLines(args[['REDUCTION']])

catalogs <- llply(catalogs, function(df) {
    return(df[df$SVTYPE == 'TRA' | df$subtype != 0, ])
})

print('Reduction classes:')
print(cbind(catalogs[[1]][, c('clustered', 'SVTYPE', 'subtype')], reduction))

new <- llply(catalogs, function(df) {
    df$concat <- reduction
    new_df <- ddply(df, 'concat', function(z) {
        data.frame(subtype = z$subtype[1], SVTYPE = z$SVTYPE[1], clustered = z$clustered[1], value = sum(z$value))
    })
    return(new_df[, c('subtype', 'SVTYPE', 'clustered', 'value')])
})

saveRDS(new, args[['OUTPUT']])

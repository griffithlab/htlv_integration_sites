library(ggplot2)
library(cowplot)
library(reshape2)
library(gtools)

basedir = '/storage1/fs1/mgriffit/Active/griffithlab/adhoc/ratner_p01/htlv_integration_sites/beds/'
bed_files = c('CTCF_1' = paste0(basedir, 'CTCF-1.markedsorted_with_hits_to_viral_filtered_merged.bed'), 
              #'CTCF_3' = paste0(basedir, 'CTCF-3.markedsorted_with_hits_to_viral_filtered_merged.bed'), skipped, no reads
              'CTCF_7' = paste0(basedir, 'CTCF-7.markedsorted_with_hits_to_viral_filtered_merged.bed'),
              'CTCF_8' = paste0(basedir, 'CTCF-8.markedsorted_with_hits_to_viral_filtered_merged.bed'),
              'P12_10B' = paste0(basedir, 'P12-10B.markedsorted_with_hits_to_viral_filtered_merged.bed'),
              'P12_14' = paste0(basedir, 'P12-14.markedsorted_with_hits_to_viral_filtered_merged.bed'),
              'P12_5' = paste0(basedir, 'P12-5.markedsorted_with_hits_to_viral_filtered_merged.bed'),
              'P12_8' = paste0(basedir, 'P12-8.markedsorted_with_hits_to_viral_filtered_merged.bed'))

samples = lapply(bed_files, read.table, 
                 col.names = c('chr', 'start', 'end', 'reads'))

samples = lapply(samples, function(x) {x[x$chr %in% paste0('chr', c(seq(1, 24), 'X', 'Y')), ]})
samples = lapply(samples, function(x) {x[order(x$reads), ]})
samples = lapply(samples, function(x) {x$reads = cumsum(x$reads); return(x)})
samples = lapply(samples, function(x) {x$loc = seq(1 : nrow(x)); return(x)})
samples = lapply(samples, function(x) {x$loc_prop = x$loc / max(x$loc); return(x)})
samples = lapply(samples, function(x) {x$diagonal = seq(0, tail(x$reads, 1), tail(x$reads, 1) / nrow(x))[1:nrow(x)]; return(x)})
samples = lapply(samples, function(x) {x$reads_prop = x$reads / max(x$reads); return(x)})
samples = lapply(samples, function(x) {x$diagonal_prop = x$diagonal / max(x$diagonal); return(x)})

plots = list()
for (sample in names(samples)) {
  x = samples[[sample]]
  subtitle_string = paste(paste0('Gini index: ', 
                          round(abs(sum(x$reads_prop - x$diagonal_prop)) / (nrow(x) / 2), 2)),
                          '\n',
                          paste0('# merged locations: ',
                          nrow(x),
                          sep = ''), sep = '')
  plots[[sample]] = 
    ggplot(x, aes(x = loc_prop, 
                  y = reads_prop, group = 1)) + geom_line() +
    theme_classic() +
    geom_line(aes(x = loc_prop, y = diagonal_prop), inherit.aes = F, linetype = 'dotted') + 
    geom_ribbon(aes(ymax = diagonal_prop, ymin = reads_prop), fill = 'pink4', alpha = 0.25) +
    ggtitle(sample, subtitle = subtitle_string) + xlab('Proportion of locations') +
    ylab('Proportion of reads')
}

png('/storage1/fs1/mgriffit/Active/griffithlab/adhoc/ratner_p01/htlv_integration_sites/slides-tables/combined_gini_plots.png', width = 1000, height = 2000)
plot_grid(plotlist = plots, ncol = 2)
dev.off()

library('ggplot2')
library('reshape2')
library('grid')
library('scales')
library('dplyr')
setwd('~/R_scripts')
cov <- read.csv('WGScoverages.csv')
#cov.m <- melt(cov) # skips removing any columns
cov.shrunk <- transmute(cov, SampleID, 
                        X10x/100, 
                        X15x/100,
                        X20x/100,
                        X40x/100,
                        X50x/100) # convert values to be between [0,1] instead of [0,100]
cov.sub = cov.shrunk[, c(1,3,4,5)] # drop 10x and 50x columns
cov.msub = melt(cov.sub)

p <- ggplot(cov.msub, aes(SampleID, value))
p <- p + geom_bar(aes(fill = variable), width = .7, stat = 'identity',
           position = position_dodge(width = .7)) +
     theme_classic() + # drop gridlines and grey background
     theme(plot.title = element_text(size = 8),
           legend.title = element_text(size = 6),
           legend.text = element_text(size = 6),
           legend.key.size = unit(.2, 'cm'),
           axis.text.x = element_text(angle = 90, vjust = .5, size = 6),
           axis.text.y = element_text(size = 6),
           axis.title.x = element_text(size = 8),
           axis.title.y = element_text(size = 8)) +
     labs(title = 'WGS', x = 'Sample ID', y = 'Coverage') +
     scale_y_continuous(labels = percent, expand = c(0,0), 
                        limits = c(0,1)) + # convert to %, push graph area out to y-axis, include 0 and 1
     scale_x_discrete(labels = abbreviate, expand = c(0,0)) + # push graph area out to x-axis
     scale_fill_discrete(name = 'Read Depth', 
                         labels = c('15x','20x','40x')) # change as needed

ggsave('WGScov_plot.png', p, width = 8, height = 3, dpi = 1200)
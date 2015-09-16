scales
ggplot2
reshape2

bfreqs <- read.table('vcounts_bin10kb.tsv', header = TRUE)
bfreqs.m <- melt(bfreqs)
bfreqs.m$bins <- as.character(bfreqs.m$variable)
bfreqs.m$bins <- gsub("^X", "", bfreqs.m$bins)
p <- ggplot(bfreqs.m, aes(bins, value))
p <- p + geom_bar(aes(fill = bin_start), width = .7, stat = 'identity',
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
  labs(title = '', x = '<-- START     Region on TTN gene     END -->', y = 'Variant frequency') +
  scale_y_continuous(labels = percent, expand = c(0,0), 
                     limits = c(0,0.1)) + # convert to %, push graph area out to y-axis, include 0 and 10%
  scale_x_discrete(expand = c(0,0)) + # push graph area out to x-axis
  scale_fill_discrete(name = 'Cohort', 
                      labels = c('ESP6500','SIM WGS')) # change as needed

ggsave('bin10kb_dodgedplot.png', p, width = 8, height = 3, dpi = 1200)
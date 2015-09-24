library("ggplot2")

setwd('~/R_scripts')
reads3 <- reads <- read.csv("WESread_depths.csv")
#reads3$SampleID <- factor(reads$SampleID, 
#                           levels=reads3[order(reads$Read_Depth), 
#                                          "SampleID"]) # orders bars by increasing values
readplot <- ggplot(reads, aes(x = SampleID, y = Read_Depth))
  
readplot <- readplot + geom_bar(stat = "identity", width = .75, 
                     color = "black", fill = "#006699") +
            theme_classic() +
            labs(title = 'WES', x = "Sample ID", y = "Average Read Depth") +
            scale_x_discrete(labels = abbreviate, expand = c(0,0)) +
            scale_y_continuous(labels = c("0","10x","20x","30x","40x","50x","60x","70x","80x","90x","100x","110x","120x"),
                               breaks = seq(0,120,10), 
                               expand = c(0,0), limit = c(0,120)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9))

ggsave('WESrd_plot.png', readplot, width = 8, height = 3, dpi = 1200)
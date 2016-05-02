create.dotmap <- function (x, bg.data = NULL, filename = NULL, main = NULL, alpha = 1, 
                           continuous.bg = TRUE,
                           legend.position = "left",
                           bg.legend.title = "P value",
                           bg.palette = NULL,
                           bg.title.position = "top",
                           bg.fill.low = "black", 
                           bg.fill.high = "white",
                           bg.legend.title.text.size = 10, 
                           bg.legend.title.fontface = "bold", 
                           bg.legend.title.colour = "black", 
                           bg.legend.title.angle = 0,
                           dot.legend.title = "Direction", 
                           dot.legend.label = c("Up", "Down"),
                           dot.legend.colour = c("red", "blue"),
                           dot.size.title = "Fold Change",
                           dot.legend.labels = c("Positive", "Negative"), 
                           dot.legend.colours = c("red", "blue"), 
                           scale_size = c(0, 15), 
                           x.axis.label = NULL, 
                           y.axis.label = NULL, 
                           height = 7, width = 7, unit = 'in', resolution = 600) {

	bg.data$name = rownames(bg.data)
	bg.data <- bg.data[rev(bg.data$name), ]
	bg.data$name <- factor(bg.data$name, levels = bg.data$name)
	bg.data.m <- melt(bg.data, id.vars = "name")
	
	fg.data <- x
	fg.data$name = rownames(fg.data)
	fg.data <- fg.data[rev(fg.data$name), ]
	fg.data.m <- melt(fg.data, id.vars = "name")
	
	po.nopanel <- list(theme(panel.background = element_blank(), 
                             panel.grid.minor = element_blank(), 
                             panel.grid.major = element_blank(), 
                             axis.text.x = element_text(vjust = 0.5, angle = 90)))

	p <- ggplot(bg.data.m, aes(variable, name), environment = environment())
	
	if(continuous.bg) {
	  fill.colour <- bg.data.m$value
		effect.size <- abs(fg.data.m$value)
		fg.colour <- ifelse(fg.data.m$value > 0, "red", "blue");
		fg.colour.na <- ifelse(is.na(fg.data.m$value), "black", NA)
		
		 
    p <- p + geom_tile(aes(fill = fill.colour), colour = "black", alpha = alpha)
		p <- p + geom_point(aes(size = effect.size, colour = fg.colour), breaks = 10)
		p <- p + geom_point(aes(shape = 4), colour = fg.colour.na, size = 12)
		p <- p + scale_shape_identity()
		p <- p + scale_size(range = scale_size, breaks = seq(0, 3, 0.5), limits = c(0, 3))
		p <- p + scale_fill_gradient(low = bg.fill.low, 
									 high = bg.fill.high, 
									 limits = c(0, 1),
									 guide = "colourbar",
									 labels = c("0.0-0.2","0.2-0.4","0.4-0.6",
                                                "0.6-0.8","0.8-1.0"))
		p <- p + scale_colour_manual(name = dot.legend.title,
									 labels = dot.legend.labels, 
									 values = dot.legend.colours)
	
	}else {
	  if (NA %in% bg.data.m$value) {
	    bg.fill.count <- length(unique(bg.data.m$value)) - 1
	  }else {
	    bg.fill.count <- length(unique(bg.data.m$value))
		}
	  
	  fg.colour <- factor(fg.data.m$value)
    	  
	  p <- p + geom_tile(aes(fill = value), colour = "black", alpha = alpha)
    if (is.null(bg.palette)) {
      p <- p + scale_fill_manual(values = brewer.pal(n = bg.fill.count, name = "Set3"))
	  }else {
      p <- p + scale_fill_manual(values = bg.palette)
	  }
    p <- p + geom_point(aes(colour = fg.colour), size = 4) 
	  p <- p + scale_colour_manual(name = dot.legend.title,
	                               breaks = c("1","0","-1"),
	                               labels = dot.legend.labels,
	                               values = rev(dot.legend.colours))
	 }

  p <- p + ggplot2::theme(legend.position = legend.position,
	                        legend.key = element_rect(colour = 'black', fill = 'white'))
  p <- p + guides(fill = guide_legend(order = 3,
                                      override.aes = list(colour = NULL),
                                      title = bg.legend.title,
                                      title.position = bg.title.position,
                                      title.theme = element_text(size =   bg.legend.title.text.size,
                                                                 face =   bg.legend.title.fontface,
                                                                 colour = bg.legend.title.colour,
                                                                 angle =  bg.legend.title.angle)))
  p <- p + labs(x = x.axis.label, y = y.axis.label, size = dot.size.title)
  p <- p + po.nopanel
  
# output to file
	if(is.null(filename)){
    return(p)
  }else{
    tiff(filename, 
         height = height, 
         width = width, 
         units = unit, 
         res = resolution, 
         compression = 'lzw', 
         type = 'cairo')
	  print(p)
	  dev.off()
	}
}

\name{create.dotmap}
\alias{create.dotmap}
\title{Make a dotmap with coloured background}
\description{This plotting function is intended to plot two data sources in a single plot. One data source is used to plot a heatmap with colours showing P values of a statistical test. The other data source is used to plot coloured dots in the foreground, with dot size showing the effect size such as fold change in a microarry experiment and colours showing the direction of change, such as up- or down-regulations. Alternatively, the heatmap can be used to represent functional consequences of genomic variants for affected genes and dots to represent copy number status of those genes.}
\usage{
  create.dotmap(x, bg.data = NULL, filename = NULL, main = NULL, alpha = 1, 
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
                height = 7, width = 7, unit = 'in', resolution = 600)
}
\arguments{
  \item{x}{A data.frame for the dotmap.}
  \item{bg.data}{A data.frame for background data. It must have the same dimension as \dQuote{x}. }
  \item{filename}{Filename for tiff output, or if NULL returns the ggplot2 plotting object.}
  \item{main}{Title for figure.}
  \item{alpha}{Alpha for background heatmap.}
  \item{continuous.bg}{Flag for plot type.}
  \item{legend.position}{Legend position.}
  \item{bg.legend.title}{Background data legend title.}
  \item{bg.palette}{Palette to use for discrete heatmap.}
  \item{bg.title.position}{Legend title position.}
  \item{bg.fill.low}{Low color for continuous background.}
  \item{bg.fill.high}{High color for continuous background.}
  \item{bg.legend.title.text.size}{Legend title text size.}
  \item{bg.legend.title.fontface}{Legend title font face.}
  \item{bg.legend.title.colour}{Legend title colour.}
  \item{bg.legend.title.angle}{Legend title angle.}
  \item{dot.legend.title}{Foreground direction-of-change data legend title.}
  \item{dot.legend.label}{Foreground direction-of-change data legend labels.}
  \item{dot.legend.colour}{Foreground direction-of-change data legend label colours. Must be equal to length of dot.legend.label.}
  \item{dot.size.title}{Foreground effect-size data legend title. Only available for continuous.}
  \item{dot.legend.labels}{Foreground effect-size data legend labels. Only available for continuous.} 
  \item{dot.legend.colours}{Foreground effect-size data legend label colours. Only available for continuous.}
  \item{scale_size}{Scale for size of effect bubbles.}
  \item{x.axis.label}{X-axis label.} 
  \item{y.axis.label}{Y-axis label.}
  \item{height}{Output file height.}
  \item{width}{Output file width.}
  \item{unit}{Units for height and width.}
  \item{resolution}{Resolution in dpi.}
}
\details{Fill details here.}
\seealso{}
\examples{
# for continuous plot:

set.seed(1000)

tmp <- runif(n = 100, min = 0, max = 1)
bg.data <- matrix(tmp, nrow = 10)
bg.data <- as.data.frame(bg.data)
colnames(bg.data) <- paste('Sample', 1:10, sep = '_')
rownames(bg.data) <- paste('Gene', 1:10, sep = '_')

tmp2 <- runif(n = 100, min = -3, max = 3)
fg.data <- matrix(tmp2, nrow = 10)
fg.data <- as.data.frame(fg.data)
fg.data[1,1] <- NA
fg.data[8,7] <- NA
colnames(fg.data) <- paste('Sample', 1:10, sep = '_')
rownames(fg.data) <- paste('Gene', 1:10, sep = '_')


create.dotmap(
  x = fg.data,
  bg.data = bg.data,
	filename = 'Dotmap_continuous.tiff',
	scale_size = c(0, 12),
	dot.legend.title = "Direction",
	dot.legend.labels = c("Up", "Down"),
	dot.legend.colour = c("red", "blue"),
	x.axis.label = "Sample",
	y.axis.label = "Gene",
	alpha = 0.8)


# for discrete plot:

set.seed(1000);
tmp <- sample(1:6, size = 100, replace = TRUE)
bg.data <- matrix(tmp,nrow=10)
bg.data <- as.data.frame(bg.data)
bg.data[bg.data=='1'] <- 'A'
bg.data[bg.data=='2'] <- 'B'
bg.data[bg.data=='3'] <- 'C'
bg.data[bg.data=='4'] <- 'D'
bg.data[bg.data=='5'] <- 'E'
bg.data[bg.data=='6'] <- NA
colnames(bg.data) <- paste('Sample', 1:10, sep = '_'); 
rownames(bg.data) <- paste('Gene', 1:10, sep = '_');

tmp2 <- c(NA, -1, 0, 1)
tmp2 <- sample(tmp2, size = 100, prob = c(0.7,0.1,0.1,0.1), replace = TRUE);
fg.data <- matrix(tmp2, nrow = 10);
fg.data <- as.data.frame(fg.data);
colnames(fg.data) <- paste('Sample', 1:10, sep = '_');
rownames(fg.data) <- paste('Gene', 1:10, sep = '_');


create.dotmap(
  x = fg.data,
  bg.data = bg.data,
  filename = 'Dotmap_discrete.tiff',
  continuous.bg = FALSE,
  scale_size = c(0, 15),
  bg.legend.title = "Var type", 
  dot.legend.title = "Status", 
  dot.legend.labels = c("AMP", "No Change", "DEL"),
  dot.legend.colours = c("red", "yellow", "blue"),
  x.axis.label = "Sample",
  y.axis.label = "Gene",
  alpha = 0.8)
}

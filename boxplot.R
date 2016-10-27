# Load required libraries
#
library(ggplot2)
#
# 
# Custom functions
#
make_contrast_coord <- function(n) {
  tmp <- do.call(rbind,lapply(1:n, (function(i){
    do.call(rbind,lapply(1:n, (function(j){
      if(j > i) {
        c(i,j)  
      }
    })))
  })))
  tmp <- data.frame(tmp)
  colnames(tmp) <- c("str", "end")
  tmp$ave <- apply(tmp, 1, mean)
  tmp$len <- apply(tmp, 1, (function(vct){ max(vct) - min(vct) }))
  return(tmp)
}
pval_to_asterisks <- function(p_vals) {
  astk <- sapply(as.numeric(as.character(p_vals)), (function(pv){
    if(pv >= 0 & pv < 0.0001) {
      "****"
    } else if (pv >= 0 & pv < 0.001) {
      "***"
    } else if (pv >= 0 & pv < 0.01) {
      "**"
    }  else if (pv >= 0 & pv < 0.05) {
      "*"
    } else {
      NA
    }
  }))
  return(astk)
}
#
#
# Generate a random dataset: a 20 row x 4 column matrix
#
set.seed(999)
my_data <- matrix(sapply(c(1,2.2,2.9,4.2), (function(i){ rnorm(20, 30 + (5*i), 7+i)})),
                  nrow = 20,
                  ncol = 4, 
                  dimnames = list(paste("r", as.character(1:20), sep =""), c("A", "B", "C", "D")))
boxplot(my_data, col = "darkorange", pch = 19, ylim = c(0,120), main = "Standard barplot")
head(my_data)
#
# Critical step: convert the matrix into a data frame
# In the resulting data frame, the condition <<A,B,C,D>> is converted in new variable
# The 20 x 4 matrix is converted to a 80 x 2 matrix
# Each row has two cells: condition {A,B,C,D} and value {measured value}
#
my_df <- data.frame(do.call(rbind, lapply(colnames(my_data), (function(clnm){
  values <- my_data[,clnm]
  group <- rep(clnm, nrow(my_data))
  cbind(values, group)
}))), row.names = NULL)
#
#
# Now plot a basic boxplot with ggplot2
#
bp <- ggplot(my_df, aes(x = group, y = as.numeric(values), fill = factor(group)))
bp <- bp + geom_boxplot(notch = F)
bp
#
#
# define ylim, set colors and impose jitter
#
bp <- bp + ylim(c(0,120))
my_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
bp <- bp + scale_fill_manual(values=my_colors)
bp <- bp + geom_point(position = position_jitter(width = 0.35))
bp
#
#
# change labels and title and theme; adjust labels margins as needed
#
bp <- bp + labs(title="My ggplot2 Boxplot", x="Sample groups", y = "Count of something")
bp <- bp + theme_bw()
bp <- bp + theme(legend.position="none") # Remove legend
bp <- bp + theme(axis.text.x = element_text(colour="grey20",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
                 axis.text.y = element_text(colour="grey20",size=11,angle=0,hjust=1,vjust=0,face="plain"),  
                 axis.title.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
                 axis.title.y = element_text(colour="black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"),
                 plot.title = element_text(size=18, face = "bold"))
bp
#
#
# Let's complete this analysis by running an ANOVA significance test
# and adding annotation to the plot: which are the significan differences?
# What about the p-values? Let's start with the ANOVA
#
my_anova <- aov(as.numeric(values)~group, data = my_df)
my_anova <- TukeyHSD(my_anova)
my_anova <- data.frame(cbind(my_anova$group, 
                             make_contrast_coord(length(levels(my_df$group)))))
my_anova$astks <- pval_to_asterisks(my_anova$p.adj)
#
#
# my_anova includes information about the ANOVA test on our dataset.
# Let's remove non-significant differences, as those won't be plotted.
# The filtered data frame is called tiny_anova. We also order the
# rows according to the distance between boxes (in the boxplot).
# Differences between boxes far away will be annotated at the top
#
tiny_anova <- my_anova[my_anova$p.adj < 0.05,]
tiny_anova <- tiny_anova[order(tiny_anova$len, decreasing = FALSE),]
#
#
# Let's set some parameters for the annotation. Change this depending to the specific boxplot
#
lowest.y <- 85
highest.y <- 110
margin.y <- 5
#
#
# Define the positions where the p-value segments will be drawn and then annotate the p-vals
#
actual.ys <- seq(lowest.y, highest.y, length.out = nrow(tiny_anova))
tiny_anova$ys <- actual.ys
bp_ask <- bp + annotate("segment", x = tiny_anova$str, y = tiny_anova$ys, 
                    xend = tiny_anova$end, yend = tiny_anova$ys, 
                    colour = "black", size = 0.95)
bp_val <- bp_ask + annotate("text", x = tiny_anova$ave, y = (tiny_anova$ys + margin.y) , 
                    xend = tiny_anova$end, yend = tiny_anova$ys,
                    label = paste ("p-val =", format(round(tiny_anova$p.adj, 4), nsmall = 4)))
bp_ask <- bp_ask + annotate("text", x = tiny_anova$ave, y = (tiny_anova$ys + (margin.y/3)) , 
                        xend = tiny_anova$end, yend = tiny_anova$ys,
                        label = tiny_anova$astks, size = 10)
#
#
# Here are three plots:
# 1) boxplot without annotations
# 2) boxplot with p-val annotations (up to 4 decimal places)
# 3) boxplot with *** annotations
#
bp
bp_val
bp_ask
#
#
# Success! Done...
# For questions, contact me at damiano.fantini@gmail.com
# Or visit: http://www.biotechworld.it/bioinf/
#
#

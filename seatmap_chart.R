# ---------------------------------------------
#        --- ggplot2 Tutorial Series ---
#        -------- SEATMAP CHART --------           
#
# by Damiano Fantini <damiano.fantini@gmail.com>
# ---------------------------------------------
#
#
# Load libraries
#
library(ggplot2)
#
#
# Generate some data (matrix-format)
#
set.seed(999)
my_opts <- c(NA, "Missense Mutation", "Frameshift Mutation", "Nonsense Mutation", "Splicing and UTR")
my_genes <- paste("gene", 1:25, sep = "_")
my_samples <- paste("SMPL", 1:10, sep = "")
my_mat <- matrix(data = sample(x = my_opts, 
                               size = (length(my_genes) * length(my_samples)), 
                               replace = TRUE,  
                               prob = c(0.8, 0.3, 0.1, 0.1, 0.2)),
                 nrow = length(my_genes), 
                 ncol = length(my_samples), 
                 dimnames = list(my_genes, my_samples))
my_mat[1:10,1:3]
#         SMPL1                 SMPL2               SMPL3              
# gene_1  NA                    NA                  NA                 
# gene_2  "Splicing and UTR"    NA                  NA                 
# gene_3  "Missense Mutation"   "Splicing and UTR"  NA                 
# gene_4  NA                    "Missense Mutation" "Missense Mutation"
# gene_5  "Frameshift Mutation" NA                  "Missense Mutation"
# gene_6  NA                    NA                  NA                 
# gene_7  NA                    NA                  NA                 
# gene_8  NA                    "Missense Mutation" NA                 
# gene_9  NA                    "Nonsense Mutation" NA                 
# gene_10 NA                    NA                  "Missense Mutation"
#
#
# Convert data to a integer matrix we can image (still matrix-format)
#
my_levels <- unique(as.vector(my_mat))
num_mat <- apply(my_mat,2,(function(clmn){
  sapply(clmn, (function(jj){
    if (is.na(jj)) 1
    else if (jj == my_levels[2]) 2
    else if (jj == my_levels[3]) 3
    else if (jj == my_levels[4]) 4
    else if (jj == my_levels[5]) 5
  }))  
}))
head(num_mat)
#       SMPL1 SMPL2 SMPL3 SMPL4 SMPL5 SMPL6 SMPL7 SMPL8 SMPL9 SMPL10
#gene_1     1     1     1     1     1     1     1     2     2      1
#gene_2     2     1     1     2     1     1     1     1     5      1
#gene_3     3     2     1     1     1     4     3     5     3      3
#gene_4     1     3     3     3     1     1     5     1     1      3
#gene_5     4     1     3     1     4     3     2     1     1      1
#gene_6     1     1     1     1     2     1     1     3     1      2
#
#
# Plot data using the R standard plotting system (naive way)
#
my_colors <- c("gray90", "#ff7f00", "#e31a1c", "#cab2d6", "#33a02c")
image(t(num_mat), frame = FALSE, axes = FALSE,
      xlim = c(-0.2,1.85), 
      ylab = "", xlab = "", 
      col = my_colors,
      main = "Gene Mutation Status")
legend("right", 
       legend = c("wild type", my_levels[2:5]), 
       fill = my_colors)
axis(1, at = seq(0, 1, along.with = my_samples), labels = my_samples, cex.axis = 0.75, font = 2, las = 2, tick = FALSE, pos = 0.0)
axis(2, at = seq(0, 1, along.with = my_genes), labels = my_genes, cex.axis = 0.75, font = 2, las = 1, tick = FALSE, pos = -0.02)
#
#
# Fix issue about gene order
#
image(t(num_mat)[,nrow(num_mat):1], frame = FALSE, axes = FALSE,
      xlim = c(-0.2,1.85), 
      ylab = "", xlab = "", 
      col = my_colors,
      main = "Gene Mutation Status")
legend("right", 
       legend = c("wild type", my_levels[2:5]), 
       fill = my_colors)
axis(1, at = seq(0, 1, along.with = my_samples), labels = my_samples, cex.axis = 0.75, font = 2, las = 2, tick = FALSE, pos = 0.0)
axis(2, at = seq(0, 1, along.with = my_genes), labels = my_genes[length(my_genes) : 1], cex.axis = 0.75, font = 2, las = 1, tick = FALSE, pos = -0.02)
#
#
# We can easily get much better results via ggplot2
# First, we need to convert the matrix data into a data frame
# Data Variables are: gene, sample and mutation status
#
my_df <- data.frame(do.call(rbind, lapply(1:nrow(my_mat), (function(i){
  t(sapply(1:ncol(my_mat), (function(j){
    c(rownames(my_mat)[i],
      colnames(my_mat)[j],
      my_mat[i,j])  
  })))
}))))
colnames(my_df) <- c("gene","sample","status")
head(my_df)
#    gene    sample  status
#  1 gene_1  SMPL1   <NA>
#  2 gene_1  SMPL2   <NA>
#  3 gene_1  SMPL3   <NA>
#  4 gene_1  SMPL4   <NA>
#  5 gene_1  SMPL5   <NA>
#  6 gene_1  SMPL6   <NA>
#
#
# Now we can plot a heatmap-like chart (seatmap chart) via ggplot2 and geom_tile
#
p <- ggplot(my_df, aes(y=gene, x=sample))
p <- p + geom_tile(aes(fill=status), width=.875, height=.875)
p
#
#
# We need to fix column and row order here as well. 
# We can easily get this done by scaling the x and the y 
# 
p <- p + scale_y_discrete(limits=rev(unique(my_df$gene)))
p <- p + scale_x_discrete(limits=unique(my_df$sample))
p
#
#
# Now we apply some formatting and we change theme and color sets
#
p <- p + theme_minimal(base_size = 11) + labs(x = "", y = "")
p <- p + labs(title = "Gene Mutation Status")
p <- p + scale_fill_manual(values = my_colors[c(4,3,5,2)], name = "Mutation Type",
                           breaks = levels(my_df$status)[c(2,3,1,4)], 
                           na.value = "grey90")
p <- p + guides(color = guide_legend(ncol = 1)) +
  theme(legend.key = element_rect(size = 2, color = "white"),
        legend.key.size = unit(1.5, 'lines'))
p
#
#
# We complete the Gene Mutation chart by formatting text and defining some margins
#
p <- p + theme(text = element_text(color = "gray20"),
               legend.position = c("right"), # position the legend in the upper left 
               legend.justification = 0, # anchor point for legend.position.
               legend.text = element_text(size = 9, color = "gray10"),
               title = element_text(size = 15, face = "bold", color = "gray10"),
               axis.text = element_text(face = "bold"),
               panel.grid.major.y = element_blank(),
               panel.grid.major.x = element_blank()
)
p <- p + theme(axis.text.x=element_text(angle = 90, hjust = 1 , vjust = 0.5, 
                                        margin=margin(-8,0,-15,0)),
               axis.text.y=element_text(hjust = 1, vjust = 0.5, 
                                        margin = margin(0,-10,0,0)))
p
#
#
# Success. We have a very nice-looking seatmap chart summarizing the type of mutations
# That were identified in 25 genes while analyzing 10 samples. Publication-grade figure!
# For questions, contact me at damiano.fantini@gmail.com
# Or visit: http://www.biotechworld.it/bioinf/
#
#


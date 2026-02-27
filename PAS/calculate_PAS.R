##load required libraries:
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(dplyr)

##define data paths:
network_metadata_path = 'data/asynpn-meta.txt'
expr_path = 'data/expr.txt' #col1:gene, col2-n:samples.

##load data:
df <- read.delim(network_metadata_path)
exprs <- read.delim(expr_path) #note: exprs should be normalised to SNCA. 
samples <- colnames(exprs[,-1])

##z-score normalization function:
z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

##PAS activity score function: 
calculate_pas <- function(data) {
  #Activity: direction of activity a gene has on aSyn aggregation.
  #Z_Score: z-score.
  data$gene_activity <- data$Activity * data$Z_Score
  M <- nrow(data) ##M: Size of network.
  activity_score <- (1 / sqrt(M)) * sum(data$gene_activity)
  return(activity_score)
}

##run:
ratios_z <- apply(exprs[,-1], 2, z_score)
row.names(ratios_z) <- exprs$gene

results = data.frame(Sample=samples)
for (i in 1:nrow(results)) {
  temp <- data.frame(Z_Score = ratios_z[,results$Sample[i]])
  temp$Gene <- rownames(temp)
  temp <- left_join(df, temp, by = "Gene")
  results$PAS[i] <- calculate_pas(data=temp)
}

##normalize PAS:
results$PAS_Norm <- z_score(results$PAS)

##save results:
write.table(results, 'results/PAS_scores.txt', sep = "\t", row.names = FALSE)
message("PAS calculation complete. Results saved to results/PAS_scores.txt")

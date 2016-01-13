# Microbial Assemblage Normalized Transcript Analysis Plotting Function
# Based on manta (Bioconductor)
#
# Author: Robert Lampe
#
# References:
# - Marchetti A, Schruth DM, et al. (2012). Comparative metatran-
#   scriptomics identifies molecular bases for the physiological 
#   responses of phytoplankton to varying iron availability. 
#   Proceedings of the National Academy of Sciences 109, no.7, E317
# 
# Required Inputs:
# - counts - data frame with columns for counts, KO, and genus fields
# - exactTest - data frame of exactTest output for KO and counts - KO, logFC, logCPM, and PValue required
# - group1 - a vector for names of counts columns for the first group of comparison
# - group2 - a vector for names counts columns for the second group of comparison
#
# Dependencies:
# - caroline 0.7.6
# - plotrix 3.6-1
#
# FIXME:
# - Support replicates
# - Support custom colors

library(caroline)
library(plotrix)

mantaPlot <- function(counts, exactTest, group1, group2, taxa = 4,
                      title=NULL, xlab = "Average LogCPM", ylab = "Log Fold Change", ylim=c(-15, 15)) {

  # Check if replicates
  if (length(group1) > 1 || length(group2) > 1) {
 
  } else {
    counts <- subset(counts, select = c('KO', 'genus', group1, group2))
    counts[,group1] <- counts[,group1] / (sum(counts[,group1]) * 1e-6)
    counts[,group2] <- counts[,group2] / (sum(counts[,group2]) * 1e-6)
    counts$sum <- counts[,group1] + counts[,group2]   
  }
  
  # Build new data frame
  # Rows: KOs
  # Columns: All genera
  kos <- unique(counts$KO)
  genera <- as.character(unique(counts$genus))
  
  df <- as.data.frame(matrix(NA, nrow=length(kos), ncol=length(genera)))
  names(df) <- genera
  rownames(df) <- kos
  df <- df[order(rownames(df)),]
  
  # Calculate the percentage of each KO for each tax group
  total_kos <- groupBy(counts, by="KO", clmns=c('sum'), aggregation='sum')
  for (genus in genera) {
    df[,genus] <- (groupBy(counts[counts$genus == genus,], by="KO", clmns = c('sum'), aggregation='sum') / total_kos)
  }
  df[is.na(df)] <- 0
  
  # Simplify the dataframe
  tax_groups <- rev(names(tail(sort(colSums(df)), taxa)))
  df <- df[,tax_groups]
  df[,"Other"] <- 1 - rowSums(df)
  df <- df + 0.0000000001 ### Floating pies can't use 0
  
  # Build plot
  plot(exactTest$logCPM, exactTest$logFC, 
       cex = 0.1, 
       main = title, 
       xlab = xlab, ylab = ylab, 
       ylim=ylim)
  abline(h = 0)
  colors <- c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#7F7F7F")
  legend('topright', inset = 0.02, names(df), fill = colors, cex = 0.75)
  
  kos <- unique(exactTest$KO) # In case your ET output has less KOs than counts file for some reason
  for (ko in kos) {
    x = exactTest[exactTest$KO == ko,]$logCPM
    y = exactTest[exactTest$KO == ko,]$logFC
    pie_percents <- as.numeric(df[ko,])
    radius <- (x + abs(y))/ 100
    colors <- c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#7F7F7F")
    border = "black"
    if (exactTest[exactTest$KO == ko,]$PValue > 0.05) {
      colors <- paste(colors, "33", sep='')
      border = "#4D4D4D"
    }
    floating.pie(x, y, pie_percents, radius = radius, col = colors, border = border)
  }
}
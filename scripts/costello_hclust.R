# Cut R's hclust at k = 2:50 for a condensed distance vector, the way the Castro-Insua et al. 2018 script does it.
# Usage: Rscript costello_hclust.R <dir>, reading dist.txt and writing labels_ward.D2.txt and labels_average.txt.
dir <- commandArgs(TRUE)[1]
d <- scan(file.path(dir, "dist.txt"), quiet = TRUE)
n <- as.integer((1 + sqrt(1 + 8 * length(d))) / 2)
attr(d, "Size") <- n
attr(d, "Diag") <- FALSE
attr(d, "Upper") <- FALSE
attr(d, "method") <- "custom"
class(d) <- "dist"
for (m in c("ward.D2", "average")) {
  write.table(cutree(hclust(d, method = m), k = 2:50), file.path(dir, paste0("labels_", m, ".txt")), row.names = FALSE, col.names = FALSE)
}

library(ggplot2)
library(dplyr)

counts <-
  read.delim(
    "https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/thursday/all_counts.txt"
  )
counts_mat <- as.data.frame(counts[, 1])
colnames(counts_mat) <- c("Count")

counts_mat2 <-
  counts_mat %>% 
  filter(Count > 100) %>% 
  arrange(desc(Count))

counts_mat2$Genes <- 
  seq(1, nrow(counts_mat2))

counts_mat2 %>% 
  ggplot(aes(y = Count)) + 
  geom_histogram(fill = "red") + 
  theme_minimal()

counts_mat2$LogCounts <- 
  log2(counts_mat2$Count + 0.5)

counts_mat2 %>% 
  ggplot(aes(x = LogCounts)) + 
  geom_histogram(fill = "Purple", 
                 col ="black") + 
  theme_minimal()

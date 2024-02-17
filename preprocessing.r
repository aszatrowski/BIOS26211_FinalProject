library(tidyverse)

gse_20680 <- read.table("GEOData/GSE20680_no-metadata.txt", sep = "\t", header = TRUE, quote = "", comment.char = "!", skip = 0, fill = TRUE, nrow = 45300)

# remove auto-added X. pattern
colnames(gse_20680) <- gsub(pattern = "X.", replacement = "", x = colnames(gse_20680))
colnames(gse_20680) <- gsub(pattern = "\\.", replacement = "", x = colnames(gse_20680))

# Transpose matrix into tidy format: rows with observations (patients) and columns with gene IDs: doesn't quite work

# gse_20680_tidy <- as_tibble(t(gse_20680[1:45015,2:196]))
# colnames(gse_20680_tidy) <- paste0("g", colnames(gse_20680_tidy))
# gse_20680_tidy <- gse_20680_tidy %>%
#     mutate(sample_id = colnames(gse_20680)) %>%
#     select(sample_id, g1:g45015)

# columns: REF_ID for numeric gene index, then 195 columns for individuals
write_csv(gse_20680, "GEOData/GSE20680.csv")
library(tidyverse)
path = "/Users/austinszatrowski/Documents/UChicago/2024.1 Winter/BIOS 26211 Math Methods for Biological Sciences II/GEOData/"
gse_20680 <- read.table(paste0(path, "GSE20680_no-metadata.txt"),
    sep = "\t",
    header = TRUE,
    quote = "\"",
    comment.char = "!",
    skip = 0,
    fill = TRUE,
    nrow = 45300)

gse_20680_labels <- read.table(paste0(path, "GSE20680_labels.txt"),
    sep = "\t",
    header = FALSE,
    #quote = "\"",
    skip = 0,
    fill = TRUE,
    col.names = "cad_status")

gse_20680_labels_vec <- gsub(
    pattern = "disease state: ",
    replacement = "",
    gse_20680_labels$cad_status
)

# remove auto-added X. pattern
colnames(gse_20680) <- gsub(pattern = "X.", replacement = "", x = colnames(gse_20680))
colnames(gse_20680) <- gsub(pattern = "\\.", replacement = "", x = colnames(gse_20680))


# Transpose matrix into tidy format: rows with observations (patients) and columns with gene IDs: doesn't quite work

gse_20680_tidy <- as_tibble(t(gse_20680[1:45015,2:196]))
colnames(gse_20680_tidy) <- paste0("g", colnames(gse_20680_tidy))
gse_20680_tidy <- gse_20680_tidy %>%
    mutate(sample_id = colnames(gse_20680)[-1], cad_status = gse_20680_labels_vec) %>%
    select(sample_id, cad_status, g1:g45015)

# columns: REF_ID for numeric gene index, then 195 columns for individuals, then cad_status for CAD status
write_csv(gse_20680_tidy, "../GEOData/GSE20680.csv")
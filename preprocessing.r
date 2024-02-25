library(tidyverse)
# everyone will need to update this path since we have local copies of the data, or update to a remote folder 
path = "/Users/austinszatrowski/Documents/UChicago/2024.1 Winter/BIOS 26211 Math Methods for Biological Sciences II/GEOData/"

preprocess_matrix <- function(path, accession) {
    # read in gene expression matrix
    expression_matrix <- read.table(paste0(path, accession, "_ExpMatrix.txt"),
        sep = "\t",
        header = TRUE,
        quote = "\"",
        comment.char = "!",
        skip = 0,
        fill = TRUE,
        nrow = 45300)

    # read in labels vector from manual preprocessing
    labels <- read.table(paste0(path, accession, "_labels.txt"),
        sep = "\t",
        header = FALSE,
        #quote = "\"",
        skip = 0,
        fill = TRUE,
        col.names = "cad_status")

    # remove extraneous text from labels
    labels_vec <- gsub(
        pattern = "disease state: ",
        replacement = "",
        labels$cad_status
    )

    # remove auto-added X. pattern
    colnames(expression_matrix) <- gsub(pattern = "X.", replacement = "", x = colnames(expression_matrix))
    colnames(expression_matrix) <- gsub(pattern = "\\.", replacement = "", x = colnames(expression_matrix))


    # Transpose matrix into tidy format: rows with observations (patients) and columns with gene IDs + CAD status

    # Columns: sample_id human sample ID from GEO, cad_status for CAD status (0: control, 1: intermediate, 2: severe disease) then 45,017 gene expression columns

    # expression_matrix_tidy <- as_tibble(t(expression_matrix[1:45015,2:196]))
    # colnames(expression_matrix_tidy) <- paste0("g", colnames(expression_matrix_tidy))
    # expression_matrix_tidy <- expression_matrix_tidy %>%
    #     mutate(sample_id = colnames(expression_matrix)[-1], cad_status = labels_vec) %>%
    #     select(sample_id, cad_status, g1:g45015)

    write_csv(expression_matrix, paste0("../GEOData/", accession, "_RBEReady.csv"))
}


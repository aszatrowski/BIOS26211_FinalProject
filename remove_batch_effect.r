library(tidyverse)
library(limma)
library(sva)

path = "../GEOData/"

gse20680 <- read_csv(paste0(path, "GSE20680_RBEReady.csv"))
gse20681 <- read_csv(paste0(path, "GSE20681_RBEReady.csv"))

df_unadjusted <- full_join(gse20680, gse20681, by = "ID_REF")

batches = c(
    "batch",
    rep("gse20680", ncol(gse20680)-1),
    rep("gse20681", ncol(gse20681)-1)
)

gse20680_labels <- scan(paste0(path, "GSE20680_labels.txt"),
                        what = "character",
                        sep = "\t")

gse20681_labels <- scan(paste0(path, "GSE20681_labels.txt"),
                        what = "character",
                        sep = "\t")

labels = c("cad_status", gse20680_labels, gse20681_labels)

df_unadjusted_covars <- rbind(df_unadjusted, labels, batches)

#df_adjusted <- removeBatchEffect(df_unadjusted[,1:394], batch=batches)
df_adj <- ComBat(df_unadjusted, batches, mod=modBatch, par.prior = TRUE, prior.plots = FALSE)
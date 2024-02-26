library(tidyverse)
library(limma)

path = "../GEOData/"

gse20680 <- read_csv(paste0(path, "GSE20680_RBEReady.csv"))
gse20681 <- read_csv(paste0(path, "GSE20681_RBEReady.csv"))

df_unadjusted <- full_join(gse20680, gse20681, by = "ID_REF")

batches = c(
    rep("gse20680", ncol(gse20680)),
    rep("gse20681", ncol(gse20681))
)

df_adjusted <- removeBatchEffect(df_unadjusted[,1:394], batch=batches)
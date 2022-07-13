update.packages(oldPkgs = c("withr", "rlang"))

if (!requireNamespace('remotes', quietly = TRUE)) {
  install.packages('remotes')
}
remotes::install_github('satijalab/azimuth', ref = 'master')

library(Seurat)
library(SeuratData)

SeuratData::AvailableData()

# Install the PBMC systematic comparative analyis (pmbcsca) dataset
InstallData("pbmc3k")

# returns a Seurat object named pbmcsca
LoadData("pbmc3k")

library(Azimuth)
saveRDS(pbmc3k,file="pbmc3k.rds")

#Launch Azimuth Locally
Azimuth::AzimuthApp()

#See All availible data in Azimuth
available_data <- AvailableData()
available_data[grep("Azimuth", available_data[, 3]), 1:3]



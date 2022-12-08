# Last checked: 7th December 2022
# R version 3.6.3

# Launch R
module load R-base/3.6.3
module load R-cbrg/202006
R

# Load packages
library(dplyr) # version 1.0.0
library(Seurat) # version 2.0.1
library(Matrix) # version 1.2-18
library(ggplot2) # version 3.3.6
library(plyr) # version 1.8.6

# Read Seurat object
path <- "/datashare/wwen/Eleni_JEM_2021/data/"
file <- "All.data.Seurat.Robj"
load(paste(path, file, sep=""))

# Plot tSNE: CD96
path <- "/project/meadlab/wwen/Eleni/Eleni_2020/Github/Elouka_JEM_2021/figures/"
file <- "figure_6B.png"
png(paste(path, file, sep=""), width=450, height=450)

FeaturePlot(all.data,
            c("CD96"),
            cols.use=c("grey","cyan","green","yellow","red"),
            pt.size=0.9,
            no.axes=TRUE) + theme_set(theme_cowplot())

dev.off()

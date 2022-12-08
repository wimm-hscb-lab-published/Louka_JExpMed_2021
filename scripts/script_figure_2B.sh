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
library(ggplot2) # version 3.3.1
library(plyr) # version 1.8.6

# Read Seurat object
path <- "/datashare/wwen/Eleni_JEM_2021/data/"
file <- "All.data.Seurat.Robj"
load(paste(path, file, sep=""))

# Convert cluster IDs from 0- to 1-based
cls <- all.data@ident
cls <- as.numeric(as.character(cls)) + 1
cls <- as.factor(as.character(cls))
names(cls) <- names(all.data@ident)
all.data@ident <- cls

# Plot tSNE: Cluster IDs
path <- "/project/meadlab/wwen/Eleni/Eleni_2020/Github/Elouka_JEM_2021/figures/"
file <- "figure_2B.png"
png(paste(path, file, sep=""), width=500, height=500)

TSNEPlot(all.data, pt.size=0.9)

dev.off()

# Save coordinates
    # Retrieve coordinates
    df.coord <- as.data.frame(all.data@dr[["tsne"]]@cell.embeddings)
    . <- data.frame("cell.id"=row.names(df.coord))
    df.coord <- cbind.data.frame(., df.coord)
    
    # Annotate cluster ID
    df.cl <- data.frame(all.data@ident)
    . <- data.frame("cell.id"=row.names(df.cl))
    df.cl <- cbind.data.frame(., df.cl)
    df.coord <- join(df.coord, df.cl, by="cell.id", type="left")

    # Save file
    #path <- "/project/meadlab/wwen/Eleni/Eleni_2020/Github/Elouka_JEM_2021/figures/"
    #file <- "figure_2B_data.txt"
    #write.table(df.coord, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

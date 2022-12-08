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

# Retrieve gene expression values
. <- all.data@data["CD96",]
. <- data.frame("Cell.Barcode"=names(.), "Gene"=., stringsAsFactors=FALSE)

# Retrieve patient ID
barcode.split <- strsplit(as.character(.$Cell.Barcode), split="\\_")
.$Patient.ID <- sapply(barcode.split, function(x) {x[1]})
.$Patient.ID <- factor(.$Patient.ID, levels=c("CB35", "CB37", "HR", "AK"))
table(.$Patient.ID)

# Definition
data <- .
x <- data$Patient.ID
y <- data$Gene
z <- data$Patient.ID
maintitle <- ""
xtitle <- ""
ytitle <- "Log UMI counts"

# Plot
plot <- ggplot() +
        geom_violin(data, mapping=aes(x=x, y=y, fill=z), scale="width", adjust=3) +
        scale_fill_manual(values=c("blue", "deepskyblue", "firebrick", "orange")) +
        labs(title=maintitle, x=xtitle, y=ytitle) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_blank(),
              axis.line = element_line(colour="black"),
              axis.text.x=element_text(colour="black"),
              axis.text=element_text(size=12),
              axis.title=element_text(size=16),
              plot.title=element_text(size=16),
              legend.position="none"
              )

# Save plot
path <- "/project/meadlab/wwen/Eleni/Eleni_2020/Github/Elouka_JEM_2021/figures/"
file <- "figure_6C.png"
ggsave(paste(path, file, sep=""), plot, width=6, height=6)
            

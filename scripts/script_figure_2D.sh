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
        
# Determine cell type frequency in each clusters
    # Extract clusters
    . <- all.data@ident
    . <- data.frame("Cell.Barcode"=names(.), "Cluster"=.)
    
    # Extract cell type
    barcode.split <- strsplit(as.character(.$Cell.Barcode), split="\\_")
    .$Cell.Type <- sapply(barcode.split, function(x) {x[1]})
    .$Cell.Type[which(.$Cell.Type=="HR")] <- "JMML"
    .$Cell.Type[which(.$Cell.Type=="AK")] <- "JMML"
    .$Cell.Type[which(.$Cell.Type=="CB35")] <- "CB"
    .$Cell.Type[which(.$Cell.Type=="CB37")] <- "CB"
    .$Cell.Type <- factor(.$Cell.Type, levels=c("CB", "JMML"))
    table(.$Cell.Type)
    
    # Reformat data frame
    . <- as.data.frame(table(.$Cluster, .$Cell.Type))
    names(.) <- c("Cluster", "Cell.Type", "Freq")
    
    # Write file
    #path <- "/project/meadlab/wwen/Eleni/Eleni_2020/Github/Elouka_JEM_2021/figures/"
    #file <- "figure_2C_data.txt"
    #write.table(., paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Stacked barplots
    # Definition
    data <- .
    x <- data$Cluster
    y <- data$Freq
    z <- data$Cell.Type
    maintitle <- ""
    ytitle <- "Cells"
    xtitle <- "Clusters"
    legendtitle <- "Cell Type"
    fivenum(y); ymin <- 0; ymax <- 5500; yinterval <- 500
    
    # Plot
    plot <- ggplot() +
            geom_bar(data, mapping=aes(x=x, y=y, fill=z), color="black", position="stack", stat="identity") +
            scale_y_continuous(limits=c(ymin, ymax), breaks=seq(ymin, ymax, by=yinterval)) +
            scale_fill_manual(values=c("darkblue", "firebrick")) +
            labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.background=element_blank(),
                  axis.line=element_line(colour="black"),
                  axis.text=element_text(size=13),
                  axis.title=element_text(size=18),
                  axis.text.x=element_text(colour="black"),
                  plot.title=element_text(size=18)
                  )
                  
    # Save plot
    path <- "/project/meadlab/wwen/Eleni/Eleni_2020/Github/Elouka_JEM_2021/figures/"
    file <- "figure_2C.png"
    ggsave(paste(path, file, sep=""), plot, width=8, height=6)

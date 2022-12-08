# Last checked: 7th December 2022
# R version 4.0.5

# Load packages
library(gplots) # version 3.1.1

# Read data object
path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/data/"
file <- "data.rdata"
load(paste(path, file, sep=""))

# Retrieve bulk genotyping data
df <- object$bulk_genotyping

# Tidy patient IDs
names(df) <- paste("ID", names(df), sep="")

# Convert data frame to matrix
df <- as.matrix(df)

# Heatmap
    # Definitions
    data <- df
    colors <- c("gray96",'khaki4','maroon','firebrick','darkorange1','orange','green4','darkorchid', "olivedrab3",'blue','blue','blue','darkorchid','dodgerblue','grey70')

    # Plot
    path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
    file <- "figure_3A.pdf"
    pdf(paste(path, file, sep=""), width=8, height=4)
    
    heatmap.2(data,
              col=colors,
              dendrogram="none",
              Rowv=FALSE,
              Colv=FALSE,
              trace="none",
              key=TRUE,
              density.info="none",
              colsep=c(0:ncol(data)),
              rowsep=c(0:nrow(data)),
              sepcolor="gray0",
              sepwidth=c(0.005,0.005),
              srtRow=0,
              offsetRow=-60,
              srtCol=90,
              offsetCol=-25,
              cexRow=1.2)
    
    dev.off()
    

# Last checked: 7th December 2022
# R version 4.0.5

# Load packages
library(ggplot2) # version 3.3.6
library(ggExtra) # version 0.10.0

# Read data object
path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/data/"
file <- "data.rdata"
load(paste(path, file, sep=""))

# Read file
df <- object$facs

# Scatterplot
    # Definition
    data <- df
    x <- data$HbF
    y <- data$CD96
    maintitle <- ""
    ytitle <- "log2(TPM + 1)"
    xtitle <- ""
    ymin <- 0; ymax <- 150
    xmin <- 0; xmax <- 65
    
    # Plot
    plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y), size=5, color='darkred') +
            geom_smooth(data, mapping=aes(x=x, y=y), method=lm, linetype='dashed', color='blue', fill='grey') +
            ylim(ymin, ymax) +
            xlim(xmin, xmax) +
            theme_bw() +
            theme(axis.title=element_blank())
            
    # Save
    path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
    file <- "figure_6I.png"
    ggsave(paste(path, file, sep=""), plot, width=7, height=5)

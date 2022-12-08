# Last checked: 7th December 2022
# R version 4.0.5

# Load packages
library(ggplot2) # version 3.3.6

# Read data object
path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/data/"
file <- "data.rdata"
load(paste(path, file, sep=""))

# Read file
df <- object$bulk_rnaseq$gsea
    
# Subset significant gene sets
df <- df[which(df$NOM.p.val < 0.10), ]

# Transform p-values
    # Transform
    df$NOM.p.val.trans <- -log10(df$NOM.p.val)
    
    # Set upper limit
    df$NOM.p.val.trans[is.infinite(df$NOM.p.val.trans)] <- 6.0

# Set factor levels
df <- df[order(df$NOM.p.val.trans),]
df <- df[c(nrow(df):1), ]
df$NAME <- gsub("^HALLMARK_", "", df$NAME)
df$NAME <- factor(df$NAME, levels=df$NAME)

# Barplot
    # Definition
    data <- df
    x <- data$NAME
    y <- data$NOM.p.val.trans
    z <- round(data$NES, digits=2)
    maintitle <- ""
    ytitle <- "-log10(Nominal p-value)"
    xtitle <- ""
    yintercept <- -log10(0.10)
    ymin <- 0
    ymax <- max(y) + 1
    yinterval <- 2
    
    # Plot
    plot <- ggplot() +
        geom_bar(data, mapping=aes(x=x, y=y), fill="blue", stat="identity") +
        geom_text(data, mapping=aes(x=x, y=y, label=z), vjust=0.5, hjust=-0.2, size=5) +
        geom_hline(yintercept=yintercept, linetype="dashed", size=1.5) +
        scale_y_continuous(limits=c(ymin, ymax), breaks=seq(ymin, ymax, by=yinterval)) +
        labs(title=maintitle, x=xtitle, y=ytitle) +
        theme_bw() +
        theme(panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.minor=element_blank(),
              panel.grid.major=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line=element_line(colour="black"),
              axis.text.x=element_text(size=16, colour="black"),
              axis.text.y=element_text(size=12, colour="black"),
              axis.text=element_text(size=16),
              axis.title=element_text(size=16),
              plot.title=element_text(size=16),
              legend.position="none"
              ) +
        coord_flip()
              
    # Save file
    path <- "/Users/seanwen/Documents/Eleni/Eleni_2020/Github/Louka_JExpMed_2021/figures/"
    file <- "figure_5I.png"
    ggsave(paste(path, file, sep=""), plot, width=10, height=10)
    

